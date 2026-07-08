#!/usr/bin/env python3
"""
Utility to replace aromatic rings in a GROMACS trajectory with their centres of mass.

This script reads a processed GROMACS topology (``.top``) alongside the binary
topology (``.tpr``) and a coordinate trajectory (``.xtc``). It detects 5– and
6-membered aromatic rings by analysing the connectivity defined in the topology.
For each frame of the trajectory it computes the centre of mass (COM) of the
heavy (non-hydrogen) atoms forming each ring.  The positions of all heavy
atoms belonging to a ring are then replaced with a single particle located
at the COM.  All other atoms (hydrogens, atoms outside rings, molecules
without rings) are left untouched.  The modified coordinates are written to
a new trajectory file (``out.xtc``) and, optionally, a corresponding
``.gro`` structure.

The COM calculation properly accounts for periodic boundary conditions (PBC):
each ring is unwrapped so that all atoms are in the same periodic image
before computing its centre of mass.  This avoids artefacts such as
unphysically small COM–atom distances in radial distribution functions.

Usage::

    python com_aromatic_rings.py topol.tpr topol.top traj.xtc out.xtc

The output structure will additionally be written as ``out.gro``.

"""

import sys
import re
import numpy as np
import networkx as nx
import MDAnalysis as mda

# Optional import: ``apply_PBC`` is only available in newer MDAnalysis releases.
try:
    from MDAnalysis.lib.distances import apply_PBC  # type: ignore
except Exception:
    apply_PBC = None  # fallback defined in compute function


def parse_topology(topfile: str):
    """Parse a processed GROMACS ``.top`` file.

    Parameters
    ----------
    topfile : str
        Path to the processed topology file.

    Returns
    -------
    tuple
        ``(moltypes, mol_order)`` where ``moltypes`` is a dictionary keyed by
        molecule type name containing keys ``natoms``, ``atoms`` and ``bonds``.
        ``mol_order`` is a list of ``(molname, count)`` pairs describing the
        order and multiplicity of molecule types in the ``[ molecules ]`` section.

    Notes
    -----
    This parser is intentionally simple: it only extracts atom names and bonds
    and ignores masses and charges because the binary topology (``.tpr``)
    contains these.  Comments starting with ``;`` and blank lines are
    stripped.  If a molecule type appears more than once in ``[ molecules ]``
    the bonds and atom definitions are reused for each instance.
    """
    moltypes = {}
    mol_order = []

    current_section = None
    current_moltype = None

    with open(topfile, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith(';'):
                continue

            # section headers
            if line.startswith('[') and line.endswith(']'):
                current_section = line.strip('[]').strip().lower()
                if current_section == 'moleculetype':
                    current_moltype = None
                continue

            if current_section == 'moleculetype':
                # Expect ``name  nrexcl`` on first non-comment line
                if current_moltype is None:
                    parts = line.split()
                    if parts:
                        molname = parts[0]
                        current_moltype = molname
                        moltypes[current_moltype] = {
                            'natoms': 0,
                            'atoms': [],
                            'bonds': []
                        }
                continue

            if current_section == 'atoms':
                parts = line.split()
                if len(parts) < 5:
                    continue
                local_id = int(parts[0])
                resnr = parts[2]
                resname = parts[3]
                atomname = parts[4]
                if current_moltype is None:
                    raise ValueError("Encountered [ atoms ] without a current moleculetype")
                moltypes[current_moltype]['atoms'].append({
                    'id': local_id,
                    'name': atomname,
                    'resname': resname,
                    'resnr': resnr,
                })
                # track highest index for natoms
                moltypes[current_moltype]['natoms'] = max(
                    moltypes[current_moltype]['natoms'], local_id
                )
                continue

            if current_section == 'bonds':
                parts = line.split()
                if len(parts) < 2:
                    continue
                try:
                    i = int(parts[0])
                    j = int(parts[1])
                except ValueError:
                    continue
                if current_moltype is None:
                    raise ValueError("Encountered [ bonds ] without a current moleculetype")
                moltypes[current_moltype]['bonds'].append((i, j))
                continue

            if current_section == 'molecules':
                parts = line.split()
                if len(parts) < 2:
                    continue
                molname = parts[0]
                try:
                    count = int(parts[1])
                except ValueError:
                    continue
                mol_order.append((molname, count))
                continue

    return moltypes, mol_order


def parse_topology_full(topfile: str):
    """Parse a standalone processed GROMACS .top preserving atom fields needed to rewrite [atoms] & [bonds].

    Returns
    -------
    tuple
        (preamble_lines, moltypes, mol_order)

    Notes
    -----
    - Keeps a copy of all lines before the first [ moleculetype ] as preamble.
    - For each moleculetype, stores raw atom fields so we can rewrite without losing columns.
    - Only parses [ moleculetype ], [ atoms ], [ bonds ], [ molecules ] but retains
      all other blocks as raw text inside each moleculetype so we can drop angles/dihedrals cleanly.
    """
    preamble = []
    moltypes = {}
    mol_order = []

    current_section = None
    current_moltype = None
    seen_first_moleculetype = False

    # Store per-moltype raw blocks so we can remove sections
    mol_blocks = {}  # molname -> list of (section_name, lines)
    section_lines = []

    def flush_section():
        nonlocal section_lines, current_section, current_moltype
        if current_moltype is None or current_section is None:
            section_lines = []
            return
        mol_blocks.setdefault(current_moltype, []).append((current_section, section_lines))
        section_lines = []

    with open(topfile, "r") as f:
        for raw in f:
            line = raw.rstrip("\n")

            # Detect section header
            m = re.match(r"^\s*\[\s*([A-Za-z0-9_]+)\s*\]\s*(;.*)?$", line)
            if m:
                # flush previous section if we were collecting
                flush_section()

                current_section = m.group(1).strip().lower()

                if current_section == "moleculetype":
                    seen_first_moleculetype = True
                    current_moltype = None  # next non-comment line defines it
                else:
                    # keep collecting lines for sections belonging to current moltype
                    pass

                # Always keep header lines as part of section_lines so we can reproduce formatting if needed
                section_lines = [line]
                continue

            # Before first moleculetype: keep as preamble verbatim
            if not seen_first_moleculetype:
                preamble.append(line)
                continue

            # After first moleculetype: we parse content
            # If we are inside [ moleculetype ] and current_moltype not set, parse it
            if current_section == "moleculetype" and current_moltype is None:
                stripped = line.strip()
                section_lines.append(line)
                if not stripped or stripped.startswith(";"):
                    continue
                parts = stripped.split()
                if len(parts) >= 1:
                    molname = parts[0]
                    current_moltype = molname
                    moltypes[current_moltype] = {
                        "nrexcl_line": line,      # keep raw line for later reuse
                        "atoms": [],              # list of dicts with raw fields
                        "bonds": [],              # list of (i,j,funct,rest_tokens,raw_comment)
                    }
                continue

            # Collect section lines for later rewrite/removal
            section_lines.append(line)

            # Parse atoms/bonds/molecules into structured form
            stripped = line.strip()
            if not stripped or stripped.startswith(";"):
                continue

            if current_section == "atoms":
                # Typical: nr type resnr resid atom cgnr charge mass [..]
                parts = stripped.split()
                if len(parts) < 8:
                    continue
                nr = int(parts[0])
                atomtype = parts[1]
                resnr = parts[2]
                resid = parts[3]
                atomname = parts[4]
                cgnr = parts[5]
                charge = parts[6]
                mass = parts[7]
                extra = parts[8:] if len(parts) > 8 else []
                moltypes[current_moltype]["atoms"].append({
                    "nr": nr,
                    "type": atomtype,
                    "resnr": resnr,
                    "resid": resid,
                    "atom": atomname,
                    "cgnr": cgnr,
                    "charge": charge,
                    "mass": mass,
                    "extra": extra,
                    "raw": line,
                })
                continue

            if current_section == "bonds":
                # Typical: i j funct [..]
                parts = stripped.split()
                if len(parts) < 3:
                    continue
                try:
                    i = int(parts[0])
                    j = int(parts[1])
                    funct = parts[2]
                except ValueError:
                    continue
                rest = parts[3:] if len(parts) > 3 else []
                moltypes[current_moltype]["bonds"].append((i, j, funct, rest, line))
                continue

            if current_section == "molecules":
                parts = stripped.split()
                if len(parts) >= 2:
                    name = parts[0]
                    try:
                        count = int(parts[1])
                    except ValueError:
                        continue
                    mol_order.append((name, count))
                continue

    # flush last
    flush_section()

    # attach raw blocks
    for molname, blocks in mol_blocks.items():
        if molname in moltypes:
            moltypes[molname]["raw_blocks"] = blocks
        else:
            moltypes[molname] = {"atoms": [], "bonds": [], "raw_blocks": blocks}

    return preamble, moltypes, mol_order


def _is_hydrogen_name(atomname: str) -> bool:
    return atomname.strip().upper().startswith("H")


def _is_carbon_name(atomname: str) -> bool:
    return atomname.strip().upper().startswith("C")


def rewrite_moleculetype_with_com_and_hetero(mol: dict, rings_local_1b: list[list[int]]):
    """Rewrite one moleculetype: remove ring heavy atoms except COM rep + heteroatoms; rewrite bonds as type 5.

    Parameters
    ----------
    mol : dict
        One entry from moltypes (must contain 'atoms' and 'bonds').
    rings_local_1b : list of rings
        Each ring is a list of 1-based atom indices (local to moleculetype) of HEAVY atoms.

    Returns
    -------
    tuple
        (new_atoms, new_bonds, old_to_new, com_reps_1b, hetero_1b_per_ring)

    Notes
    -----
    - COM representative for a ring = first atom in ring with 1 < mass < 13 (in the ring ordering). 
    - Heteroatoms in ring = atoms in ring with mass > 13. 
    - Atoms removed = ring atoms except (rep + heteroatoms).
    - Bonds:
        * every bond funct set to '5'
        * if a bond endpoint is removed ring atom, it is redirected to the ring rep
        * bonds internal to removed atoms collapse accordingly (duplicates removed)
        * if a heteroatom is in the ring, add fictitious bond hetero–rep (type 5)
    """
    atoms = mol["atoms"]
    bonds = mol["bonds"]

    # Map local index -> atom record
    atom_by_nr = {a["nr"]: a for a in atoms}

    # Determine which atoms belong to rings, reps, heteros
    ring_atoms = set()
    rep_atoms = set()
    hetero_atoms = set()
    hetero_per_ring = []

    for ring in rings_local_1b:
        for i in ring:
            ring_atoms.add(i)

        # choose rep = first carbon atom in ring (lowest nr with name starting C)
        rep = None
        for i in ring:
            try:
                m = float(atom_by_nr[i]["mass"])
            except Exception:
                m = 0.0
            if (m > 1.0) and (m < 13.0):
                rep = i
                break
        if rep is None:
            rep = min(ring)
        rep_atoms.add(rep)
        
        # heteroatoms = atoms in ring with mass > 13 (excluding representative)
        heteros = []
        for i in ring:
            if i == rep:
                continue
            try:
                m = float(atom_by_nr[i]["mass"])
            except Exception:
                m = 0.0
            if m > 13.0:
                heteros.append(i)

        for h in heteros:
            hetero_atoms.add(h)
        hetero_per_ring.append((rep, heteros, ring))

    # Remove ring atoms except reps + heteros
    remove_atoms = set()
    for i in ring_atoms:
        if i in rep_atoms:
            continue
        if i in hetero_atoms:
            continue
        remove_atoms.add(i)

    # Build old->new map (renumber local indices)
    keep_atoms = [a for a in atoms if a["nr"] not in remove_atoms]
    keep_atoms.sort(key=lambda x: x["nr"])

    old_to_new = {}
    new_atoms = []
    for new_nr, a in enumerate(keep_atoms, start=1):
        old_to_new[a["nr"]] = new_nr
        # copy record; nr updated
        rec = dict(a)
        rec["nr"] = new_nr
        new_atoms.append(rec)

    # For each ring rep: make it represent the COM.
    # We keep its atom name/type/resid etc, BUT update charge+mass to be sums of the ring heavy atoms used for COM:
    # ring heavy atoms = ring list (already heavy)
    # (If you want to keep original charge/mass, remove the next block.)
    for ring in rings_local_1b:
        rep_old = None
        for i in ring:
            try:
                m = float(atom_by_nr[i]["mass"])
            except Exception:
                m = 0.0
            if (m > 1.0) and (m < 13.0):
                rep_old = i
                break
        if rep_old is None:
            rep_old = min(ring)
        rep_new = old_to_new[rep_old]


        # Sum over ring atoms (including heteroatoms) for COM pseudo-atom
        charges = []
        masses = []
        for i in ring:
            ai = atom_by_nr[i]
            try:
                charges.append(float(ai["charge"]))
            except Exception:
                charges.append(0.0)
            try:
                masses.append(float(ai["mass"]))
            except Exception:
                masses.append(0.0)

        qsum = sum(charges)
        msum = sum(masses)

        # Update rep record
        new_atoms[rep_new - 1]["charge"] = f"{qsum:.6f}"
        new_atoms[rep_new - 1]["mass"] = f"{msum:.6f}"

        # IMPORTANT: name of COM in .gro should be the representative atom kept. We therefore do not rename here.

    removed_to_rep = {}
    for rep_old, heteros, ring in hetero_per_ring:
        # Map only atoms that belong to THIS ring
        for i in ring:
            if i in remove_atoms:
               removed_to_rep[i] = min(removed_to_rep.get(i, rep_old), rep_old) 


    def map_endpoint(old_idx: int) -> int:
        """Map old atom index to new atom index, redirecting removed ring atoms to the rep."""
        if old_idx in remove_atoms:
            # redirect
            rep_old = removed_to_rep.get(old_idx)
            if rep_old is None:
                # Should not happen; drop
                return -1
            return old_to_new[rep_old]
        # keep as-is
        return old_to_new.get(old_idx, -1)

    # Rewrite bonds
    new_bond_set = set()
    new_bonds = []

    for (i, j, funct, rest, rawline) in bonds:
        ni = map_endpoint(i)
        nj = map_endpoint(j)
        if ni < 1 or nj < 1:
            continue
        if ni == nj:
            continue
        a, b = (ni, nj) if ni < nj else (nj, ni)
        if (a, b) in new_bond_set:
            continue
        new_bond_set.add((a, b))
        # force funct = 5
        new_bonds.append((a, b, "5", []))

    # Add fictitious hetero–COM (rep) bonds
    for rep_old, heteros, _ring in hetero_per_ring:
        rep_new = old_to_new[rep_old]
        for h_old in heteros:
            h_new = old_to_new[h_old]
            a, b = (rep_new, h_new) if rep_new < h_new else (h_new, rep_new)
            if a == b:
                continue
            if (a, b) in new_bond_set:
                continue
            new_bond_set.add((a, b))
            new_bonds.append((a, b, "5", []))

    # Sort bonds for tidy output
    new_bonds.sort(key=lambda x: (x[0], x[1]))

    return new_atoms, new_bonds, old_to_new


def write_modified_top(
    in_top: str,
    out_top: str,
    rings_per_moltype: dict,
    remove_sections: set[str] | None = None,
):
    """Write a new .top reflecting COM/hetero changes and bond redirection.

    Parameters
    ----------
    in_top : str
        Input standalone processed topology.
    out_top : str
        Output topology path.
    rings_per_moltype : dict
        {molname: [rings]} with rings as local 1-based indices (heavy atoms).
    remove_sections : set[str] or None
        Moleculetype-local sections to drop. By default drops angles/dihedrals.
        You can extend it e.g. {'angles','dihedrals','pairs','constraints','exclusions'}.
    """
    if remove_sections is None:
        remove_sections = {"angles", "dihedrals"}

    preamble, moltypes, mol_order = parse_topology_full(in_top)

    with open(out_top, "w") as out:
        # write preamble verbatim
        for line in preamble:
            out.write(line + "\n")

        # For each moleculetype block in the original, rewrite atoms/bonds and drop sections
        for molname, mol in moltypes.items():
            # If this moltype has no atoms, just skip
            if "atoms" not in mol or len(mol["atoms"]) == 0:
                # write raw blocks if present
                if "raw_blocks" in mol:
                    for sec, lines in mol["raw_blocks"]:
                        for l in lines:
                            out.write(l + "\n")
                continue

            rings_local_1b = rings_per_moltype.get(molname, [])
            # Rewrite only if rings exist; if none, keep atoms/bonds but still set bond funct=5 and drop sections
            if rings_local_1b:
                new_atoms, new_bonds, _ = rewrite_moleculetype_with_com_and_hetero(mol, rings_local_1b)
            else:
                # keep all atoms, but ensure bonds funct=5
                new_atoms = []
                for i, a in enumerate(sorted(mol["atoms"], key=lambda x: x["nr"]), start=1):
                    rec = dict(a)
                    rec["nr"] = i
                    new_atoms.append(rec)
                new_bonds = []
                bond_set = set()
                for (i, j, funct, rest, rawline) in mol["bonds"]:
                    a, b = (i, j) if i < j else (j, i)
                    if (a, b) in bond_set:
                        continue
                    bond_set.add((a, b))
                    new_bonds.append((a, b, "5", []))
                new_bonds.sort(key=lambda x: (x[0], x[1]))

            # Now write a clean moleculetype with only moleculetype/atoms/bonds plus any other kept blocks
            # We try to preserve the original [ moleculetype ] header and line.
            # Find the original nrexcl line (already stored as raw)
            out.write("\n[ moleculetype ]\n")
            # Use original nrexcl line if present, else default
            nrexcl_line = mol.get("nrexcl_line")
            if nrexcl_line and not nrexcl_line.strip().startswith("["):
                out.write(nrexcl_line.strip() + "\n")
            else:
                out.write(f"{molname} 3\n")

            # Write atoms
            out.write("\n[ atoms ]\n")
            out.write("; nr  type  resnr  resid  atom  cgnr  charge  mass\n")
            for a in new_atoms:
                extra = " " + " ".join(a["extra"]) if a.get("extra") else ""
                out.write(
                    f"{a['nr']:5d} {a['type']:>8s} {a['resnr']:>5s} {a['resid']:>6s} {a['atom']:>6s} "
                    f"{a['cgnr']:>5s} {a['charge']:>10s} {a['mass']:>10s}{extra}\n"
                )

            # Write bonds with funct=5
            out.write("\n[ bonds ]\n")
            out.write("; i  j  funct\n")
            for (i, j, funct, rest) in new_bonds:
                out.write(f"{i:5d} {j:5d}  5\n")

            # Optionally: keep other sections from original moleculetype except the removed ones
            # We drop whole blocks if their section name is in remove_sections.
            # IMPORTANT: many other blocks may reference removed atoms; extend remove_sections if needed.
            raw_blocks = mol.get("raw_blocks", [])
            for sec, lines in raw_blocks:
                secname = sec.lower()
                if secname in ("moleculetype", "atoms", "bonds", "molecules"):
                    continue
                if secname in remove_sections:
                    continue
                # write the block unchanged
                out.write("\n")
                for l in lines:
                    out.write(l + "\n")

        # Finally, write [ system ] and [ molecules ] from the original file:
        # They are present in raw_blocks of the last parsed state, but easier is:
        # just append a minimal [ system ] and [ molecules ] using mol_order.
        out.write("\n[ system ]\n")
        out.write("; generated by com_aromatic_rings.py\n")
        out.write("Reduced system\n")

        out.write("\n[ molecules ]\n")
        out.write("; name  count\n")
        for name, count in mol_order:
            out.write(f"{name:20s} {count:d}\n")


def find_rings_in_moltypes(moltypes: dict, min_size: int = 5, max_size: int = 6):
    """Identify cycles of heavy atoms (non-hydrogen) in each moleculetype.

    Graph theory is used to find simple cycles of size between ``min_size`` and
    ``max_size`` in the heavy-atom graph.  Each ring is returned as a sorted list
    of 1-based atom indices local to the molecule.

    Parameters
    ----------
    moltypes : dict
        Molecule definitions as returned by :func:`parse_topology`.
    min_size, max_size : int, optional
        Minimum and maximum ring sizes to consider.  Defaults to 5 and 6.

    Returns
    -------
    dict
        ``{molname: [rings]}`` mapping each molecule name to a list of rings
        (each a list of local atom indices).
    """
    rings_per_moltype = {}
    for molname, data in moltypes.items():
        atoms = data['atoms']
        bonds = data['bonds']
        # Build graph of heavy atoms
        G = nx.Graph()
        heavy_ids = set()
        for at in atoms:
            i = at['id']
            name = at['name'].strip().upper()
            if name.startswith('H'):
                continue
            heavy_ids.add(i)
            G.add_node(i)
        for i, j in bonds:
            if i in heavy_ids and j in heavy_ids:
                G.add_edge(i, j)
        if not heavy_ids:
            rings_per_moltype[molname] = []
            continue
        # Find cycles in the graph
        raw_cycles = nx.cycle_basis(G)
        rings = []
        for cycle in raw_cycles:
            size = len(cycle)
            if min_size <= size <= max_size:
                rings.append(sorted(cycle))
        # Remove duplicates (cycles can appear in different orientations)
        unique = []
        seen = set()
        for ring in rings:
            key = tuple(ring)
            if key not in seen:
                seen.add(key)
                unique.append(ring)
        rings_per_moltype[molname] = unique
    return rings_per_moltype


def expand_rings_to_global_indices(moltypes: dict, mol_order: list, rings_per_moltype: dict):
    """Expand local ring indices to global atom indices (0-based).

    Parameters
    ----------
    moltypes : dict
        Molecule definitions from :func:`parse_topology`.
    mol_order : list
        List of ``(molname, count)`` describing the order of molecules in the system.
    rings_per_moltype : dict
        Mapping of molecule type to lists of rings (1-based local indices).

    Returns
    -------
    list
        A list of rings, each defined as a list of global 0-based atom indices.
    """
    global_rings = []
    atom_offset = 0
    for molname, count in mol_order:
        if molname not in moltypes:
            raise ValueError(f"Unknown molecule type '{molname}' in [ molecules ] section")
        natoms = moltypes[molname]['natoms']
        local_rings = rings_per_moltype.get(molname, [])
        for _ in range(count):
            for ring in local_rings:
                global_ring = [atom_offset + (i - 1) for i in ring]
                global_rings.append(global_ring)
            atom_offset += natoms
    return global_rings


def _unwrap_coords(coords: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Unwrap a set of coordinates into a contiguous image using the minimum image convention.

    The first coordinate is taken as a reference; all other coordinates are moved
    by integer multiples of the box vectors so that they remain closest to the
    reference.  Only orthogonal boxes are supported by this helper; for
    triclinic boxes, use the built-in MDAnalysis unwrapping (see notes).

    Parameters
    ----------
    coords : ndarray, shape (N, 3)
        Coordinates of the atoms to unwrap.
    box : ndarray, shape (3,)
        Box lengths along x, y and z (orthogonal cell).

    Returns
    -------
    ndarray
        Unwrapped coordinates of the same shape as ``coords``.

    Notes
    -----
    This function is only used as a last resort if MDAnalysis does not provide
    an ``unwrap`` option or ``apply_PBC``.  For triclinic boxes the built-in
    MDAnalysis unwrapping should be preferred.
    """
    ref = coords[0]
    unwrapped = coords.copy()
    unwrapped[0] = ref
    for i in range(1, coords.shape[0]):
        delta = coords[i] - ref
        # apply minimum image convention on each axis
        delta -= box * np.round(delta / box)
        unwrapped[i] = ref + delta
    return unwrapped


def _unwrap_coords_by_bonds(coords: np.ndarray, box: np.ndarray, edges: list[tuple[int, int]]) -> np.ndarray:
    """Unwrap a set of coordinates into a contiguous image by propagating along a bond graph.

    This is a more robust alternative to :func:`_unwrap_coords` for whole molecules:
    starting from atom 0 as anchor, it traverses the bond graph and places each
    neighbour in the closest periodic image.  Only orthogonal boxes are supported.

    Parameters
    ----------
    coords : ndarray, shape (N, 3)
        Coordinates of the atoms to unwrap (local indexing 0..N-1).
    box : ndarray, shape (3,)
        Box lengths along x, y and z (orthogonal cell).
    edges : list of tuple(int, int)
        Bond list in local indices (0-based) describing connectivity.

    Returns
    -------
    ndarray
        Unwrapped coordinates of the same shape as ``coords``.

    Notes
    -----
    This helper is used to "make the molecule whole" in the output trajectory
    (i.e., to avoid a molecule split across periodic boundaries).  It does not
    change internal geometry: it only chooses the periodic image for each atom.
    """
    if coords.shape[0] == 0:
        return coords

    adj = [[] for _ in range(coords.shape[0])]
    for i, j in edges:
        if 0 <= i < coords.shape[0] and 0 <= j < coords.shape[0]:
            adj[i].append(j)
            adj[j].append(i)

    unwrapped = coords.copy()
    visited = np.zeros(coords.shape[0], dtype=bool)
    visited[0] = True
    stack = [0]

    while stack:
        i = stack.pop()
        ri = unwrapped[i]
        for j in adj[i]:
            if visited[j]:
                continue
            # place j in the closest image relative to i
            delta = coords[j] - coords[i]
            delta -= box * np.round(delta / box)
            unwrapped[j] = ri + delta
            visited[j] = True
            stack.append(j)

    # Safety: if connectivity is incomplete and some atoms remain unvisited,
    # fall back to minimum-image vs atom 0.
    if not np.all(visited):
        ref = unwrapped[0]
        for k in range(coords.shape[0]):
            if visited[k]:
                continue
            delta = coords[k] - coords[0]
            delta -= box * np.round(delta / box)
            unwrapped[k] = ref + delta

    return unwrapped


def _wrap_point_to_primary_cell(point: np.ndarray, box6: np.ndarray) -> np.ndarray:
    """Wrap a single point into the primary unit cell (orthorhombic boxes)."""
    if apply_PBC is not None:
        try:
            return apply_PBC(point, box6)
        except Exception:
            pass
    lx, ly, lz = box6[:3]
    # Avoid divide-by-zero if box is missing
    if lx == 0 or ly == 0 or lz == 0:
        return point
    return np.array([
        point[0] - lx * np.floor(point[0] / lx),
        point[1] - ly * np.floor(point[1] / ly),
        point[2] - lz * np.floor(point[2] / lz),
    ])


def _choose_ring_representative_first_carbon(universe: mda.Universe, ring_global: list[int]) -> int:
    """Choose the representative atom for the ring COM.

    The representative is the first atom in the ring whose mass satisfies:
        1 < mass < 13

    This excludes H-like masses (~1) and excludes heteroatoms with masses > 13.
    If no such atom is present, fall back to the lowest global index.

    Notes
    -----
    We keep the representative's original atom name in the output .gro and .top.
    """
    for idx in ring_global:
        m = float(universe.atoms[idx].mass)
        if (m > 1.0) and (m < 13.0):
            return idx
    return min(ring_global)


def _ring_heteroatoms(universe: mda.Universe, ring_global: list[int], rep_idx: int) -> list[int]:
    """Return heteroatoms in a ring (mass > 13), excluding the representative."""
    heteros = []
    for idx in ring_global:
        if idx == rep_idx:
            continue
        m = float(universe.atoms[idx].mass)
        if m > 13.0:
            heteros.append(idx)
    return heteros


def _rings_in_molecule_from_tpr(mol_atoms: mda.core.groups.AtomGroup, min_size: int = 5, max_size: int = 6):
    """Detect rings (cycles) of heavy atoms within a molecule using only bond information.

    This is the TPR-based replacement for parsing the .top: we build a heavy-atom
    graph from bonds and look for cycles of size 5-6 (configurable).
    """
    # Build graph of heavy atoms
    G = nx.Graph()
    heavy_globals = []
    heavy_set = set()
    for at in mol_atoms:
        name = at.name.strip().upper()
        if name.startswith('H'):
            continue
        heavy_globals.append(at.index)
        heavy_set.add(at.index)
        G.add_node(at.index)

    # No heavy atoms -> no rings
    if not heavy_globals:
        return []

    # Add edges for heavy-heavy bonds
    # Bond information is available through the Universe; here we pull bonds from the molecule slice
    # by filtering universe bonds that touch only atoms in this molecule and are heavy-heavy.
    # If the topology has no bonds, no rings can be detected.
    try:
        mol_bonds = mol_atoms.bonds
    except Exception:
        mol_bonds = []

    for b in mol_bonds:
        i, j = b.indices
        if i in heavy_set and j in heavy_set:
            G.add_edge(i, j)

    if G.number_of_edges() == 0:
        return []

    raw_cycles = nx.cycle_basis(G)
    rings = []
    for cyc in raw_cycles:
        if min_size <= len(cyc) <= max_size:
            rings.append(sorted(cyc))

    # Remove duplicates
    unique = []
    seen = set()
    for ring in rings:
        key = tuple(ring)
        if key not in seen:
            seen.add(key)
            unique.append(ring)

    return unique


def build_molecule_instances_from_tpr(universe: mda.Universe, min_size: int = 5, max_size: int = 6):
    """Build molecule instances directly from the TPR (no .top needed).

    This replaces the .top-driven expansion logic and instead relies on the
    molecule identities present in the TPR (molnums) and the bond graph.

    Returns
    -------
    tuple
        ``(mol_instances, global_rings)`` where:
          - mol_instances is a list of dicts describing each molecule
          - global_rings is a flattened list of all rings in the system
    """
    mol_instances = []
    global_rings = []

    # Prefer the TPR-provided molecule grouping (molnums). If unavailable, fall back to fragments.
    try:
        molecules = universe.atoms.molecules
    except Exception:
        molecules = universe.atoms.fragments

    for mol in molecules:
        # Use a stable, topology-like ordering for "local indices":
        # in GROMACS TPR, atoms of a molecule are typically contiguous and ordered;
        # sorting by global index reproduces that order robustly.
        idxs = np.asarray(sorted(mol.atoms.indices), dtype=int)

        # Build global->local mapping
        global_to_local = {g: i for i, g in enumerate(idxs.tolist())}

        # Build local bond list
        edges = []
        try:
            for b in mol.bonds:
                gi, gj = b.indices
                # Only include bonds fully inside this molecule (should be true for mol.bonds)
                if gi in global_to_local and gj in global_to_local:
                    edges.append((global_to_local[gi], global_to_local[gj]))
        except Exception:
            edges = []

        # Detect rings as global indices directly
        rings_global_list = _rings_in_molecule_from_tpr(mol.atoms, min_size=min_size, max_size=max_size)

        # Build per-ring metadata:
        #  - local indices (for COM-from-unwrapped-molecule coordinates)
        #  - representative atom = first carbon in ring
        #  - heteroatoms kept in output if present
        rings = []
        for ring_g in rings_global_list:
            ok = True
            ring_l = []
            for g in ring_g:
                if g not in global_to_local:
                    ok = False
                    break
                ring_l.append(global_to_local[g])
            if not ok:
                continue

            rep_g = _choose_ring_representative_first_carbon(universe, ring_g)
            rep_l = global_to_local[rep_g]
            hetero_g = _ring_heteroatoms(universe, ring_g, rep_g)
            hetero_l = [global_to_local[h] for h in hetero_g if h in global_to_local]

            rings.append({
                'ring_global': ring_g,
                'ring_local': ring_l,
                'rep_global': rep_g,
                'rep_local': rep_l,
                'hetero_globals': hetero_g,
                'hetero_locals': hetero_l,
            })

        # Store molecule instance
        mol_instances.append({
            'molname': getattr(mol, 'resnames', ['MOL'])[0] if len(getattr(mol, 'resnames', [])) > 0 else 'MOL',
            'global_indices': idxs,
            'global_to_local': global_to_local,
            'edges': edges,
            'rings': rings,
        })

        # Collect system rings
        for r in rings:
            global_rings.append(r['ring_global'])

    # Remove duplicates globally (safe)
    unique = []
    seen = set()
    for ring in global_rings:
        key = tuple(ring)
        if key not in seen:
            seen.add(key)
            unique.append(ring)

    return mol_instances, unique

def _canonical_moleculetype_signature(universe: mda.Universe, idxs_sorted: np.ndarray, edges_local: list[tuple[int, int]]) -> tuple:
    """Build a canonical signature that approximates a 'moleculetype' using TPR-only data.

    The signature is composed of:
      - ordered atom names (local order defined by sorted global indices)
      - ordered bond list in local indices (sorted and normalised)

    This lets us group molecule instances that are the same chemical topology.
    """
    # Atom names in the chosen local order
    names = tuple(universe.atoms[idxs_sorted].names.tolist())

    # Canonicalise edges: (min, max), sort
    canon_edges = []
    for i, j in edges_local:
        a, b = (i, j) if i < j else (j, i)
        canon_edges.append((a, b))
    canon_edges.sort()
    edges = tuple(canon_edges)

    return (names, edges)


def print_detected_rings_per_moleculetype(universe: mda.Universe, mol_instances: list, max_examples_per_type: int = 2):
    """Print detected rings ONCE per moleculetype-like group (not per molecule instance).

    This restores the original intent: show the cycles detected for each molecule type
    as local 1-based indices within the moleculetype, printing only representative examples.

    Notes
    -----
    Because we are in TPR-only mode, 'moleculetype' is approximated by grouping
    molecules with the same ordered atom names and the same bond list in local indices.
    """
    # Group molecule instances by a canonical signature
    groups = {}  # signature -> list of indices into mol_instances
    for mi, inst in enumerate(mol_instances):
        idxs = inst['global_indices']
        # idxs are already sorted in build_molecule_instances_from_tpr (per our earlier recommendation)
        sig = _canonical_moleculetype_signature(universe, idxs, inst['edges'])
        groups.setdefault(sig, []).append(mi)

    # Print once per group
    for gi, (sig, members) in enumerate(groups.items()):
        # pick representative instance
        rep_inst = mol_instances[members[0]]
        molname = rep_inst.get('molname', 'MOL')

        rings = rep_inst.get('rings', [])
        if not rings:
            continue

        print(f"  Molécula {molname}: {len(rings)} anillo(s) detectado(s), ejemplos:")

        for r_i, r in enumerate(rings[:max_examples_per_type]):
            ring_l0 = r['ring_local']
            ring_l1 = [x + 1 for x in ring_l0]
            print(f"    {ring_l1}")

        # Optional: also show representative/hetero info for first examples
        for r_i, r in enumerate(rings[:max_examples_per_type]):
            ring_l0 = r['ring_local']
            ring_l1 = [x + 1 for x in ring_l0]
            rep_l1 = r['rep_local'] + 1
            hetero_l1 = [x + 1 for x in r.get('hetero_locals', [])]

            print(f"    Anillo {r_i} detalle (local 1-based): ring={ring_l1}, COM_rep={rep_l1}, hetero={hetero_l1 if hetero_l1 else '[]'}")

def print_detected_rings_summary(universe: mda.Universe, mol_instances: list, max_examples_per_mol: int = 2):
    """Print detected rings per molecule with 1-based indices within the moleculetype.

    We define "local index within moleculetype" as the 1-based position of an atom
    within the molecule instance (ordered by global atom index). For GROMACS TPRs
    this matches the typical moleculetype atom numbering.

    It prints:
      - per molecule instance: number of rings detected
      - per ring example: the local 1-based indices (moleculetype-like)
      - per ring atom: local index + global index + resid/resname/atomname
      - representative atom and heteroatoms (useful for the heteroaromatic rule)
    """
    printed_any = False

    for mi, inst in enumerate(mol_instances):
        rings = inst.get('rings', [])
        if not rings:
            continue

        printed_any = True
        molname = inst.get('molname', 'MOL')
        print(f"  Molécula {molname} (instancia {mi}): {len(rings)} anillo(s) detectado(s)")

        # show up to max_examples_per_mol rings
        for r_i, r in enumerate(rings[:max_examples_per_mol]):
            ring_g = r['ring_global']
            ring_l0 = r['ring_local']          # 0-based local
            ring_l1 = [x + 1 for x in ring_l0] # 1-based local (moleculetype-like)

            rep_g = r.get('rep_global', min(ring_g))
            rep_l1 = r.get('rep_local', inst['global_to_local'][rep_g]) + 1
            hetero_g = r.get('hetero_globals', [])
            hetero_l1 = [inst['global_to_local'][h] + 1 for h in hetero_g if h in inst['global_to_local']]

            print(f"    Anillo ejemplo {r_i}: índices locales (1-based, moleculetype) = {ring_l1}")

            # Print a per-atom line showing local and global ids and residue identity
            for loc1, g in sorted(zip(ring_l1, ring_g), key=lambda t: t[0]):
                a = universe.atoms[g]
                print(f"      local={loc1:4d}  global={g:6d}  resid={a.resid:5d}  resname={a.resname:>6s}  atom={a.name:>6s}")

            rep_a = universe.atoms[rep_g]
            print(f"      -> COM representante: local={rep_l1}  global={rep_g}  resid={rep_a.resid}  resname={rep_a.resname}  atom={rep_a.name}")

            if hetero_g:
                het_items = []
                for h, hl1 in sorted(zip(hetero_g, hetero_l1), key=lambda t: t[1]):
                    ha = universe.atoms[h]
                    het_items.append(f"local={hl1},global={h}:{ha.name}(resid {ha.resid})")
                print(f"      -> Heteroátomos conservados: " + "; ".join(het_items))
            else:
                print(f"      -> Heteroátomos conservados: ninguno")

    if not printed_any:
        print("  (No se han detectado anillos en ninguna molécula.)")

def _build_ring_keep_remove_sets(universe: mda.Universe, global_rings: list[list[int]], mol_instances: list | None = None):
    """Build atom keep/remove sets for ring reduction with heteroatom preservation.

    Notes
    -----
    Requirement implemented:
      - For heteroaromatic rings, keep BOTH:
          * the ring COM (written into the representative atom)
          * the original heteroatom position(s)
      - The representative atom for the COM is the first carbon in the ring so
        the COM atom name in the .gro output is carbon-like and does not duplicate
        the heteroatom name.

    This function is careful in case rings share atoms: an atom is removed only
    if it is not required to be kept by any ring.
    """
    keep_set = set()
    remove_candidates = set()
    ring_infos = []

    if mol_instances is not None and len(mol_instances) > 0:
        # Use the richer ring metadata if available
        for inst in mol_instances:
            for r in inst['rings']:
                ring_g = r['ring_global']
                rep_g = r['rep_global']
                hetero_g = r['hetero_globals']
                keep_set.add(rep_g)
                for h in hetero_g:
                    keep_set.add(h)
                for idx in ring_g:
                    if idx != rep_g and idx not in hetero_g:
                        remove_candidates.add(idx)
                ring_infos.append({
                    'ring_global': ring_g,
                    'rep_global': rep_g,
                    'hetero_globals': hetero_g,
                })
    else:
        # Derive ring representative/heteroatoms directly from masses 
        for ring_g in global_rings:
            rep_g = _choose_ring_representative_first_carbon(universe, ring_g)
            hetero_g = _ring_heteroatoms(universe, ring_g, rep_g)
            keep_set.add(rep_g)
            for h in hetero_g:
                keep_set.add(h)
            for idx in ring_g:
                if idx != rep_g and idx not in hetero_g:
                    remove_candidates.add(idx)
            ring_infos.append({
                'ring_global': ring_g,
                'rep_global': rep_g,
                'hetero_globals': hetero_g,
            })

    # Final remove set excludes anything that must be kept
    remove_set = set([idx for idx in remove_candidates if idx not in keep_set])

    return ring_infos, keep_set, remove_set


def compute_and_write_reduced(universe: mda.Universe, global_rings: list, out_xtc: str, out_gro: str | None = None, mol_instances: list | None = None):
    """Compute COMs of aromatic rings and write a reduced trajectory.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        The universe containing both topology and trajectory.
    global_rings : list of list of int
        A list of rings, each defined by global atom indices (0-based).
    out_xtc : str
        Path to the output trajectory file (XTC format).
    out_gro : str or None, optional
        Path to the output structure file (GRO format).  If ``None``, the
        structure is not written.
    mol_instances : list or None, optional
        Per-molecule metadata built from the topology parser.  If provided, the
        code will also "make whole" the full molecule (orthorhombic PBC) in the
        output trajectory, i.e. it reconstructs the complete molecule as a
        contiguous object in a single periodic image, not just the ring COM.

    Notes
    -----
    For each ring, the heavy atoms are unwrapped (if possible) before computing
    the centre of mass.  The resulting COM is then wrapped back into the
    primary unit cell.  Only a single representative atom (the one with the
    smallest index in the ring) is retained to carry the COM; all other heavy
    atoms belonging to the ring are removed from the output.  Hydrogens and
    atoms outside rings are untouched.
    """
    n_atoms = len(universe.atoms)

    # Determine which atoms to keep/remove (one representative per ring, keep heteroatoms)
    ring_infos, keep_ring_atoms, non_rep_atoms = _build_ring_keep_remove_sets(
        universe, global_rings, mol_instances=mol_instances
    )

    keep_indices = np.array([i for i in range(n_atoms) if i not in non_rep_atoms], dtype=int)
    keep_indices.sort()

    print(f"Número total de átomos original : {n_atoms}")
    print(f"Átomos eliminados               : {len(non_rep_atoms)}")
    print(f"Átomos en la trayectoria final  : {len(keep_indices)}")

    ag_keep = universe.atoms[keep_indices]

    # If no molecule instances were provided, we fall back to the original behaviour:
    # only ring COMs are computed with PBC awareness, and the rest of the system is
    # written exactly as read.
    if mol_instances is None:
        mol_instances = []

    # Write reduced trajectory
    with mda.Writer(out_xtc, n_atoms=len(keep_indices)) as W:
        for ts in universe.trajectory:
            # copy full positions for this frame
            new_pos = universe.atoms.positions.copy()
            box = ts.dimensions  # [lx, ly, lz, alpha, beta, gamma]
            ortho_box = box[:3]

            # If mol_instances is provided, we will also "make whole" the full
            # molecule (orthorhombic PBC) before inserting COMs.
            if mol_instances:
                for inst in mol_instances:
                    # Molecules without rings do not need to be reconstructed here.
                    if not inst['rings']:
                        continue

                    idxs = inst['global_indices']
                    coords_mol = new_pos[idxs].astype(float)

                    # Unwrap the full molecule using the bond graph (orthorhombic only).
                    coords_unwrapped = _unwrap_coords_by_bonds(coords_mol, ortho_box, inst['edges'])

                    # Anchor the whole molecule translation to the first ring COM:
                    # compute COM in unwrapped coordinates, then shift the whole molecule
                    # so that this COM lies in the primary unit cell.
                    r0 = inst['rings'][0]
                    ring0_local = r0['ring_local']
                    ring0_global = r0['ring_global']
                    masses0 = universe.atoms[ring0_global].masses
                    com0_unwrapped = np.average(coords_unwrapped[ring0_local], axis=0, weights=masses0)
                    com0_wrapped = _wrap_point_to_primary_cell(com0_unwrapped, box)
                    shift = com0_wrapped - com0_unwrapped
                    coords_unwrapped += shift

                    # Now compute COM for each ring in this molecule from the unwrapped coords
                    # and write it into the representative atom position (global). Heteroatoms
                    # are kept with their (reconstructed) positions.
                    for r in inst['rings']:
                        ring_local = r['ring_local']
                        ring_global = r['ring_global']
                        masses = universe.atoms[ring_global].masses
                        com_unwrapped = np.average(coords_unwrapped[ring_local], axis=0, weights=masses)
                        # Keep COM in primary unit cell for consistency
                        com_wrapped = _wrap_point_to_primary_cell(com_unwrapped, box)
                        rep_idx = r['rep_global']
                        rep_local = r['rep_local']
                        coords_unwrapped[rep_local] = com_wrapped

                    # Store the reconstructed molecule coordinates back into the full coordinate array
                    new_pos[idxs] = coords_unwrapped
            else:
                # Original behaviour: only compute COM for each ring
                for info in ring_infos:
                    ring = info['ring_global']
                    rep_idx = info['rep_global']

                    ag_ring = universe.atoms[ring]
                    try:
                        # Try to unwrap the ring before computing the COM.  The
                        # ``unwrap`` flag became available in MDAnalysis 1.0.0; if
                        # unsupported a TypeError is raised.
                        com = ag_ring.center_of_mass(unwrap=True)
                    except Exception:
                        # Fallback: manually unwrap using minimum image convention
                        coords = ag_ring.positions.astype(float)
                        masses = ag_ring.masses
                        unwrapped = _unwrap_coords(coords, ortho_box)
                        com = np.average(unwrapped, axis=0, weights=masses)

                    # Wrap the COM back into the primary unit cell
                    com = _wrap_point_to_primary_cell(np.asarray(com, dtype=float), box)

                    # Set representative atom position to COM
                    new_pos[rep_idx] = com

            # Build reduced coordinate array and write
            reduced_pos = new_pos[keep_indices]
            ag_keep.positions = reduced_pos
            W.write(ag_keep)

    # Optionally write structure
    if out_gro is not None:
        # reset to first frame
        universe.trajectory[0]
        ts = universe.trajectory.ts
        box = ts.dimensions
        ortho_box = box[:3]
        new_pos = universe.atoms.positions.copy()

        if mol_instances:
            for inst in mol_instances:
                if not inst['rings']:
                    continue

                idxs = inst['global_indices']
                coords_mol = new_pos[idxs].astype(float)
                coords_unwrapped = _unwrap_coords_by_bonds(coords_mol, ortho_box, inst['edges'])

                r0 = inst['rings'][0]
                ring0_local = r0['ring_local']
                ring0_global = r0['ring_global']
                masses0 = universe.atoms[ring0_global].masses
                com0_unwrapped = np.average(coords_unwrapped[ring0_local], axis=0, weights=masses0)
                com0_wrapped = _wrap_point_to_primary_cell(com0_unwrapped, box)
                shift = com0_wrapped - com0_unwrapped
                coords_unwrapped += shift

                for r in inst['rings']:
                    ring_local = r['ring_local']
                    ring_global = r['ring_global']
                    masses = universe.atoms[ring_global].masses
                    com_unwrapped = np.average(coords_unwrapped[ring_local], axis=0, weights=masses)
                    com_wrapped = _wrap_point_to_primary_cell(com_unwrapped, box)
                    rep_local = r['rep_local']
                    coords_unwrapped[rep_local] = com_wrapped

                new_pos[idxs] = coords_unwrapped
        else:
            for info in ring_infos:
                ring = info['ring_global']
                rep_idx = info['rep_global']

                ag_ring = universe.atoms[ring]
                try:
                    com = ag_ring.center_of_mass(unwrap=True)
                except Exception:
                    coords = ag_ring.positions.astype(float)
                    masses = ag_ring.masses
                    unwrapped = _unwrap_coords(coords, ortho_box)
                    com = np.average(unwrapped, axis=0, weights=masses)

                com = _wrap_point_to_primary_cell(np.asarray(com, dtype=float), box)
                new_pos[rep_idx] = com

        reduced_pos = new_pos[keep_indices]
        ag_keep.positions = reduced_pos
        with mda.Writer(out_gro, n_atoms=len(keep_indices)) as Wgro:
            Wgro.write(ag_keep)
        print(f"Estructura reducida escrita en: {out_gro}")


def main():
    # New usage (no .top required):
    #     python com_aromatic_rings.py topol.tpr traj.xtc out.xtc
    #
    # Backwards-compatible usage (second argument is ignored, kept only so old calls don't break):
    #     python com_aromatic_rings.py topol.tpr topol.top traj.xtc out.xtc
    if len(sys.argv) not in (4, 5):
        print("Uso: python com_aromatic_rings.py topol.tpr traj.xtc out.xtc")
        print(" (compatibilidad) python com_aromatic_rings.py topol.tpr topol.top traj.xtc out.xtc")
        sys.exit(1)

    tprfile = sys.argv[1]

    if len(sys.argv) == 5:
        # Keep accepting a .top argument but do not read it anymore.
        topfile = sys.argv[2]
        xtcfile = sys.argv[3]
        out_xtc = sys.argv[4]
        print(f"NOTA: argumento .top recibido pero NO se usa: {topfile}")
    else:
        xtcfile = sys.argv[2]
        out_xtc = sys.argv[3]

    out_gro = out_xtc.replace('.xtc', '.gro')

    print(f"Leyendo coordenadas y trayectoria desde:")
    print(f"  TPR : {tprfile}")
    print(f"  XTC : {xtcfile}")
    u = mda.Universe(tprfile, xtcfile)

    # Build per-molecule instances directly from the TPR, and find all rings globally
    print("Detectando anillos de átomos pesados (5-6 miembros) desde el TPR...")
    mol_instances, global_rings = build_molecule_instances_from_tpr(u, min_size=5, max_size=6)
    print(f"Total de anillos en el sistema: {len(global_rings)}")

    print("Detectando anillos de átomos pesados (5-6 miembros) por moleculetype (TPR-only)...")
    print_detected_rings_per_moleculetype(u, mol_instances, max_examples_per_type=2)

    print(f"Escribiendo trayectoria reducida en: {out_xtc}")
    print(f"Escribiendo estructura reducida en: {out_gro}")
    compute_and_write_reduced(u, global_rings, out_xtc, out_gro, mol_instances=mol_instances)
    print("Listo.")

    # OPTIONAL: write a modified .top (requires an input standalone .top)
    if len(sys.argv) == 5:
        # Here topfile is available as sys.argv[2] (even though you don't use it for rings anymore)
        # We can still use it to generate a consistent reduced topology.
        print("Escribiendo topología reducida (.top) coherente con la trayectoria...")
        # Build rings_per_moltype from the standalone .top
        moltypes_text, mol_order_text = parse_topology(topfile)
        rings_per_moltype_text = find_rings_in_moltypes(moltypes_text, min_size=5, max_size=6)

        # Sanity check: ring detection in .top vs TPR-only may differ; warn/abort if so.
        # This is a minimal guard to avoid producing a .top inconsistent with the reduced trajectory.
        if len(rings_per_moltype_text) == 0:
            print("ADVERTENCIA: no se han detectado anillos en el .top; el .top reducido puede no ser coherente con el .xtc reducido.")
        

        out_top = out_xtc.replace(".xtc", ".top")
        write_modified_top(
            in_top=topfile,
            out_top=out_top,
            rings_per_moltype=rings_per_moltype_text,
            remove_sections={"angles", "dihedrals"},
        )
        print(f"Topología reducida escrita en: {out_top}")



if __name__ == '__main__':
    main()

