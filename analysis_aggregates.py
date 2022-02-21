#!/usr/bin/env python3

import MDAnalysis
import numpy as np
import sys
import re
import time

def angle(pos1, pos2, pos3):
    '''Angle 1-2-3'''
    
    v21 = pos1 - pos2
    v23 = pos3 - pos2
    
    cosA = v21.dot(v23)/np.linalg.norm(v21)/np.linalg.norm(v23)
    
    return np.arccos(cosA)

def identify_blocks(bmatrix):
    '''Identify blocks in a block matrix'''
    
    blocks = []
    
    i = 0
    while i < len(bmatrix):
        arrow = bmatrix[i]
        if bmatrix[i,i] == 0:
            shift = 1
        else:
            shift = 0
        isize = len(arrow[arrow != 0]) + shift
        blocks.append(isize)
        i += isize
        
    return blocks

def analyze_cluster(cluster,verbose=False):
    
    from MDAnalysis.analysis.distances import distance_array 
    
    bmatrix = np.zeros([cluster.n_residues,cluster.n_residues])
    resids  = [ res.resid for res in cluster.residues ]
    base    = []
    for i,resid in enumerate(resids):
        if resid not in base:
            base.append(resid)
        if i>=len(resids):
            break
        for resid2 in resids[i+1:]:
            d = distance_array(u.residues[resid].atoms[11].position,
                            u.residues[resid2].atoms[11].position)[0,0]
            if d<3.5:
                if verbose:
                    print(f'Possible connection: {resid:<3}, {resid2:<3}, d={d:6.2f}')
                a1 = angle(u.residues[resid].atoms[11].position,
                        u.residues[resid2].atoms[11].position,
                        u.residues[resid2].atoms[12].position)
                a2 = angle(u.residues[resid2].atoms[11].position,
                        u.residues[resid].atoms[11].position,
                        u.residues[resid].atoms[12].position)
                if a1*180/np.pi < 30.:
                    if verbose:
                        print(f'   Donor {resid2:3} -- Acceptor: {resid :3}  |  Angle: {a1*180/np.pi:6.2f}')
                    # Add resid to base if not present
                    if resid2 not in base:
                        base.append(resid2)
                    # Get base indices for each res
                    ib1 = base.index(resid)
                    ib2 = base.index(resid2)
                    bmatrix[ib1,ib2] = 2
                    bmatrix[ib2,ib1] = 1
                elif a2*180/np.pi < 30.:
                    if verbose:
                        print(f'   Donor {resid :3} -- Acceptor: {resid2:3}  |  Angle: {a2*180/np.pi:6.2f}')
                    # Add resid to base if not present
                    if resid2 not in base:
                        base.append(resid2)
                    # Get base indices for each res
                    ib1 = base.index(resid)
                    ib2 = base.index(resid2)
                    bmatrix[ib1,ib2] = 1
                    bmatrix[ib2,ib1] = 2
                else:
                    if verbose:
                        print(f'   No Hbon  |  Angles: {a1*180/np.pi:6.2f}, {a2*180/np.pi:6.2f}')
                        
            
    iblocks = identify_blocks(bmatrix)
        
    if len(iblocks) > 1:
        #print('Split block of size', cluster.n_residues, 'into bloks of size',iblocks)
        i0 = 0
        blocks = []
        for iblock in iblocks:
            blocks.append(base[i0:i0+iblock])
            i0 += iblock
    
        # Split cluster into new blocks
        cluster_list = cluster.split('residue')
        resids_list = [ res.residues[0].resid for res in cluster_list ]
        
        clusters = []
        for block in blocks:
            i = resids_list.index(block[0])
            cluster = cluster_list[i]
            for r in block[1:]:
                i = resids_list.index(r)
                cluster |= cluster_list[i]
            clusters.append(cluster)
    else:
        clusters = [cluster]
    
    
    return clusters




def atname2element(name, mass):
    """Turn atomname from MD topology into element name, using the mass to guide to assigment
    
    Input:
     name, string: input name (from MD topology)
     mass, float : atomic mass

    Output:
     name, string: output name (converted if possible, or same as input otherwise)
     
    Note:
     It uses an array with atomic masses per mayor isotopes, instead of average atomic masses. 
     Since it finds the closest mass to the one given, it should no be an issue.
     Requires qcelemental
     """

    import qcelemental as qcel

    # Now uses masses to confirm the name assigment

    # Name from mass
    masses = np.array([qcel.periodictable.to_mass(i) for i in range(118)])
    # Find element with closer mass and compare with input name
    name_found = False
    while not name_found:
        Z = (np.abs(masses - mass)).argmin()
        if abs(masses[Z] - mass) > 10:
            print("Name not found. Returning input")
            return name
        name_qcel = qcel.periodictable.to_E(Z)
        if len(name_qcel) > len(name):
            # This name is not possible: skip picking from masses
            masses[Z] = -1.0
        elif name_qcel.upper() == name[: len(name_qcel)].upper():
            name_found = True
        else:
            # This name is not possible: skip picking from masses
            masses[Z] = -1.0

    return name_qcel

def find_cluster_atoms(i,u,sel_command):
    '''Find cluster around residue i. 
Connected residues are obtained with select_atoms using sel_command
    
    Input:
     i: residue index
     u, MDAnalysis.Universe: the universe where the residue lives
     sel_command, string: selection command to get the neighbors. 
                          The residue index is replaced by i
       Examples:
        sel_command = f'same residue as name O1 and (around {rcutoff} (name O1 and resid {i}))'
        sel_command = f'same residue (around {rcutoff} (resid {i}))'

    Output:
     [...], list of ints: list of residue indices 
'''

    sel_command_ = sel_command.replace('{i}',str(i))
    cluster = u.select_atoms(sel_command_)
    cluster_size = -1
    resids = [r.resid for r in cluster.residues]
    resids_ = [i]
    while cluster.n_residues > cluster_size:
        cluster_size = cluster.n_residues
        for r in resids:
            if r in resids_:
                continue
            sel_command_ = sel_command.replace('{i}',str(r))
            cluster |= u.select_atoms(sel_command_)
        resids_ = resids
        resids  = [r.resid for r in cluster.residues]
        
    return cluster


def find_cluster(i,u,sel_command):
    '''Find cluster around residue i. 
Connected residues are obtained with select_atoms using sel_command
    
    Input:
     i: residue index
     u, MDAnalysis.Universe: the universe where the residue lives
     sel_command, string: selection command to get the neighbors. 
                          The residue index is replaced by i
       Examples:
        sel_command = f'same residue as name O1 and (around {rcutoff} (name O1 and resid {i}))'
        sel_command = f'same residue (around {rcutoff} (resid {i}))'

    Output:
     [...], list of ints: list of residue indices 
'''

    # Initialize cluster with the residue index
    cluster = [i]
    # Initialize connetions with neighbors around the residue
    sel_command_ = sel_command.replace('{i}',str(i))
    sel = u.select_atoms(sel_command_)
    connections = [rx.resid for rx in sel.residues]
    # Underscores variables used to save info of the previous step
    connections_ = set()
    while len(connections) > 0:
        cluster_ = set(cluster)
        # Iterate over all connections
        for r in connections:
            if r in connections_:
                continue
            # Get connected residues
            sel_command_ = sel_command.replace('{i}',str(r))
            sel = u.select_atoms(sel_command_)
            cluster += [rx.resid for rx in sel.residues] + [r]
        # with set(), only unique indices retained
        cluster = list(set(cluster))
        # Indices in the cluster are removed from connections
        connections_ = set(connections)
        connections = set(cluster) - cluster_
        
    return sorted(cluster)


def cluster_pop(clusters,MAX_SIZE=9):
    '''Get the population of clusters with different size'''
    
    # Initialize population array
    populations = np.zeros(MAX_SIZE,dtype=int)
    # Generate array with all cluster sizes (dimers and beyond)
    cluster_type = [ cluster.n_residues for cluster in clusters ]
    for isize in sorted(set(cluster_type)):
        if isize>MAX_SIZE:
            print(f'WARNNG: MAX_SIZE exceeded ({isize})')
            populations[MAX_SIZE-1] = cluster_type.count(isize)
        else:
            # Note that smaller cluster has size=2
            populations[isize-2] = cluster_type.count(isize)
        
    return populations


def compact_atoms(box,atoms):

    # Fist, wrap
    atoms.wrap()

    # Ref position is the center_of_geometry of QM layer
    Rref = atoms[0].position #center_of_geometry()
    
    # Iterate over all atoms
    for atom in atoms:
        for ix in range(3):
            d = atom.position[ix] - Rref[ix]
            if abs(d) > box[ix]/2:
                v = atom.position
                if (d>0):
                    v[ix] -= box[ix]
                else:
                    v[ix] += box[ix]
                atom.position =  v
                    
    return None



if __name__ == "__main__":
    import argparse
    ti = time.time()

    # Input parser. Set flags
    parser = argparse.ArgumentParser(
        description="Analysis of aggregates through an MD trajectory."
    )
    parser.add_argument("-f", metavar="file.trr", help="Trajectory file", required=True)
    parser.add_argument(
        "-f_fmt", metavar="trr", help="Format of trajectory file", required=False
    )
    parser.add_argument(
        "-s", metavar="file.tpr", help="Binary topoly file", required=True
    )
    parser.add_argument(
        "-s_fmt", metavar="tpr", help="Format of topoly file", required=False
    )
    parser.add_argument(
        "-sel",
        metavar="select_string",
        help="Selection command to identify neighbors",
        default='same residue as name O1 and (around 3.5 (name O1 and resid {i})) or resid {i}',
        required=False,
    )
    parser.add_argument(
        "-b",
        metavar="<time>",
        help="First frame (ps) to read from trajectory",
        type=float,
        default=-1.0,
    )
    parser.add_argument(
        "-e",
        metavar="<time>",
        help="Last frame (ps) to read from trajectory",
        type=float,
        default=-1.0,
    )
    parser.add_argument(
        "-dt",
        metavar="<time>",
        help="Only use frame when t MOD dt = first time (ps)",
        type=float,
        default=-1.0,
    )
    parser.add_argument(
        "-fixnames",
        action="store_true",
        help="Try to convert atomnames into element names (WARNING: this might work unexpectedly)",
        default=False,
    )
    # Parse input
    args = parser.parse_args()

    # Get topology and coordinates
    u = MDAnalysis.Universe(
        args.s, args.f, topology_format=args.s_fmt, format=args.f_fmt
    )

    # Fix names if requested (do on universe only once)
    if args.fixnames:
        for atom in u.atoms:
            atom.name = atname2element(atom.name, atom.mass)

    # Set some defaults to slice the traj
    if args.b == -1.0:
        tini = u.trajectory.time
    else:
        tini = args.b
    if args.e == -1.0:
        tfin = u.trajectory.totaltime
    else:
        tfin = args.e
    if args.dt == -1.0:
        dt = u.trajectory.dt
    else:
        dt = args.dt

    clusters_traj = []
    pops_traj = []
    for istp in range(u.trajectory.n_frames):
        # Determine whether using the frame or not
        if u.trajectory.n_frames == 1:
            # Always use the frame when there is only one
            pass
        elif u.trajectory.time > tfin + 0.0001:
            break
        elif u.trajectory.time % dt > 0.0001:
            continue
        elif u.trajectory.time < tini:
            continue

        u.trajectory.next()

        # Generate all clusters for this snapshot
        clusters = []
        selected_residue = u.select_atoms('resname X')
        for res in u.residues:
            i = res.resid
            if res in selected_residue:
                continue
            cluster = find_cluster_atoms(i,u,args.sel)
            if cluster.n_residues > 1 and cluster not in clusters:
                #clusters.append(refine_cluster(cluster))
                clusters.append(cluster)
            selected_residue |= cluster
            
        # Accumulate clusters along the trajectory    
        clusters_traj.append(clusters)
        
        clusters_add = []
        clusters_rm  = []
        for cluster in clusters:
            compact_atoms(u.dimensions[:3],cluster.atoms)
            refined_cluster = analyze_cluster(cluster)
            if len(refined_cluster) == 1:
                continue
            else:
                clusters_rm.append(cluster)
                clusters_add += refined_cluster
        for cluster_rm in clusters_rm:
            clusters.remove(cluster_rm)
        for cluster_add in clusters_add:
            clusters.append(cluster_add)
        #for i,cluster in enumerate(clusters):
            #if cluster.n_residues == 2:
                #compact_atoms(u.dimensions[:3],cluster.atoms)
                #cluster.atoms.write(f'Step{istp:03g}_Dimer{i:03g}.gro')
                
        # Get populations after refining blocks
        pops = cluster_pop(clusters)
        pops_traj.append(pops)
            
        
    # Analyze populations
    pops_traj = np.array(pops_traj)
    pops_traj = pops_traj.transpose()
    for i,ctype in enumerate(pops_traj[:-1]):
        p_av  = ctype.mean()
        p_dev = np.std(ctype)
        n = len(ctype)
        p1_av = ctype[:int(n/2)].mean()
        p2_av = ctype[int(n/2):].mean()
        dev = abs(p1_av-p2_av)/2
        print(f'Size: {i+2} -- Average pop: {p_av:6.2f} +-{p_dev:6.2f}  | Blocks +-{dev:5.2f}  (+-{dev/p_av*100:5.1f}%))')
    
    # Analyze clusters with size>MAX_SIZE
    ctype = pops_traj[-1]
    p_av  = ctype.mean()
    p_dev = np.std(ctype)
    n = len(ctype)
    p1_av = ctype[:int(n/2)].mean()
    p2_av = ctype[int(n/2):].mean()
    dev = abs(p1_av-p2_av)/2
    print(f'Size >{len(pops_traj)} -- Average pop: {p_av:6.2f} +-{p_dev:6.2f}  | Blocks +-{dev:5.2f}  (+-{dev/p_av*100:5.1f}%))')

    
    # Report elapsed time
    tf = time.time()
    print(f'Elapsed time       : {tf-ti:8.2f} s')
        


