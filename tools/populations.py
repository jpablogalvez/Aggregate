import sys
import os
import networkx as nx
from networkx.algorithms.isomorphism.vf2pp import vf2pp_is_isomorphic
import matplotlib.pyplot as plt

##############################################################################
#  Función para clasificar grafos no dirigidos, omitiendo el primer bloque
##############################################################################
def clasificar_grafos_streaming_nodir(ruta_fichero: str):
    """Clasifica grafos no dirigidos en clases de isomorfía.

    Se salta la primera matriz (grafo base) pero la conserva para
    identificar las aristas intramoleculares.  Para matrices completas no
    se realiza la simetrización por OR, de modo que únicamente se añaden
    aristas cuando la posición (i,j) de la matriz contiene 'T' y i < j.
    Después de construir cada grafo, se crea un atributo 'tag' que contiene
    (etiqueta, grado) para poder usar `vf2pp_is_isomorphic` con un único
    atributo `node_label`.

    Args:
        ruta_fichero: ruta al fichero de entrada

    Returns:
        tuple (clases, asignaciones), donde `clases` es una lista de
        diccionarios con llaves `nx_graph`, `labels`, `matriz`,
        `base_graph` y `count`, y `asignaciones` es una lista de tuplas
        (id_grafo, id_clase).
    """
    clases = []
    asignaciones = []
    grafo_id = 0

    with open(ruta_fichero, "r") as f:
        # 1) Leer las etiquetas de la PRIMERA línea
        labels = f.readline().strip().split()
        n = len(labels)
        if n == 0:
            raise ValueError("No se encontraron etiquetas en la primera línea.")

        primer_matriz_leida = False
        grafo_base = None

        # 2) Leer las matrices de adyacencia hasa llegar al final del fichero y clasificarlas on-the-fly
        while True:
            # Leer la primera línea del siguiente bloque.  Si no hay más líneas,
            # terminamos la lectura.  Independientemente de si la matriz es
            # triangular o completa, se tratará como triangular tomando sólo
            # los valores de la parte inferior.
            first = f.readline()
            if not first:
                break
            tokens0 = first.strip().split()
            # Si encontramos una línea vacía, salimos del bucle
            if not tokens0:
                break

            # Siempre leemos la matriz en formato triangular inferior.  Incluso
            # si la fila contiene n columnas (matriz completa), sólo se usa la
            # primera columna para i=0, las dos primeras para i=1, etc.  Esto
            # reproduce el comportamiento de la versión antigua y elimina la
            # necesidad de un if/else.
            matriz_tri = []
            # Fila 0: tomar solo el primer valor, ignorando columnas extra
            matriz_tri.append([x == "T" for x in tokens0[:1]])
            # Leer el resto de filas
            for i in range(1, n):
                line = f.readline()
                # Si no hay más líneas (fin de fichero), salimos del bucle principal
                if not line:
                    return clases, asignaciones
                # Parsear la línea
                toks = line.strip().split()
                if len(toks) < i + 1:
                    raise ValueError(
                        f"Se esperaban al menos {i+1} valores en fila {i}, pero hubo {len(toks)}"
                    )
                # Tomar sólo los primeros i+1 valores (parte triangular)
                matriz_tri.append([x == "T" for x in toks[: i + 1]])
            # Reconstruir la matriz completa copiando el valor triangular a ambos lados
            complete_matrix = [[False] * n for _ in range(n)]
            for i in range(n):
                for j in range(i + 1):
                    val = matriz_tri[i][j]
                    complete_matrix[i][j] = val
                    complete_matrix[j][i] = val

            # Construir el grafo usando sólo los valores de la matriz para i < j
            G = nx.Graph()
            # Usamos identificadores únicos por índice, y guardamos la etiqueta en un atributo
            for idx, lbl in enumerate(labels):
                G.add_node(idx, label=lbl)

            # Añadir aristas
            for i in range(n):
                for j in range(i + 1, n):
                    if complete_matrix[i][j]:
                        G.add_edge(i, j)

            # Añadir el atributo combinado (etiqueta, grado)
            for node in G.nodes:
                label = G.nodes[node]["label"]
                deg = G.degree[node]
                G.nodes[node]["tag"] = (label, deg)

            # La primera matriz se usa como grafo base intramolecular
            if not primer_matriz_leida:
                grafo_base = G.copy()
                primer_matriz_leida = True
                continue

            grafo_id += 1
            asignada = False
            # Comparar con las clases existentes usando el atributo 'tag'
            for j, clase in enumerate(clases):
                rep = clase["nx_graph"]
                if vf2pp_is_isomorphic(G, rep, node_label="tag"):
                    # Es isomorfo a la clase j
                    clase["count"] += 1
                    asignaciones.append((grafo_id, j + 1))
                    asignada = True
                    break

            if not asignada:
                # No es isomorfo a ninguna de las clases existentes
                # Guardamos la información
                clases.append({
                    "nx_graph": G,
                    "labels": labels,
                    "matriz": complete_matrix,
                    "base_graph": grafo_base,
                    "count": 1
                })
                # La nueva clase tiene índice = len(clases)
                asignaciones.append((grafo_id, len(clases)))
    return clases, asignaciones

##############################################################################
#  Función para clasificar grafos dirigidos, omitiendo el primer bloque
##############################################################################
def clasificar_grafos_streaming_dir(ruta_fichero: str):
    """Clasifica grafos dirigidos en clases de isomorfía.

    Se salta la primera matriz (grafo base) para utilizarla como referencia
    intramolecular.  Cada nodo recibe un atributo `tag` que contiene la
    etiqueta, el grado de entrada y el grado de salida.  Este atributo se
    utiliza como `node_label` en `vf2pp_is_isomorphic`.

    Args:
        ruta_fichero: ruta al fichero de entrada

    Returns:
        tuple (clases, asignaciones) igual que la función de no dirigidos.
    """
    clases = []
    asignaciones = []
    grafo_id = 0
    with open(ruta_fichero, "r") as f:
        # 1) Leer las etiquetas de la PRIMERA línea
        labels = f.readline().strip().split()
        n = len(labels)
        if n == 0:
            raise ValueError("No se encontraron etiquetas en la primera línea.")

        primer_matriz_leida = False
        grafo_base = None

        # 2) Leer los bloques de NxN líneas hasta llegar al final del fichero
        while True:
            matriz = []
            for _ in range(n):
                line = f.readline()
                # Si no hay más líneas (fin de fichero), salimos del bucle principal
                if not line:
                    break
                # Parsear la línea
                tokens = line.strip().split()
                if len(tokens) != n:
                   # Si se encuentra una fila vacía o con columnas distintas, cortamos
                   if not tokens:
                       break
                   raise ValueError(
                        f"Número de columnas distinto de {n} en una de las filas.\n"
                        f"Fila leída: {tokens}"
                    )
                row = [True if x == "T" else False for x in tokens]
                matriz.append(row)

            # Si la matriz quedó incompleta (por haber salido del bucle antes de leer NxN líneas),
            # significa que ya no hay más bloques completos y termina la lectura
            if len(matriz) < n:
                break

            # Construir grafo dirigido
            G = nx.DiGraph()
            # Usamos identificadores únicos por índice, y guardamos la etiqueta en un atributo
            for idx, lbl in enumerate(labels):
                G.add_node(idx, label=lbl)

            # Añadir aristas
            for i in range(n):
                for j in range(n):
                    if i != j and matriz[i][j]:
                        G.add_edge(i, j)

            # Calcular grados e insertar atributo combinado
            for node in G.nodes:
                label = G.nodes[node]["label"]
                in_deg = G.in_degree[node]
                out_deg = G.out_degree[node]
                G.nodes[node]["tag"] = (label, in_deg, out_deg)

            # La primera matriz es intramolecular
            if not primer_matriz_leida:
                grafo_base = G.copy()
                primer_matriz_leida = True
                continue

            grafo_id += 1
            asignada = False
            # Comparamos con cada clase existente
            for j, clase in enumerate(clases):
                rep = clase["nx_graph"]
                if vf2pp_is_isomorphic(G, rep, node_label="tag"):
                    # Es isomorfo a la clase j
                    clase["count"] += 1
                    asignaciones.append((grafo_id, j + 1))
                    asignada = True
                    break

            if not asignada:
                # No es isomorfo a ninguna de las clases existentes
                # Guardamos la información
                clases.append({
                    "nx_graph": G,
                    "labels": labels,
                    "matriz": matriz,
                    "base_graph": grafo_base,
                    "count": 1
                })
                # La nueva clase tiene índice = len(clases)
                asignaciones.append((grafo_id, len(clases)))

    return clases, asignaciones

##############################################################################
#  Función principal (‘main’) con visualización de colores para cada tipo
##############################################################################
def main():
    # 1) Obtener argumentos y tipo de grafo
    if len(sys.argv) < 2:
        print("Uso: python tu_script.py fichero.txt [--directed]")
        print("Si no se especifica '--directed', se asumirá que el grafo es no dirigido.")
        return

    ruta_fichero = sys.argv[1]
    dirigido_flag = "--directed" in sys.argv[2:]

    # Basename del fichero de entrada (sin ruta ni extensión)
    basename = os.path.splitext(os.path.basename(ruta_fichero))[0]

    print(f"Leyendo grafos desde: {ruta_fichero}")
    print(f"Grafo dirigido: {dirigido_flag}")

    # 2) Clasificar en isomorfos on the fly 
    fn_stream = (
        clasificar_grafos_streaming_dir
        if dirigido_flag
        else clasificar_grafos_streaming_nodir
    )
    lista_grafos_unicos, asignaciones = fn_stream(ruta_fichero)
    total_count = sum(clase["count"] for clase in lista_grafos_unicos)

    # ============================
    # Reordenar clases por count (descendente)
    # ============================

    # Guardar el índice original de cada clase
    for idx, clase in enumerate(lista_grafos_unicos):
        clase["_old_index"] = idx + 1  # índices 1-based

    # Ordenar por count de mayor a menor
    lista_grafos_unicos.sort(key=lambda c: c["count"], reverse=True)

    # Construir un mapa old_index -> new_index
    old_to_new = {clase["_old_index"]: new_idx + 1 for new_idx, clase in enumerate(lista_grafos_unicos)}

    # Mantener el orden de aparición (id_grafo) y solo remapear la clase
    asignaciones = [(id_grafo, old_to_new[id_clase]) for (id_grafo, id_clase) in asignaciones]

    # Limpiar el campo auxiliar
    for clase in lista_grafos_unicos:
        del clase["_old_index"]


    # 3) Imprimir resumen por consola
    print("Resumen de grafos (clases de isomorfía) encontrados:")
    for i, clase in enumerate(lista_grafos_unicos, start=1):
        cociente = clase["count"] / total_count * 100
        print(f"  Clase #{i}:")
        print(f"    Labels: {clase['labels']}")
        print(f"    Matriz: {clase['matriz']}")
        print(f"    Count: {clase['count']}")
        print(f"    Population: {cociente:.4f}")
        print("")

    # 4) Escribir las asignaciones
    salida = f"{basename}_clasificacion.txt" 
    with open(salida, "w") as f_out:
        for (id_grafo, id_clase) in asignaciones:
            f_out.write(f"{id_grafo} {id_clase}\n")
    print(f"Se ha escrito la clasificación en '{salida}'.")

    # ============================
    # Visualización de cada clase
    # ============================
    
    # 1) Construir un mapeo global etiqueta → color
    all_labels = []
    for clase in lista_grafos_unicos:
        all_labels.extend(clase["labels"])
    etiquetas_unicas = sorted(set(all_labels))
    cmap = plt.get_cmap("tab20")  # colormap con hasta 20 colores
    etiqueta2color = {
        etiqueta: cmap(i / max(1, len(etiquetas_unicas) - 1))
        for i, etiqueta in enumerate(etiquetas_unicas)
    }
    
    # 2) Dibujar cada clase con colores según etiqueta
    for idx, clase in enumerate(lista_grafos_unicos, start=1):
        G = clase["nx_graph"]
        porcentaje = (clase["count"] / total_count) * 100

        # Grafo base (intramolecular) leído del primer bloque
        base = clase.get("base_graph", None)
        if base is None:
            raise ValueError(
                "No se encontró 'base_graph' en la clase. "
                "Asegúrate de guardar base_graph al leer la primera matriz."
            )
    
        # Calcular posiciones de nodos
        pos = nx.spring_layout(G, seed=idx)
    
        # Asignar color a cada nodo según su etiqueta
        node_colors = [etiqueta2color[G.nodes[n]["label"]] for n in G.nodes]
    
        # Dibujar nodos
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=300)

        # ============================
        # Separar enlaces intra/inter
        # ============================
        intramol_edges = []
        intermol_edges = []
        for u, v in G.edges():
            # En DiGraph la dirección importa: base.has_edge(u,v)
            # En Graph da igual el orden, has_edge(u,v) también funciona
            if base.has_edge(u, v):
                intramol_edges.append((u, v))
            else:
                intermol_edges.append((u, v))

        # Dibujar aristas intramoleculares en negro
        if dirigido_flag:
            nx.draw_networkx_edges(
                G, pos,
                edgelist=intramol_edges,
                edge_color="black",
                width=2.0,
                arrows=True,
                arrowstyle='-|>',
                arrowsize=12
            )
            # Dibujar aristas intermoleculares en gris
            nx.draw_networkx_edges(
                G, pos,
                edgelist=intermol_edges,
                edge_color="gray",
                style="dashed",
                width=1.5,
                arrows=True,
                arrowstyle='-|>',
                arrowsize=12
            )
        else:
            nx.draw_networkx_edges(
                G, pos,
                edgelist=intramol_edges,
                edge_color="black",
                width=2.0
            )
            nx.draw_networkx_edges(
                G, pos,
                edgelist=intermol_edges,
                edge_color="gray",
                style="dashed",
                width=1.5
            )

        # Etiquetas de nodos
        nx.draw_networkx_labels(G, pos)
    
        # Crear y mostrar leyenda (nodos por etiqueta + enlaces intra/inter)
        handles = []
        for etiqueta, color in etiqueta2color.items():
            handles.append(
                plt.Line2D([0], [0], marker="o", color="w",
                           markerfacecolor=color, markersize=8,
                           label=etiqueta)
            )

        # Handles extra para enlaces intra/inter (NetworkX no los añade a la leyenda)
#        handles.append(plt.Line2D([0], [0], color="black", lw=2.0))
#        handles.append(plt.Line2D([0], [0], color="gray", lw=1.5, linestyle="--"))

        plt.legend(handles=handles[::-1], title="Etiqueta") 
    
        # Título con el porcentaje (población)
        plt.title(f"Clase #{idx} — {porcentaje:.2f}%")
    
        # Guardar la figura
        filename = f"{basename}_conformer_{idx}.png"
        plt.savefig(filename, bbox_inches="tight")
        plt.clf()  # limpiar para la siguiente figura
    
    print("Se han generado las imágenes de cada clase.")

if __name__ == "__main__":
    main()

