import numpy as np
import networkx as nx
from itertools import combinations
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# -----------------------------
# Build NetworkX graphs for loops
# -----------------------------

def build_loop_graphs_from_loops_dict(loops, close=True, directed=False):
    """
    loops: dict {loopID: [ {"nodeID": int, "P": np.array([x,y,z]), ...}, ... ] }
    Returns:
      loop_graphs: dict loopID -> nx.Graph/DiGraph
      G_all: combined graph with edge attr loopID
      pos: dict nodeID -> position (P), for drawing layouts if you want
    """
    loop_graphs = {}
    pos = {}  # nodeID -> P
    G_all = nx.DiGraph() if directed else nx.Graph()

    for lid, nodes in loops.items():
        if len(nodes) < 2:
            continue
        G = nx.DiGraph() if directed else nx.Graph()
        ids = [int(n["nodeID"]) for n in nodes]
        Ps  = [np.array(n["P"], dtype=float) for n in nodes]

        # store positions (last one wins if repeated across loops)
        for nid, P in zip(ids, Ps):
            pos[nid] = P
            G.add_node(nid, P=P, loopID=int(lid))
            G_all.add_node(nid, P=P)

        # connect consecutive nodes
        for a, b in zip(ids[:-1], ids[1:]):
            G.add_edge(a, b, loopID=int(lid))
            G_all.add_edge(a, b, loopID=int(lid))

        # close the loop
        if close and ids[0] != ids[-1]:
            G.add_edge(ids[-1], ids[0], loopID=int(lid))
            G_all.add_edge(ids[-1], ids[0], loopID=int(lid))

        loop_graphs[int(lid)] = G
    return loop_graphs, G_all, pos


def plot_all_loops_connectivity_2d(G_loops_all, seed=0, node_size=60, font_size=6, show_labels=False):
    plt.figure(figsize=(10, 8))
    pos2d = nx.spring_layout(G_loops_all, seed=seed)

    # Collect loopIDs for edges
    e_loop = []
    for u, v, d in G_loops_all.edges(data=True):
        e_loop.append(d.get("loopID", -1))
    unique = sorted(set(e_loop))
    lid_to_i = {lid:i for i, lid in enumerate(unique)}
    cmap = cm.get_cmap("tab20", max(1, len(unique)))
    edge_colors = [cmap(lid_to_i[lid]) for lid in e_loop]

    # Draw
    nx.draw_networkx_nodes(G_loops_all, pos2d, node_size=node_size, node_color=edge_colors)
    nx.draw_networkx_edges(G_loops_all, pos2d, width=1.5, edge_color=edge_colors)
    if show_labels:
        nx.draw_networkx_labels(G_loops_all, pos2d, font_size=font_size)

    plt.title("All loops-nodes connectivity")
    plt.axis("off")
    plt.tight_layout()
    plt.show()



# -----------------------------
# Build NetworkX graphs for junctions
# -----------------------------

def loop_centroids_2d(loops, plane="xy"):
    """
    loops: dict loopID -> list of dicts with key 'P' (3D)
    plane: 'xy' / 'xz' / 'yz'
    Returns: pos2d dict loopID -> (x,y) for NetworkX
    """
    pos = {}
    for lid, nodes in loops.items():
        P = np.array([n["P"] for n in nodes], float)
        if P.size == 0:
            continue
        c = P.mean(axis=0)
        if plane == "xy":
            pos[int(lid)] = (float(c[0]), float(c[1]))
        elif plane == "xz":
            pos[int(lid)] = (float(c[0]), float(c[2]))
        elif plane == "yz":
            pos[int(lid)] = (float(c[1]), float(c[2]))
        else:
            raise ValueError("plane must be 'xy', 'xz', or 'yz'")
    return pos


def build_loop_junction_graph(junctions, all_loopIDs):
    """
    Nodes = all_loopIDs (even isolated)
    Edges = junction connections inferred from junctions
    """
    G = nx.Graph()
    for lid in all_loopIDs:
        G.add_node(int(lid))

    pair_to_links = defaultdict(list)
    for j in junctions:
        lids = sorted(set(int(x) for x in j["loopIDs"]))
        # only keep loopIDs that exist in all_loopIDs
        lids = [lid for lid in lids if lid in set(all_loopIDs)]
        if len(lids) < 2:
            continue
        for a, b in combinations(lids, 2):
            pair_to_links[(a, b)].append(j["link_key"])

    for (a, b), link_keys in pair_to_links.items():
        G.add_edge(a, b, count=len(link_keys), link_keys=link_keys)

    return G


def plot_loops_and_junctions_2d(G, pos, node_size=450, font_size=9, show_edge_counts=True):
    """
    G: loop-junction graph (nodes=loopIDs)
    pos: dict loopID -> (x,y) (centroid-projected positions)
    """
    if G.number_of_nodes() == 0:
        print("No loop-loop junction connections to plot (graph is empty).")
        return

    # Ensure all nodes have a position; fallback to spring layout for missing ones
    missing = [n for n in G.nodes() if n not in pos]
    if len(missing) > 0:
        pos_fallback = nx.spring_layout(G.subgraph(missing), seed=0)
        pos = {**pos, **pos_fallback}

    loop_ids = sorted(G.nodes())
    cmap = cm.get_cmap("tab20", max(1, len(loop_ids)))
    color_map = {lid: cmap(i) for i, lid in enumerate(loop_ids)}
    node_colors = [color_map[lid] for lid in loop_ids]

    plt.figure(figsize=(10, 8))
    nx.draw_networkx_edges(G, pos, edge_color="k", width=2.0, alpha=0.8)     # black edges = junctions
    nx.draw_networkx_nodes(G, pos, nodelist=loop_ids, node_color=node_colors, node_size=node_size)     # colored nodes = loops
    nx.draw_networkx_labels(G, pos, labels={lid: str(lid) for lid in loop_ids}, font_size=font_size)     # labels = loopID

    if show_edge_counts:
        edge_labels = {(u, v): d.get("count", 1) for u, v, d in G.edges(data=True)}
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=8)     # optional edge labels = number of junction segments

    plt.title("Loop-junction connections in 2D")
    plt.axis("equal")
    plt.axis("off")
    plt.tight_layout()
    plt.show()
