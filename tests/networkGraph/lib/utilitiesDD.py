import numpy as np
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import matplotlib.cm as cm

# -----------------------------
# Segment Information Extraction
# -----------------------------
def build_segments(DN, skip_boundary=True, skip_zeroB=True):
    """
    Returns list of dicts:
      {
        "link_key": lk,
        "source_id": sid, "sink_id": tid,
        ...
        ...
        "loopIDs": set(...)
      }
    """
    out = []
    skipped = {"null": 0, "zeroB": 0, "boundary": 0}
    for lk in DN.networkLinks().keys():
        link = DN.networkLinks().getRef(lk)
        if link is None:
            skipped["null"] += 1
            continue
        if skip_zeroB and link.hasZeroBurgers():
            skipped["zeroB"] += 1
            continue
        if skip_boundary and (link.isBoundarySegment() or link.isGrainBoundarySegment()):
            skipped["boundary"] += 1
            continue
        sid = int(lk[0]); tid = int(lk[1])
        sid_position = DN.networkNodes().getRef(sid).position()
        tid_position = DN.networkNodes().getRef(tid).position()
        sid_isBN = DN.networkNodes().getRef(sid).isBoundaryNode()
        tid_isBN = DN.networkNodes().getRef(tid).isBoundaryNode()
        loop_ids = link.loopIDs() 
        out.append({
            "link_key": lk,
            "source_id": sid, "sink_id": tid,
            "source_pos": sid_position, "sink_pos": tid_position,
            "source_is_boundary": sid_isBN, "sink_is_boundary": tid_isBN,
            "burgers": np.array(link.burgers(), dtype=float),
            "glidePlaneNormal": np.array(link.glidePlaneNormal(), dtype=float),
            "is_boundary_link": bool(link.isBoundarySegment()),
            "is_glissile": bool(link.isGlissile()),
            "is_sessile": bool(link.isSessile()),
            "loopIDs": set(loop_ids) if loop_ids is not None else set(),
        })
    return out, skipped


# -----------------------------
# Loop Information Extraction
# -----------------------------
def build_loop_nodes(ddio):
    loop_nodes = defaultdict(dict)
    for ln in ddio.loopNodes:  # DislocationLoopNodeIO
        lid = int(ln.loopID)
        nid = int(ln.sID)
        loop_nodes[lid][nid] = {
            "networkNodeID": int(ln.networkNodeID),
            "P": np.array(ln.P, dtype=float),
            "periodicShift": np.array(ln.periodicShift, dtype=float),
        }
    return loop_nodes

def build_connectivity(ddio):
    nxt = defaultdict(dict)
    for ll in ddio.loopLinks:
        lid = int(ll.loopID)
        s = int(ll.sourceID)
        t = int(ll.sinkID)
        nxt[lid][s] = t
    return nxt


def order_loop_nodes(loop_nodes_dict, next_map_for_loop):
    """
    Given loop_nodes_dict: {nodeID: nodeinfo}, and next_map_for_loop: {sourceID: sinkID},
    returns an ordered list of nodeIDs following the loopLinks connectivity.
    Falls back to sorted nodeIDs if connectivity is incomplete.
    """
    if not next_map_for_loop:
        return sorted(loop_nodes_dict.keys())
    # pick a deterministic start node
    start = min(loop_nodes_dict.keys())
    ordered = [start]
    visited = {start}
    cur = start
    for _ in range(len(loop_nodes_dict) + 2):
        if cur not in next_map_for_loop:
            break
        cur = next_map_for_loop[cur]
        if cur in visited:
            break
        ordered.append(cur)
        visited.add(cur)
    # if we didn't visit them all, fall back to sorted (or append remaining)
    if len(visited) != len(loop_nodes_dict):
        remaining = [n for n in sorted(loop_nodes_dict.keys()) if n not in visited]
        ordered.extend(remaining)
    return ordered


def build_loops(ddio, filter_loopIDs=None):
    loop_nodes = build_loop_nodes(ddio)
    nxt = build_connectivity(ddio)
    loops = {}
    for lid, nodes_dict in loop_nodes.items():
        if filter_loopIDs is not None and lid not in filter_loopIDs:
            continue
        order = order_loop_nodes(nodes_dict, nxt.get(lid, {}))
        loops[lid] = [{"nodeID": nid, **nodes_dict[nid]} for nid in order]
    return loops


# -----------------------------
# Junctions Detection
# -----------------------------
def find_junction_network_links(DN, skip_zeroB=False, skip_boundary=True):
    """
    Returns:
      junctions: list of dicts for each junction NetworkLink
      loop_to_junction_links: loopID -> list of (link_key)
    """
    junctions = []
    loop_to_junction_links = defaultdict(list)
    for lk in DN.networkLinks().keys():
        link = DN.networkLinks().getRef(lk)
        if link is None:
            continue
        if skip_zeroB and link.hasZeroBurgers():
            continue
        if skip_boundary and (link.isBoundarySegment() or link.isGrainBoundarySegment()):
            continue

        lids = list(link.loopIDs())  # bound
        if len(lids) >= 2:
            # endpoints in NetworkNode IDs
            sid = int(lk[0])
            tid = int(lk[1])

            junctions.append({
                "link_key": lk,
                "source_nID": sid,
                "sink_nID": tid,
                "chord": np.array(link.chord(), float),
                "burgers": np.array(link.burgers(), float),
                "loopIDs": [int(x) for x in lids],
            })
            for lid in lids:
                loop_to_junction_links[int(lid)].append(lk)
    return junctions, loop_to_junction_links



# -----------------------------
# Plotting Stuff
# -----------------------------
def draw_box(ax, xmin, xmax, lw=2.0, alpha=0.4):
    X0, X1 = np.array(xmin, float), np.array(xmax, float)
    xs = [X0[0], X1[0]]
    ys = [X0[1], X1[1]]
    zs = [X0[2], X1[2]]
    edges = [
        ([xs[0], xs[1]], [ys[0], ys[0]], [zs[0], zs[0]]),
        ([xs[0], xs[1]], [ys[1], ys[1]], [zs[0], zs[0]]),
        ([xs[0], xs[1]], [ys[0], ys[0]], [zs[1], zs[1]]),
        ([xs[0], xs[1]], [ys[1], ys[1]], [zs[1], zs[1]]),

        ([xs[0], xs[0]], [ys[0], ys[1]], [zs[0], zs[0]]),
        ([xs[1], xs[1]], [ys[0], ys[1]], [zs[0], zs[0]]),
        ([xs[0], xs[0]], [ys[0], ys[1]], [zs[1], zs[1]]),
        ([xs[1], xs[1]], [ys[0], ys[1]], [zs[1], zs[1]]),

        ([xs[0], xs[0]], [ys[0], ys[0]], [zs[0], zs[1]]),
        ([xs[1], xs[1]], [ys[0], ys[0]], [zs[0], zs[1]]),
        ([xs[0], xs[0]], [ys[1], ys[1]], [zs[0], zs[1]]),
        ([xs[1], xs[1]], [ys[1], ys[1]], [zs[0], zs[1]]),
    ]
    for ex, ey, ez in edges:
        ax.plot(ex, ey, ez, linewidth=lw, alpha=alpha, color="k")


def plot_loops_3d(loops, xMin, xMax, outpath=None, title=None, lw=2.0, ms=2.5, alpha=1.0, plotBox=False):
    """
    loops: dict {loopID: [ {nodeID:..., "P": np.array([x,y,z]), ...}, ... ] }
    """
    fig = plt.figure(figsize=(7, 12))
    ax = fig.add_subplot(111, projection="3d")
    loop_ids = sorted(list(loops.keys()))
    if len(loop_ids) == 0:
        print("No loops to plot.")
        return
    cmap = cm.get_cmap("tab20", max(1, len(loop_ids)))  
    for i, lid in enumerate(loop_ids):
        nodes = loops[lid]
        P = np.array([nd["P"] for nd in nodes], dtype=float)  # (N,3)
        # close the loop visually
        if len(P) > 1:
            Pclosed = np.vstack([P, P[0]])
        else:
            Pclosed = P
        ax.plot(Pclosed[:,0], Pclosed[:,1], Pclosed[:,2],
                lw=lw, alpha=alpha, color=cmap(i), label=f"loop {lid}")
        # show nodes as points
        ax.scatter(P[:,0], P[:,1], P[:,2], s=ms, alpha=alpha, color=cmap(i))

    if plotBox:
        draw_box(ax, xMin, xMax, lw=1.5, alpha=0.45)
        ax.set_box_aspect((xMax - xMin))
    ax.set_axis_off()

    if title is not None:
        ax.set_title(title)

    plt.tight_layout()
    if outpath is not None:
        plt.savefig(outpath, dpi=200, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()