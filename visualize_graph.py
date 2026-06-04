"""
Sheet Metal Graph Visualizer
Reads a Salvagnini SMSerializer JSON file and renders the topological graph in 3D.

Usage:
    python visualize_graph.py [file] [--labels]

    file     Path to the grafo.txt JSON file (default: grafo.txt)
    --labels Show vertex ID labels in the plot
"""
from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from typing import Optional

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


# ---------------------------------------------------------------------------
# Layer 1: Data model & parser
# ---------------------------------------------------------------------------

def _decode_arch(arch: list[list]) -> dict[int, list[int]]:
    """Decode two-parallel-array arch encoding into {source_id: [target_ids]}."""
    result: dict[int, list[int]] = defaultdict(list)
    for src, tgt in zip(arch[0], arch[1]):
        result[int(src)].append(int(tgt))
    return dict(result)


@dataclass
class SheetMetalGraph:
    vertices: dict[int, tuple[float, float, float]]   # {id: (x, y, z)}
    edges: dict[int, tuple[int, int]]                  # {id: (v1, v2)}
    edge_weights: dict[int, float]                     # {id: weight}
    lamina_face_edges: dict[int, list[int]]            # {face_id: [edge_ids]}
    bend_lamina_faces: dict[int, list[int]]            # {bend_id: [face_ids]}
    bend_angles: dict[int, float]                      # {bend_id: angle_rad}
    coarse_face_lfaces: dict[int, list[int]]           # {coarse_face_id: [lamina_face_ids]}

    # --- static parsers ---

    @staticmethod
    def _parse_vertices(nodes: list[list]) -> dict[int, tuple[float, float, float]]:
        return {int(n[0]): (float(n[1]), float(n[2]), float(n[3])) for n in nodes}

    @staticmethod
    def _parse_edges(arch: list[list]) -> dict[int, tuple[int, int]]:
        """Each edge appears twice in the arch; collect both vertex endpoints."""
        acc: dict[int, list[int]] = defaultdict(list)
        for eid, vid in zip(arch[0], arch[1]):
            acc[int(eid)].append(int(vid))
        return {
            eid: (verts[0], verts[1] if len(verts) > 1 else verts[0])
            for eid, verts in acc.items()
        }

    @staticmethod
    def _parse_edge_weights(nodes: list[list]) -> dict[int, float]:
        return {i: float(n[1]) for i, n in enumerate(nodes)}

    @staticmethod
    def _parse_bend_angles(nodes: list[list]) -> dict[int, float]:
        # BendNode layout: [id, ?, ?, ?, angle_rad, ?, ?]
        return {int(n[0]): float(n[4]) for n in nodes}

    # --- factory ---

    @classmethod
    def from_dict(cls, graph: dict) -> "SheetMetalGraph":
        # BendToLaminaFaceArchs uses sequential positions (0,1,2…) not actual bend IDs.
        # Build the positional → actual-ID mapping from BendNodes order.
        bend_seq_to_id = {i: int(n[0]) for i, n in enumerate(graph["BendNodes"])}
        raw_bend_lfaces = _decode_arch(graph["BendToLaminaFaceArchs"])
        bend_lamina_faces = {
            bend_seq_to_id[seq]: faces
            for seq, faces in raw_bend_lfaces.items()
            if seq in bend_seq_to_id
        }
        return cls(
            vertices=cls._parse_vertices(graph["VertexNodes"]),
            edges=cls._parse_edges(graph["EdgeToVertexArchs"]),
            edge_weights=cls._parse_edge_weights(graph["EdgeNodes"]),
            lamina_face_edges=_decode_arch(graph["LaminaFaceToEdgeArchs"]),
            bend_lamina_faces=bend_lamina_faces,
            bend_angles=cls._parse_bend_angles(graph["BendNodes"]),
            coarse_face_lfaces=_decode_arch(graph["FaceToLaminaFaceArchs"]),
        )

    # --- derived accessors ---

    @property
    def bend_edge_ids(self) -> set[int]:
        """Edge IDs that bound any lamina face adjacent to a bend."""
        result: set[int] = set()
        for face_ids in self.bend_lamina_faces.values():
            for fid in face_ids:
                result.update(self.lamina_face_edges.get(fid, []))
        return result

    def bend_centroid(self, bend_id: int) -> Optional[tuple[float, float, float]]:
        """Average position of all vertices in faces adjacent to this bend."""
        face_ids = self.bend_lamina_faces.get(bend_id, [])
        vids: set[int] = set()
        for fid in face_ids:
            for eid in self.lamina_face_edges.get(fid, []):
                v1, v2 = self.edges.get(eid, (None, None))
                if v1 is not None:
                    vids.add(v1)
                    vids.add(v2)
        valid = [self.vertices[v] for v in vids if v in self.vertices]
        if not valid:
            return None
        xs, ys, zs = zip(*valid)
        return (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))


# ---------------------------------------------------------------------------
# Layer 2: Vertex clustering by z-value
# ---------------------------------------------------------------------------

def classify_vertices(
    vertices: dict[int, tuple[float, float, float]],
    n_clusters: int = 4,
) -> dict[int, int]:
    """
    Assign each vertex a cluster index (0 … n_clusters-1) based on z percentile.
    Works for any number of vertices and any z distribution.
    """
    if not vertices:
        return {}
    sorted_z = sorted((z, vid) for vid, (_, _, z) in vertices.items())
    n = len(sorted_z)
    breaks = [sorted_z[min(int(n * (i + 1) / n_clusters), n - 1)][0]
              for i in range(n_clusters - 1)]
    result: dict[int, int] = {}
    for z, vid in sorted_z:
        c = sum(z > b for b in breaks)
        result[vid] = c
    return result


# ---------------------------------------------------------------------------
# Layer 3: Renderer
# ---------------------------------------------------------------------------

_CLUSTER_COLORS = ["#4878CF", "#6ACC65", "#D65F5F", "#B47CC7",
                   "#E78AC3", "#A6D854", "#FFD92F", "#E5C494"]
_CLUSTER_LABELS = ["z-cluster-0", "z-cluster-1", "z-cluster-2", "z-cluster-3",
                   "z-cluster-4", "z-cluster-5", "z-cluster-6", "z-cluster-7"]


def visualize(
    graph: SheetMetalGraph,
    title: str = "Sheet Metal Graph",
    show_labels: bool = False,
) -> None:
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection="3d")

    # Store original and current vertex positions
    original_vertices = dict(graph.vertices)
    current_vertices = dict(graph.vertices)
    bad_faces = []

    # Pre-compute axis bounds
    xs, ys, zs = zip(*original_vertices.values())
    ranges = [max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs)]
    max_range = max(ranges) / 2
    cx, cy, cz = (min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2, (min(zs) + max(zs)) / 2

    zoom_level = [1.0]  # Default zoom level

    def set_view_limits(zoom=None):
        """Set axis limits for uniform scaling."""
        if zoom is None:
            zoom = zoom_level[0]
        range_scaled = max_range / zoom
        ax.set_xlim([cx - range_scaled, cx + range_scaled])
        ax.set_ylim([cy - range_scaled, cy + range_scaled])
        ax.set_zlim([cz - range_scaled, cz + range_scaled])
        ax.set_box_aspect([1, 1, 1])

    def extract_face_polygon(face_id: int) -> tuple[list, bool]:
        """Extract vertices of a lamina face by tracing edges. Returns (polygon, is_valid)."""
        edge_ids = graph.lamina_face_edges.get(face_id, [])
        if not edge_ids:
            return [], True

        # Build adjacency map from edges
        adjacency: dict[int, list[int]] = defaultdict(list)
        for eid in edge_ids:
            v1, v2 = graph.edges.get(eid, (None, None))
            if v1 is not None and v2 is not None:
                adjacency[v1].append(v2)
                adjacency[v2].append(v1)

        if not adjacency:
            return [], True

        # Check for issues: vertices with wrong degree
        for vid, neighbors in adjacency.items():
            if len(neighbors) != 2:
                return [], False  # Invalid: vertex should have exactly 2 edges

        # Trace the polygon by following edges
        start = next(iter(adjacency.keys()))
        polygon = [start]
        current = start
        prev = None

        while True:
            neighbors = adjacency[current]
            next_vert = neighbors[0] if neighbors[0] != prev else neighbors[1]

            if next_vert == start:
                # Closed loop - valid polygon
                return [current_vertices[v] for v in polygon], True

            polygon.append(next_vert)
            prev = current
            current = next_vert

            if len(polygon) > len(adjacency) + 1:
                # Too many vertices - probably a self-intersecting or malformed face
                return [], False

        return [], False

    def redraw_plot():
        """Redraw the entire plot with current vertex positions."""
        nonlocal bad_faces
        ax.clear()
        bend_eids = graph.bend_edge_ids
        legend_handles_local: list = []
        bad_faces = []

        # --- Render solid lamina faces ---
        face_colors = ["#E8F4F8", "#D4E8F0", "#F0E8F8", "#F8F0E8",
                       "#E8F8E8", "#F8E8E8", "#E8E8F8", "#F8F8E8"]
        for face_id in graph.lamina_face_edges.keys():
            poly, is_valid = extract_face_polygon(face_id)
            if is_valid and len(poly) >= 3:
                face_color = face_colors[face_id % len(face_colors)]
                poly_collection = Poly3DCollection([poly], alpha=0.7, facecolor=face_color,
                                                    edgecolor="#333333", linewidth=0.5)
                ax.add_collection3d(poly_collection)
            elif not is_valid:
                bad_faces.append(face_id)

        # --- Highlight bend edges (batch render) ---
        bend_edge_segs = []
        for eid in bend_eids:
            if eid not in graph.edges:
                continue
            v1, v2 = graph.edges[eid]
            if v1 not in current_vertices or v2 not in current_vertices:
                continue
            x1, y1, z1 = current_vertices[v1]
            x2, y2, z2 = current_vertices[v2]
            bend_edge_segs.append([[x1, x2], [y1, y2], [z1, z2]])

        for seg in bend_edge_segs:
            ax.plot(seg[0], seg[1], seg[2], color="#E87722", linewidth=3, zorder=5)

        if bend_eids:
            legend_handles_local.append(
                mlines.Line2D([], [], color="#E87722", linewidth=3,
                              label=f"Bend edge ({len(graph.bend_angles)} bends)")
            )

        # --- Vertices (single scatter call for efficiency) ---
        xs, ys, zs = zip(*current_vertices.values())
        ax.scatter(xs, ys, zs, c="#333333", s=8, zorder=10, depthshade=False, alpha=0.6)

        # --- Vertex ID labels (optional) ---
        if show_labels:
            for vid, (x, y, z) in current_vertices.items():
                ax.text(x, y, z, str(vid), fontsize=5, color="#333333", zorder=15)

        # --- Bend annotations ---
        for bid, angle_rad in graph.bend_angles.items():
            c = graph.bend_centroid(bid)
            if c:
                angle_deg = math.degrees(angle_rad)
                ax.text(c[0], c[1], c[2], f" B{bid}\n {angle_deg:.0f}°",
                        fontsize=7, color="#C05000", fontweight="bold", zorder=20)

        # --- Set axes ---
        set_view_limits()
        ax.set_xlabel("X (mm)")
        ax.set_ylabel("Y (mm)")
        ax.set_zlabel("Z (mm)")
        nv = len(graph.vertices)
        ne = len(graph.edges)
        nb = len(graph.bend_angles)
        nf = len(graph.coarse_face_lfaces)
        nlf = len(graph.lamina_face_edges)
        ax.set_title(
            f"{title}\n"
            f"{nv} vertices · {ne} edges · {nb} bends · {nf} faces ({nlf} lamina faces)",
            fontsize=9,
        )
        ax.legend(handles=legend_handles_local, loc="upper left", fontsize=7, framealpha=0.75)
        fig.canvas.draw_idle()

    # View preset functions
    def set_view(elev, azim, label=""):
        """Set view angle."""
        ax.view_init(elev=elev, azim=azim)
        fig.canvas.draw_idle()

    # Keyboard shortcuts for view angles
    def on_key(event):
        if event.key == '1':  # Top view
            set_view(90, 0, "Top")
        elif event.key == '2':  # Front view
            set_view(0, 0, "Front")
        elif event.key == '3':  # Side view
            set_view(0, 90, "Side")
        elif event.key == '4':  # Isometric view
            set_view(20, 45, "Isometric")
        elif event.key == '5':  # Back view
            set_view(0, 180, "Back")
        elif event.key == '0':  # Home view (default isometric)
            set_view(20, 45, "Home")
        elif event.key == 'r':  # Reset zoom
            zoom_level[0] = 1.0
            slider_zoom.set_val(1.0)
            ax.view_init(elev=20, azim=45)
            set_view_limits(1.0)
            fig.canvas.draw_idle()

    fig.canvas.mpl_connect('key_press_event', on_key)

    # Add view preset buttons
    view_buttons = [
        ([0.05, 0.05, 0.06, 0.04], "Top", (90, 0)),
        ([0.12, 0.05, 0.06, 0.04], "Front", (0, 0)),
        ([0.19, 0.05, 0.06, 0.04], "Side", (0, 90)),
        ([0.26, 0.05, 0.06, 0.04], "Iso", (20, 45)),
        ([0.33, 0.05, 0.06, 0.04], "Back", (0, 180)),
    ]

    for pos, label, view_angles in view_buttons:
        ax_btn = fig.add_axes(pos)

        def make_callback(angles):
            def on_click(event):
                set_view(angles[0], angles[1], label)
            return on_click

        btn = Button(ax_btn, label, color="#D0D0D0", hovercolor="#E0E0E0")
        btn.on_clicked(make_callback(view_angles))

    # Add zoom slider
    ax_zoom = fig.add_axes([0.45, 0.05, 0.25, 0.04])
    slider_zoom = Slider(ax_zoom, "Zoom", 0.1, 5.0, valinit=1.0, color="#6699DD")

    def on_zoom(val):
        zoom_level[0] = val
        set_view_limits(val)
        fig.canvas.draw_idle()

    slider_zoom.on_changed(on_zoom)

    # Initial draw
    redraw_plot()

    # Add warning if bad faces detected
    if bad_faces:
        print(f"WARNING: {len(bad_faces)} malformed face(s) detected in graph input: {bad_faces}")
        print("   These faces have vertices with incorrect edge connectivity (not exactly 2 edges per vertex).")
        print("   This is likely a graph input issue from the SMSerializer.")

    print("\n--- Controls ---")
    print("  Keyboard: 1=Top, 2=Front, 3=Side, 4=Isometric, 5=Back, 0=Home, R=Reset zoom")
    print("  Mouse: Drag to rotate, Scroll to zoom, Middle-click to pan")

    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Visualize a Salvagnini SMSerializer sheet metal graph in 3D."
    )
    parser.add_argument(
        "file", nargs="?", default="grafo.txt",
        help="Path to the JSON graph file (default: grafo.txt)"
    )
    parser.add_argument(
        "--labels", action="store_true",
        help="Annotate each vertex with its ID"
    )
    args = parser.parse_args()

    with open(args.file, encoding="utf-8") as f:
        data = json.load(f)

    graph = SheetMetalGraph.from_dict(data["Graph"])
    visualize(graph, title=args.file, show_labels=args.labels)
