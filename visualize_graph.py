"""
Sheet Metal Graph Visualizer
Reads a Salvagnini SMSerializer JSON file and renders the topological graph in 3D.

Usage:
    python visualize_graph.py [file] [--labels]   # CLI
    python launch.pyw                              # double-click app (no console)

    file     Path to the grafo.txt JSON file. If omitted, a welcome screen opens
             where you can drag a file onto the window or press O to browse.

    --labels Show vertex ID labels in the plot

Controls:
    Keys : 1=Top  2=Front  3=Side  4=Iso  5=Back  0=Home  R=Reset  S=Screenshot
           O = open file dialog
    Mouse: Left-drag=Orbit  Middle-drag=Pan  Scroll=Zoom
    Top-left checkboxes: toggle layer visibility
    Hover over geometry to inspect elements
    Drag a new graph file onto the window to reload
"""
from __future__ import annotations

import json
import math
import os
from collections import defaultdict
from dataclasses import dataclass
from typing import Optional

try:
    import numpy as np
    import pyvista as pv
except ImportError as exc:
    raise SystemExit(
        f"Missing dependency: {exc}\n"
        "Run:  pip install \"pyvista[all]\" numpy"
    )


# ---------------------------------------------------------------------------
# Win32 drag-and-drop helper (no extra dependencies – pure ctypes)
# ---------------------------------------------------------------------------

def _setup_win32_dnd(render_window, on_file_callback) -> bool:
    """
    Subclass the VTK render window's Win32 message handler so that dropping
    a file onto the window calls on_file_callback(path: str).
    Must be called after the window is realized (i.e. from a RenderEvent cb).
    Returns True on success, False on non-Windows or on error.
    """
    import sys
    if sys.platform != "win32":
        return False
    try:
        import ctypes

        WM_DROPFILES = 0x0233
        GWL_WNDPROC  = -4

        shell32 = ctypes.windll.shell32
        user32  = ctypes.windll.user32

        # VTK returns the HWND as an opaque string '_HEXADDR_p_void', not an int
        hwnd_raw = render_window.GetGenericWindowId()
        if hwnd_raw is None:
            return False
        if isinstance(hwnd_raw, str):
            parts = hwnd_raw.split("_")
            hwnd = int(parts[1], 16) if len(parts) >= 2 and parts[1] else 0
        else:
            hwnd = int(hwnd_raw)
        if not hwnd:
            return False

        WNDPROC = ctypes.WINFUNCTYPE(
            ctypes.c_ssize_t,
            ctypes.c_void_p, ctypes.c_uint,
            ctypes.c_size_t,  ctypes.c_ssize_t,
        )
        user32.GetWindowLongPtrW.restype = ctypes.c_ssize_t

        old_ptr = user32.GetWindowLongPtrW(hwnd, GWL_WNDPROC)

        shell32.DragAcceptFiles(hwnd, True)

        def _proc(h, msg, wp, lp):
            if msg == WM_DROPFILES:
                buf = ctypes.create_unicode_buffer(260)
                if shell32.DragQueryFileW(ctypes.c_void_p(wp), 0, buf, 260):
                    try:
                        on_file_callback(buf.value)
                    except Exception:
                        pass
                shell32.DragFinish(ctypes.c_void_p(wp))
                return 0
            return user32.CallWindowProcW(old_ptr, h, msg, wp, lp)

        cb = WNDPROC(_proc)
        render_window.__dnd_cb = cb   # keep a reference to prevent GC
        user32.SetWindowLongPtrW(hwnd, GWL_WNDPROC, cb)
        return True
    except Exception:
        return False


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
# Layer 2.5: Topology helper
# ---------------------------------------------------------------------------

def _extract_face_polygon(
    graph: SheetMetalGraph, face_id: int
) -> tuple[list[tuple[float, float, float]], bool]:
    """Trace boundary edges of a lamina face into an ordered vertex list."""
    edge_ids = graph.lamina_face_edges.get(face_id, [])
    if not edge_ids:
        return [], True
    adjacency: dict[int, list[int]] = defaultdict(list)
    for eid in edge_ids:
        v1, v2 = graph.edges.get(eid, (None, None))
        if v1 is not None and v2 is not None:
            adjacency[v1].append(v2)
            adjacency[v2].append(v1)
    if not adjacency:
        return [], True
    for neighbors in adjacency.values():
        if len(neighbors) != 2:
            return [], False
    start = next(iter(adjacency.keys()))
    polygon = [start]
    current = start
    prev = None
    while True:
        neighbors = adjacency[current]
        next_vert = neighbors[0] if neighbors[0] != prev else neighbors[1]
        if next_vert == start:
            return [graph.vertices[v] for v in polygon], True
        polygon.append(next_vert)
        prev = current
        current = next_vert
        if len(polygon) > len(adjacency) + 1:
            return [], False
    return [], False


# ---------------------------------------------------------------------------
# Layer 3: Renderer (PyVista / VTK – GPU-accelerated)
# ---------------------------------------------------------------------------

_WIN_W, _WIN_H = 1400, 900  # render window size in pixels


def visualize(
    graph: SheetMetalGraph,
    title: str = "Sheet Metal Graph",
    show_labels: bool = False,
) -> Optional[str]:
    """Render *graph* in a 3D window.  Returns a new file path if the user
    drops a file onto the window to reload, or None if they just closed it."""
    import vtk  # bundled with pyvista[all]

    _next_file: list[Optional[str]] = [None]

    pv.global_theme.background = "white"
    plotter = pv.Plotter(title=title, window_size=[_WIN_W, _WIN_H])
    plotter.enable_trackball_style()  # CAD-style: left-drag=orbit, middle=pan, scroll=zoom

    actors: dict[str, list] = {
        "faces": [], "edges": [], "bend_edges": [],
        "vertices": [], "labels": [], "bend_annots": [],
    }
    bad_faces: list[int] = []
    bend_eids = graph.bend_edge_ids

    # ------------------------------------------------------------------ #
    #  Faces                                                               #
    # ------------------------------------------------------------------ #
    face_pts_all: list = []
    face_cells_all: list[int] = []
    face_id_per_pt: list[int] = []
    face_id_per_cell: list[int] = []
    pt_offset = 0

    for face_id in graph.lamina_face_edges:
        poly, is_valid = _extract_face_polygon(graph, face_id)
        if is_valid and len(poly) >= 3:
            n = len(poly)
            face_pts_all.extend(poly)
            face_id_per_pt.extend([face_id] * n)
            face_cells_all.extend([n] + list(range(pt_offset, pt_offset + n)))
            face_id_per_cell.append(face_id)
            pt_offset += n
        elif not is_valid:
            bad_faces.append(face_id)

    if face_pts_all:
        face_mesh = pv.PolyData(
            np.array(face_pts_all, dtype=float),
            np.array(face_cells_all, dtype=int),
        )
        face_mesh.point_data["face_id"] = np.array(face_id_per_pt, dtype=int)
        # Per-cell scalar drives the color cycling (8 pastel shades)
        face_mesh.cell_data["color_key"] = np.array(
            [fid % 8 for fid in face_id_per_cell], dtype=float
        )
        actor = plotter.add_mesh(
            face_mesh,
            scalars="color_key",
            cmap="Pastel1",
            show_scalar_bar=False,
            opacity=0.75,
            show_edges=False,
            lighting=True,
            pickable=True,
            name="faces",
            clim=[0, 7],
        )
        actors["faces"].append(actor)

    # ------------------------------------------------------------------ #
    #  Edges                                                               #
    # ------------------------------------------------------------------ #
    reg_pts: list = []
    reg_lines: list[int] = []
    bend_pts: list = []
    bend_lines: list[int] = []

    for eid, (v1, v2) in graph.edges.items():
        if v1 not in graph.vertices or v2 not in graph.vertices:
            continue
        p1 = list(graph.vertices[v1])
        p2 = list(graph.vertices[v2])
        if eid in bend_eids:
            i = len(bend_pts)
            bend_pts.extend([p1, p2])
            bend_lines.extend([2, i, i + 1])
        else:
            i = len(reg_pts)
            reg_pts.extend([p1, p2])
            reg_lines.extend([2, i, i + 1])

    if reg_pts:
        mesh = pv.PolyData(np.array(reg_pts, dtype=float))
        mesh.lines = np.array(reg_lines, dtype=int)
        actor = plotter.add_mesh(
            mesh, color="#555555", line_width=1.5, pickable=False, name="edges"
        )
        actors["edges"].append(actor)

    if bend_pts:
        mesh = pv.PolyData(np.array(bend_pts, dtype=float))
        mesh.lines = np.array(bend_lines, dtype=int)
        actor = plotter.add_mesh(
            mesh, color="#E87722", line_width=3.0, pickable=False, name="bend_edges"
        )
        actors["bend_edges"].append(actor)

    # ------------------------------------------------------------------ #
    #  Vertices                                                            #
    # ------------------------------------------------------------------ #
    cluster_map = classify_vertices(graph.vertices)
    vid_list = sorted(graph.vertices.keys())
    coords = np.array([graph.vertices[v] for v in vid_list], dtype=float)

    pt_cloud = pv.PolyData(coords)
    pt_cloud.point_data["cluster"] = np.array(
        [cluster_map[v] for v in vid_list], dtype=float
    )
    pt_cloud.point_data["vertex_id"] = np.array(vid_list, dtype=int)

    actor = plotter.add_mesh(
        pt_cloud,
        scalars="cluster",
        cmap="tab10",
        clim=[0, 9],
        point_size=8,
        render_points_as_spheres=True,
        show_scalar_bar=False,
        pickable=True,
        name="vertices",
    )
    actors["vertices"].append(actor)

    # ------------------------------------------------------------------ #
    #  Vertex labels (optional)                                            #
    # ------------------------------------------------------------------ #
    if show_labels and len(vid_list) > 0:
        actor = plotter.add_point_labels(
            pt_cloud,
            [str(v) for v in vid_list],
            font_size=8,
            text_color="black",
            show_points=False,
            always_visible=False,
            name="labels",
        )
        actors["labels"].append(actor)

    # ------------------------------------------------------------------ #
    #  Bend annotations                                                    #
    # ------------------------------------------------------------------ #
    annot_pts: list = []
    annot_texts: list[str] = []
    for bid, angle_rad in graph.bend_angles.items():
        c = graph.bend_centroid(bid)
        if c:
            annot_pts.append(list(c))
            annot_texts.append(f"B{bid}  {math.degrees(angle_rad):.0f}°")

    if annot_pts:
        annot_cloud = pv.PolyData(np.array(annot_pts, dtype=float))
        actor = plotter.add_point_labels(
            annot_cloud,
            annot_texts,
            font_size=10,
            text_color="#C05000",
            bold=True,
            show_points=False,
            always_visible=True,
            name="bend_annots",
        )
        actors["bend_annots"].append(actor)

    # ------------------------------------------------------------------ #
    #  Layer toggle checkboxes (top-left, pixel coords from bottom-left)  #
    # ------------------------------------------------------------------ #
    TOGGLES = [
        ("faces",       "Faces",       True,        "#D4E8F0"),
        ("edges",       "Edges",       True,        "#555555"),
        ("bend_edges",  "Bend Edges",  True,        "#E87722"),
        ("vertices",    "Vertices",    True,        "#3366AA"),
        ("labels",      "Labels",      show_labels, "#333333"),
        ("bend_annots", "Bend Annots", True,        "#C05000"),
    ]

    CB_SIZE, CB_GAP = 25, 5
    _start_y = _WIN_H - 100  # leave room for title bar; 6 checkboxes × 30px = 180px used

    def _make_toggle(key: str):
        def _cb(state: bool):
            for a in actors[key]:
                try:
                    a.SetVisibility(1 if state else 0)
                except AttributeError:
                    try:
                        a.visibility = state
                    except AttributeError:
                        pass
            plotter.render()
        return _cb

    for idx, (layer_key, layer_label, initial, color) in enumerate(TOGGLES):
        y_px = _start_y - idx * (CB_SIZE + CB_GAP)
        plotter.add_checkbox_button_widget(
            _make_toggle(layer_key),
            value=initial,
            position=(10, y_px),
            size=CB_SIZE,
            border_size=2,
            color_on=color,
            color_off="#CCCCCC",
        )
        plotter.add_text(
            layer_label,
            position=(45, y_px + 4),   # pixel coords, same system as checkbox
            font_size=9,
            color="black",
        )

    # ------------------------------------------------------------------ #
    #  Hover inspection (mouse-move VTK observer)                         #
    # ------------------------------------------------------------------ #
    # Use a raw vtkTextActor so SetInput is guaranteed to work regardless
    # of which PyVista Actor wrapper add_text returns.
    _hover_text_actor = vtk.vtkTextActor()
    _hover_text_actor.SetInput("")
    _hover_text_actor.GetTextProperty().SetFontSize(13)
    _hover_text_actor.GetTextProperty().SetColor(0.05, 0.05, 0.05)
    _hover_text_actor.GetTextProperty().SetBackgroundColor(1.0, 1.0, 0.85)
    _hover_text_actor.GetTextProperty().SetBackgroundOpacity(0.85)
    _hover_text_actor.GetTextProperty().SetShadow(True)
    # Position: normalized viewport coords (0-1), lower-right
    _hover_text_actor.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
    _hover_text_actor.SetPosition(0.62, 0.02)
    plotter.renderer.AddActor2D(_hover_text_actor)

    # Picker created once — avoids per-frame allocation / GC churn
    _hover_picker = vtk.vtkPointPicker()
    _hover_picker.SetTolerance(0.01)
    _hover_renderer = plotter.renderer
    _hover_last: list[str] = [""]   # cache: only call SetInput when text changes
    _hover_next_t: list[float] = [0.0]  # 60 fps gate: earliest time for next pick
    _HOVER_INTERVAL = 1.0 / 60.0        # ~16.67 ms

    def _on_mouse_move(obj, event):
        import time
        now = time.perf_counter()
        if now < _hover_next_t[0]:
            return                       # drop frame — within the 60 fps window
        _hover_next_t[0] = now + _HOVER_INTERVAL

        try:
            x, y = plotter.iren.get_event_position()
        except AttributeError:
            x, y = plotter.iren.GetEventPosition()

        _hover_picker.Pick(x, y, 0, _hover_renderer)
        pid  = _hover_picker.GetPointId()
        dset = _hover_picker.GetDataSet()

        new_text = ""
        if dset is not None and pid >= 0:
            pd      = dset.GetPointData()
            vid_arr = pd.GetArray("vertex_id")
            fid_arr = pd.GetArray("face_id")

            if vid_arr is not None and pid < vid_arr.GetNumberOfTuples():
                vid = int(vid_arr.GetValue(pid))
                x_, y_, z_ = graph.vertices.get(vid, (0.0, 0.0, 0.0))
                new_text = (f"Vertex {vid}\n"
                            f"  x = {x_:.3f} mm\n"
                            f"  y = {y_:.3f} mm\n"
                            f"  z = {z_:.3f} mm")
            elif fid_arr is not None and pid < fid_arr.GetNumberOfTuples():
                fid    = int(fid_arr.GetValue(pid))
                nedges = len(graph.lamina_face_edges.get(fid, []))
                new_text = f"Lamina Face {fid}\n  {nedges} boundary edges"

        if new_text != _hover_last[0]:   # only update actor when text changes
            _hover_last[0] = new_text
            _hover_text_actor.SetInput(new_text)

    # Register observer on the VTK interactor
    try:
        plotter.iren.AddObserver("MouseMoveEvent", _on_mouse_move)
    except AttributeError:
        _raw = (getattr(plotter.iren, "_iren", None)
                or getattr(plotter.iren, "interactor", None))
        if _raw is not None:
            _raw.AddObserver("MouseMoveEvent", _on_mouse_move)

    # ------------------------------------------------------------------ #
    #  Drag-and-drop: drop a new file onto the window to reload           #
    # ------------------------------------------------------------------ #
    def _get_iren():
        iren = plotter.iren
        return getattr(iren, "interactor", None) or getattr(iren, "_iren", None) or iren

    _dnd_ready = [False]

    def _on_render_start(obj, event):
        if _dnd_ready[0]:
            return
        def _on_dropped(path):
            _next_file[0] = path
            _get_iren().TerminateApp()
        if _setup_win32_dnd(plotter.render_window, _on_dropped):
            _dnd_ready[0] = True

    plotter.render_window.AddObserver("StartEvent", _on_render_start)

    # ------------------------------------------------------------------ #
    #  View preset keyboard shortcuts                                      #
    # ------------------------------------------------------------------ #
    _VIEW_KEYS: dict[str, tuple] = {
        "1": ((0.0,  0.0,  1.0), (0.0, 1.0, 0.0)),  # Top
        "2": ((0.0, -1.0,  0.0), (0.0, 0.0, 1.0)),  # Front
        "3": ((1.0,  0.0,  0.0), (0.0, 0.0, 1.0)),  # Side
        "4": ((1.0, -1.0,  1.0), (0.0, 0.0, 1.0)),  # Iso
        "5": ((0.0,  1.0,  0.0), (0.0, 0.0, 1.0)),  # Back
        "0": ((1.0, -1.0,  1.0), (0.0, 0.0, 1.0)),  # Home
    }

    def _make_view_cb(pos_dir: tuple, up: tuple):
        def _cb():
            plotter.view_vector(pos_dir, viewup=up)
            plotter.reset_camera()
        return _cb

    for key, (pos_dir, up) in _VIEW_KEYS.items():
        plotter.add_key_event(key, _make_view_cb(pos_dir, up))

    def _reset_view():
        plotter.view_vector((1.0, -1.0, 1.0), viewup=(0.0, 0.0, 1.0))
        plotter.reset_camera()

    plotter.add_key_event("r", _reset_view)

    def _on_open_file():
        path = _open_file_dialog()
        if path:
            _next_file[0] = path
            _get_iren().TerminateApp()

    plotter.add_key_event("o", _on_open_file)
    plotter.add_key_event("O", _on_open_file)

    # ------------------------------------------------------------------ #
    #  Screenshot (S key → PNG timestamped in working directory)          #
    # ------------------------------------------------------------------ #
    import datetime

    def _save_screenshot():
        ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        fname = f"sheet_metal_{ts}.png"
        plotter.screenshot(fname, transparent_background=False)
        print(f"Screenshot saved: {fname}")

    plotter.add_key_event("s", _save_screenshot)

    # ------------------------------------------------------------------ #
    #  Axes indicator, title, stats                                        #
    # ------------------------------------------------------------------ #
    plotter.add_axes(line_width=2)

    nv  = len(graph.vertices)
    ne  = len(graph.edges)
    nb  = len(graph.bend_angles)
    nf  = len(graph.coarse_face_lfaces)
    nlf = len(graph.lamina_face_edges)
    plotter.add_text(
        f"{os.path.basename(title)}  |  {nv}V  {ne}E  {nb} bends  "
        f"{nf} faces ({nlf} laminas)",
        position="upper_edge",
        font_size=10,
        color="black",
    )

    if bad_faces:
        print(f"WARNING: {len(bad_faces)} malformed face(s) skipped: {bad_faces}")

    print("\n--- Controls ---")
    print("  Keys : 1=Top  2=Front  3=Side  4=Iso  5=Back  0=Home  R=Reset  S=Screenshot")
    print("         O=Open file")
    print("  Mouse: Left-drag=Orbit  Middle-drag=Pan  Scroll=Zoom")
    print("  Hover over geometry to inspect elements (vertex coords, face info)")
    print("  Top-left checkboxes: toggle layer visibility")
    print("  Drag a graph file onto the window to reload")

    plotter.show(auto_close=False)
    try:
        plotter.close()
    except Exception:
        pass
    return _next_file[0]


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

_REQUIRED_GRAPH_KEYS = [
    "VertexNodes", "EdgeToVertexArchs", "EdgeNodes",
    "LaminaFaceToEdgeArchs", "BendToLaminaFaceArchs",
    "BendNodes", "FaceToLaminaFaceArchs",
]


def _load_graph(path: str) -> SheetMetalGraph:
    try:
        with open(path, encoding="utf-8") as f:
            raw = json.load(f)
    except FileNotFoundError:
        raise SystemExit(f"ERROR: File not found: '{path}'")
    except json.JSONDecodeError as exc:
        raise SystemExit(f"ERROR: Invalid JSON in '{path}': {exc}")

    if "Graph" not in raw:
        raise SystemExit(
            f"ERROR: JSON root must contain a 'Graph' key.  "
            f"Keys found: {list(raw.keys())}"
        )
    graph_dict = raw["Graph"]
    missing = [k for k in _REQUIRED_GRAPH_KEYS if k not in graph_dict]
    if missing:
        raise SystemExit(
            f"ERROR: 'Graph' section is missing required keys: {missing}"
        )
    return SheetMetalGraph.from_dict(graph_dict)


def _open_file_dialog() -> Optional[str]:
    """Open a native OS file picker; return selected path or None."""
    try:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        path = filedialog.askopenfilename(
            title="Open Graph File",
            filetypes=[
                ("Graph files", "*.txt *.json"),
                ("All files", "*.*"),
            ],
        )
        root.destroy()
        return path.strip() or None
    except Exception:
        return None


def _welcome_screen() -> Optional[str]:
    """Splash window shown when no file is provided at startup.
    Returns the file path the user selected/dropped, or None if they quit."""
    pv.global_theme.background = "#1C1C2E"
    pl = pv.Plotter(
        title="Sheet Metal Graph Visualizer",
        window_size=[640, 400],
    )
    pl.enable_trackball_style()

    pl.add_text(
        "Sheet Metal Graph Visualizer",
        position="upper_edge",
        font_size=16,
        color="white",
    )
    pl.add_text(
        "Drop a graph file (.txt / .json) onto this window\n\n"
        "Press  O  to browse for a file\n\n"
        "Press  Q  to quit",
        position="lower_edge",
        font_size=12,
        color="#AAAAAA",
    )

    result: list[Optional[str]] = [None]

    def _get_iren():
        iren = pl.iren
        return getattr(iren, "interactor", None) or getattr(iren, "_iren", None) or iren

    _dnd_ready = [False]

    def _on_render_start(obj, event):
        if _dnd_ready[0]:
            return
        def _on_dropped(path):
            result[0] = path
            _get_iren().TerminateApp()
        if _setup_win32_dnd(pl.render_window, _on_dropped):
            _dnd_ready[0] = True

    pl.render_window.AddObserver("StartEvent", _on_render_start)

    def _on_open():
        path = _open_file_dialog()
        if path:
            result[0] = path
            _get_iren().TerminateApp()

    pl.add_key_event("o", _on_open)
    pl.add_key_event("O", _on_open)

    pl.show(auto_close=False)
    try:
        pl.close()
    except Exception:
        pass

    pv.global_theme.background = "white"
    return result[0]


def main() -> None:
    import sys

    show_labels = "--labels" in sys.argv
    cli_files = [a for a in sys.argv[1:] if not a.startswith("-")]

    path: Optional[str] = cli_files[0] if cli_files else _welcome_screen()

    while path:
        try:
            graph = _load_graph(path)
        except SystemExit as exc:
            print(exc)
            path = _welcome_screen()
            continue
        path = visualize(graph, title=path, show_labels=show_labels)


if __name__ == "__main__":
    main()
