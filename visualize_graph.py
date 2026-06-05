"""
Sheet Metal Graph Visualizer
Reads a Salvagnini SMSerializer JSON file and renders the topological graph in 3D.

Usage:
    python visualize_graph.py [file_or_folder] [--labels]   # CLI
    python launch.pyw                                        # double-click app (no console)

    file_or_folder  Path to a .graph/.txt/.json file, or a folder containing .graph
                    files.  If omitted, a welcome screen opens where you can drag a
                    file or folder onto the window or press O to browse.

    --labels  Show vertex ID labels in the plot

Controls:
    Keys : 1=Top  2=Front  3=Side  4=Iso  5=Back  0=Home  R=Reset  S=Screenshot
           O=Open file  F=Open folder  +/-=Zoom
           Up/Down = navigate playlist (when a folder is open)
    Mouse: Left-drag=Orbit  Middle-drag=Pan  Scroll=Zoom
    Top-left checkboxes: toggle layer visibility
    Hover over geometry to inspect elements
    Drag a file or folder onto the window to reload / open playlist
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
# Debug logger — writes to debug.log (visible even in .pyw / no-console mode)
# ---------------------------------------------------------------------------

def _dbg(msg: str) -> None:
    import datetime
    _path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "debug.log")
    try:
        with open(_path, "a", encoding="utf-8") as _f:
            _f.write(f"[{datetime.datetime.now():%H:%M:%S}] {msg}\n")
    except Exception:
        pass


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
            ctypes.c_ssize_t,                           # LRESULT
            ctypes.c_void_p,                            # HWND
            ctypes.c_uint,                              # UINT  msg
            ctypes.c_size_t,                            # WPARAM (unsigned ptr-size)
            ctypes.c_ssize_t,                           # LPARAM (signed ptr-size)
        )

        # Set argtypes so ctypes uses 64-bit integers on 64-bit Windows
        user32.GetWindowLongPtrW.restype  = ctypes.c_ssize_t
        user32.GetWindowLongPtrW.argtypes = [ctypes.c_void_p, ctypes.c_int]

        user32.SetWindowLongPtrW.restype  = ctypes.c_ssize_t
        user32.SetWindowLongPtrW.argtypes = [ctypes.c_void_p, ctypes.c_int,
                                             WNDPROC]

        user32.CallWindowProcW.restype  = ctypes.c_ssize_t
        user32.CallWindowProcW.argtypes = [ctypes.c_ssize_t,   # lpPrevWndFunc
                                           ctypes.c_void_p,    # hWnd
                                           ctypes.c_uint,      # Msg
                                           ctypes.c_size_t,    # wParam
                                           ctypes.c_ssize_t]   # lParam

        shell32.DragAcceptFiles.argtypes  = [ctypes.c_void_p, ctypes.c_bool]
        shell32.DragQueryFileW.argtypes   = [ctypes.c_void_p, ctypes.c_uint,
                                             ctypes.c_wchar_p, ctypes.c_uint]
        shell32.DragQueryFileW.restype    = ctypes.c_uint
        shell32.DragFinish.argtypes       = [ctypes.c_void_p]

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


def _scan_folder(folder: str) -> list[str]:
    """Recursively find all .graph files under *folder*, sorted by path."""
    matches: list[str] = []
    for root, _dirs, files in os.walk(folder):
        for f in sorted(files):
            if f.lower().endswith(".graph"):
                matches.append(os.path.join(root, f))
    return matches


# ---------------------------------------------------------------------------
# Layer 3: Renderer (PyVista / VTK – GPU-accelerated)
# ---------------------------------------------------------------------------

_WIN_W, _WIN_H = 1400, 900  # render window size in pixels


class _FileListPanel:
    """Right-side scrollable file list shown when a playlist is active."""

    PANEL_W  = 250   # panel width in pixels
    ITEM_H   = 20    # height per item row
    HEADER_H = 32    # header area height
    PAD      = 8     # left padding for text
    MAX_VIS  = 30    # maximum simultaneously visible items

    def __init__(self, renderer, vtk_mod, files: list[str],
                 idx: int, win_w: int, win_h: int) -> None:
        self._r     = renderer
        self._vtk   = vtk_mod
        self.files  = files
        self.idx    = max(0, min(idx, len(files) - 1))
        self._scroll = 0
        self._n_vis = min(len(files), self.MAX_VIS)
        self._actors:      list = []
        self._item_actors: list = []

        margin    = 10
        self._x1  = win_w - self.PANEL_W - margin
        self._x2  = win_w - margin
        self._y2  = win_h - 95                                          # top edge
        self._y1  = (self._y2
                     - self.HEADER_H
                     - self._n_vis * self.ITEM_H
                     - 22)                                               # bottom edge

        self._build()
        self._refresh()

    # ------------------------------------------------------------------ helpers

    def _rect(self, x1: float, y1: float, x2: float, y2: float,
              rgb: tuple, alpha: float = 1.0):
        vtk = self._vtk
        pts = vtk.vtkPoints()
        for x, y in ((x1, y1), (x2, y1), (x2, y2), (x1, y2)):
            pts.InsertNextPoint(x, y, 0.0)
        cells = vtk.vtkCellArray()
        cells.InsertNextCell(4)
        for i in range(4):
            cells.InsertCellPoint(i)
        pd = vtk.vtkPolyData()
        pd.SetPoints(pts)
        pd.SetPolys(cells)
        mapper = vtk.vtkPolyDataMapper2D()
        mapper.SetInputData(pd)
        # GetTransformCoordinate() returns None until one is assigned; create it first
        _coord = vtk.vtkCoordinate()
        _coord.SetCoordinateSystemToDisplay()
        mapper.SetTransformCoordinate(_coord)
        actor = vtk.vtkActor2D()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(*rgb)
        actor.GetProperty().SetOpacity(alpha)
        self._r.AddActor2D(actor)
        self._actors.append(actor)
        return actor

    def _text(self, label: str, x: float, y: float,
              size: int, rgb: tuple) -> object:
        a = self._vtk.vtkTextActor()
        a.SetInput(label)
        a.GetPositionCoordinate().SetCoordinateSystemToDisplay()
        a.SetPosition(x, y)
        tp = a.GetTextProperty()
        tp.SetFontSize(size)
        tp.SetColor(*rgb)
        self._r.AddActor2D(a)
        self._actors.append(a)
        return a

    # ------------------------------------------------------------------ build

    def _build(self) -> None:
        # Dark background panel
        self._rect(self._x1, self._y1, self._x2, self._y2,
                   (0.11, 0.13, 0.22), alpha=0.91)
        # Thin accent line below header
        sep_y = self._y2 - self.HEADER_H
        self._rect(self._x1, sep_y - 1, self._x2, sep_y,
                   (0.28, 0.44, 0.78))
        # Header text
        self._text(
            f"  Graph Files  ({len(self.files)})",
            self._x1, self._y2 - self.HEADER_H + 9,
            11, (0.80, 0.90, 1.00),
        )
        # Item slots (one text actor per visible row)
        for _ in range(self._n_vis):
            a = self._vtk.vtkTextActor()
            a.GetPositionCoordinate().SetCoordinateSystemToDisplay()
            a.GetTextProperty().SetFontSize(10)
            self._r.AddActor2D(a)
            self._item_actors.append(a)
            self._actors.append(a)
        # Navigation hint at bottom
        self._text(
            "  [^] [v] navigate",
            self._x1, self._y1 + 4,
            9, (0.38, 0.50, 0.65),
        )

    # ------------------------------------------------------------------ state

    def _item_y(self, slot: int) -> float:
        """Bottom-left y pixel of item row *slot* (0 = top row)."""
        return self._y2 - self.HEADER_H - (slot + 1) * self.ITEM_H + 3

    def _refresh(self) -> None:
        """Recompute scroll offset and repaint all item actors."""
        if self.idx < self._scroll:
            self._scroll = self.idx
        elif self.idx >= self._scroll + self._n_vis:
            self._scroll = max(0, self.idx - self._n_vis + 1)

        for slot, actor in enumerate(self._item_actors):
            fi = self._scroll + slot
            if fi >= len(self.files):
                actor.SetInput("")
                continue
            name = os.path.basename(self.files[fi])
            if len(name) > 29:
                name = name[:27] + ".."
            selected = (fi == self.idx)
            actor.SetInput(f"  {'> ' if selected else '  '}{name}")
            actor.SetPosition(self._x1, self._item_y(slot))
            tp = actor.GetTextProperty()
            if selected:
                tp.SetColor(1.0, 1.0, 1.0)
                tp.SetBackgroundColor(0.18, 0.38, 0.78)
                tp.SetBackgroundOpacity(0.88)
                tp.SetBold(True)
            else:
                tp.SetColor(0.68, 0.82, 0.96)
                tp.SetBackgroundOpacity(0.0)
                tp.SetBold(False)

    def select(self, idx: int) -> None:
        """Move selection to *idx* and repaint."""
        self.idx = max(0, min(idx, len(self.files) - 1))
        self._refresh()

    def remove(self) -> None:
        """Remove all panel actors from the renderer."""
        for a in self._actors:
            try:
                self._r.RemoveActor2D(a)
            except Exception:
                pass
        self._actors.clear()
        self._item_actors.clear()


def visualize(
    graph: SheetMetalGraph,
    title: str = "Sheet Metal Graph",
    show_labels: bool = False,
    playlist: Optional[list[str]] = None,
    playlist_idx: int = 0,
) -> tuple[Optional[str], Optional[list[str]]]:
    """Render *graph* in a 3D window.

    If *playlist* is provided the right-side panel is shown and Up/Down
    arrows navigate through files in-place (no window close/reopen).

    Returns (next_path, next_playlist) when the user drops / opens a new
    file or folder, or (None, None) when they just close the window.
    """
    import vtk

    _next_file:     list[Optional[str]]             = [None]
    _next_playlist: list[Optional[list[str]]]       = [None]
    _cur_graph:     list[SheetMetalGraph]            = [graph]

    pv.global_theme.background = "white"
    plotter = pv.Plotter(title="Sheet Metal Graph Visualizer",
                         window_size=[_WIN_W, _WIN_H])
    plotter.enable_trackball_style()

    # Remap zoom from Up/Down (needed for playlist nav) to +/- keys
    def _zoom_in():
        plotter.camera.zoom(1.1)
        plotter.render()
    def _zoom_out():
        plotter.camera.zoom(0.9)
        plotter.render()
    plotter.add_key_event("plus", _zoom_in)
    plotter.add_key_event("equal", _zoom_in)   # +/= are same key
    plotter.add_key_event("minus", _zoom_out)
    plotter.add_key_event("underscore", _zoom_out)

    # Actors dict — checkbox callbacks close over the *lists* so we always
    # clear() in-place on rebuild (never reassign the list references).
    actors: dict[str, list] = {
        "faces": [], "edges": [], "bend_edges": [],
        "vertices": [], "labels": [], "bend_annots": [],
        "x_axis": [],
    }
    _lbl_actors: list = []   # add_point_labels actors removed manually on rebuild

    # ------------------------------------------------------------------ #
    #  Updatable title actor                                              #
    # ------------------------------------------------------------------ #
    _title_actor = vtk.vtkTextActor()
    _title_actor.GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
    _title_actor.SetPosition(0.5, 0.975)
    _tp = _title_actor.GetTextProperty()
    _tp.SetFontSize(10)
    _tp.SetColor(0.0, 0.0, 0.0)
    _tp.SetJustificationToCentered()
    plotter.renderer.AddActor2D(_title_actor)

    def _set_title(g: SheetMetalGraph, path: str) -> None:
        nv  = len(g.vertices)
        ne  = len(g.edges)
        nb  = len(g.bend_angles)
        nf  = len(g.coarse_face_lfaces)
        nlf = len(g.lamina_face_edges)
        _title_actor.SetInput(
            f"{os.path.basename(path)}  |  "
            f"{nv}V  {ne}E  {nb} bends  {nf} faces ({nlf} laminas)"
        )

    # ------------------------------------------------------------------ #
    #  Hover text actor (created early so _build_graph can reset it)     #
    # ------------------------------------------------------------------ #
    _hover_actor = vtk.vtkTextActor()
    _hover_actor.SetInput("")
    _hover_actor.GetTextProperty().SetFontSize(13)
    _hover_actor.GetTextProperty().SetColor(0.05, 0.05, 0.05)
    _hover_actor.GetTextProperty().SetBackgroundColor(1.0, 1.0, 0.85)
    _hover_actor.GetTextProperty().SetBackgroundOpacity(0.85)
    _hover_actor.GetTextProperty().SetShadow(True)
    _hover_actor.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
    _hover_actor.SetPosition(0.62, 0.02)
    plotter.renderer.AddActor2D(_hover_actor)
    _hover_last: list[str] = [""]

    # ------------------------------------------------------------------ #
    #  Always-visible 3-D cursor coordinate display (bottom-left)        #
    # ------------------------------------------------------------------ #
    _coord_actor = vtk.vtkTextActor()
    _coord_actor.SetInput("X: ---  Y: ---  Z: --- mm")
    _coord_actor.GetTextProperty().SetFontSize(12)
    _coord_actor.GetTextProperty().SetColor(0.05, 0.05, 0.05)
    _coord_actor.GetTextProperty().SetBackgroundColor(0.96, 0.98, 1.0)
    _coord_actor.GetTextProperty().SetBackgroundOpacity(0.80)
    _coord_actor.GetPositionCoordinate().SetCoordinateSystemToDisplay()
    _coord_actor.SetPosition(10, 10)
    plotter.renderer.AddActor2D(_coord_actor)

    # ------------------------------------------------------------------ #
    #  Graph build / rebuild (called once at start and on navigation)    #
    # ------------------------------------------------------------------ #
    def _build_graph(g: SheetMetalGraph, path: str) -> None:
        # Remove named mesh actors
        for _name in ("faces", "edges", "bend_edges", "vertices", "x_axis"):
            try:
                plotter.remove_actor(_name)
            except Exception:
                pass
        # Remove point-label actors (no name= in add_point_labels)
        for _a in _lbl_actors:
            try:
                plotter.remove_actor(_a)
            except Exception:
                pass
        _lbl_actors.clear()
        for _key in actors:
            actors[_key].clear()   # in-place clear so checkbox closures still work

        _cur_graph[0] = g
        _set_title(g, path)
        _hover_last[0] = ""
        _hover_actor.SetInput("")

        bend_eids = g.bend_edge_ids
        bad_faces: list[int] = []

        # -- Faces --
        face_pts:   list = []
        face_cells: list[int] = []
        fid_per_pt: list[int] = []
        fid_per_cell: list[int] = []
        pt_off = 0
        for face_id in g.lamina_face_edges:
            poly, ok = _extract_face_polygon(g, face_id)
            if ok and len(poly) >= 3:
                n = len(poly)
                face_pts.extend(poly)
                fid_per_pt.extend([face_id] * n)
                face_cells.extend([n] + list(range(pt_off, pt_off + n)))
                fid_per_cell.append(face_id)
                pt_off += n
            elif not ok:
                bad_faces.append(face_id)

        if face_pts:
            fm = pv.PolyData(np.array(face_pts, dtype=float),
                             np.array(face_cells, dtype=int))
            fm.point_data["face_id"]  = np.array(fid_per_pt,  dtype=int)
            fm.cell_data["color_key"] = np.array([i % 8 for i in fid_per_cell], dtype=float)
            actors["faces"].append(plotter.add_mesh(
                fm, scalars="color_key", cmap="Pastel1", show_scalar_bar=False,
                opacity=0.75, show_edges=False, lighting=True,
                pickable=True, name="faces", clim=[0, 7],
            ))

        # -- Edges --
        reg_pts: list = [];  reg_lines: list[int] = []
        bnd_pts: list = [];  bnd_lines: list[int] = []
        for eid, (v1, v2) in g.edges.items():
            if v1 not in g.vertices or v2 not in g.vertices:
                continue
            p1 = list(g.vertices[v1]);  p2 = list(g.vertices[v2])
            if eid in bend_eids:
                i = len(bnd_pts);  bnd_pts.extend([p1, p2]);  bnd_lines.extend([2, i, i + 1])
            else:
                i = len(reg_pts);  reg_pts.extend([p1, p2]);  reg_lines.extend([2, i, i + 1])

        if reg_pts:
            m = pv.PolyData(np.array(reg_pts, dtype=float))
            m.lines = np.array(reg_lines, dtype=int)
            actors["edges"].append(
                plotter.add_mesh(m, color="#555555", line_width=1.5,
                                 pickable=False, name="edges"))
        if bnd_pts:
            m = pv.PolyData(np.array(bnd_pts, dtype=float))
            m.lines = np.array(bnd_lines, dtype=int)
            actors["bend_edges"].append(
                plotter.add_mesh(m, color="#E87722", line_width=3.0,
                                 pickable=False, name="bend_edges"))

        # -- Vertices --
        cluster_map = classify_vertices(g.vertices)
        vid_list    = sorted(g.vertices.keys())
        coords      = np.array([g.vertices[v] for v in vid_list], dtype=float)
        pt_cloud    = pv.PolyData(coords)
        pt_cloud.point_data["cluster"]   = np.array([cluster_map[v] for v in vid_list], dtype=float)
        pt_cloud.point_data["vertex_id"] = np.array(vid_list, dtype=int)
        actors["vertices"].append(plotter.add_mesh(
            pt_cloud, scalars="cluster", cmap="tab10", clim=[0, 9],
            point_size=8, render_points_as_spheres=True,
            show_scalar_bar=False, pickable=True, name="vertices",
        ))

        # -- Vertex labels (optional) --
        if show_labels and vid_list:
            a = plotter.add_point_labels(
                pt_cloud, [str(v) for v in vid_list],
                font_size=8, text_color="black",
                show_points=False, always_visible=False, name="labels",
            )
            actors["labels"].append(a);  _lbl_actors.append(a)

        # -- Bend annotations --
        ann_pts:   list = [];  ann_txts: list[str] = []
        for bid, angle_rad in g.bend_angles.items():
            c = g.bend_centroid(bid)
            if c:
                ann_pts.append(list(c))
                ann_txts.append(f"B{bid}  {math.degrees(angle_rad):.0f}°")
        if ann_pts:
            ac = pv.PolyData(np.array(ann_pts, dtype=float))
            a  = plotter.add_point_labels(
                ac, ann_txts, font_size=10, text_color="#C05000", bold=True,
                show_points=False, always_visible=True, name="bend_annots",
            )
            actors["bend_annots"].append(a);  _lbl_actors.append(a)

        if bad_faces:
            print(f"WARNING: {len(bad_faces)} malformed face(s) skipped: {bad_faces}")

        # -- X-axis reference line (3 m = 3000 mm, centered at origin) --
        x_line = pv.Line((-1500.0, 0.0, 0.0), (1500.0, 0.0, 0.0))
        actors["x_axis"].append(plotter.add_mesh(
            x_line, color="#CC2222", line_width=2.5,
            pickable=False, name="x_axis",
        ))

        plotter.reset_camera()

    # ---- Initial build ----
    _build_graph(graph, title)

    # ------------------------------------------------------------------ #
    #  Layer toggle checkboxes                                            #
    # ------------------------------------------------------------------ #
    TOGGLES = [
        ("faces",       "Faces",       True,        "#D4E8F0"),
        ("edges",       "Edges",       True,        "#555555"),
        ("bend_edges",  "Bend Edges",  True,        "#E87722"),
        ("vertices",    "Vertices",    True,        "#3366AA"),
        ("labels",      "Labels",      show_labels, "#333333"),
        ("bend_annots", "Bend Annots", True,        "#C05000"),
        ("x_axis",      "X-Axis Line", True,        "#CC2222"),
    ]
    CB_SIZE, CB_GAP = 25, 5
    _start_y = _WIN_H - 100

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

    for _idx, (layer_key, layer_label, initial, color) in enumerate(TOGGLES):
        y_px = _start_y - _idx * (CB_SIZE + CB_GAP)
        plotter.add_checkbox_button_widget(
            _make_toggle(layer_key), value=initial,
            position=(10, y_px), size=CB_SIZE, border_size=2,
            color_on=color, color_off="#CCCCCC",
        )
        plotter.add_text(layer_label, position=(45, y_px + 4),
                         font_size=9, color="black")

    # ------------------------------------------------------------------ #
    #  Hover inspection                                                   #
    # ------------------------------------------------------------------ #
    _hover_picker    = vtk.vtkPointPicker()
    _hover_picker.SetTolerance(0.01)
    _hover_renderer  = plotter.renderer
    _hover_next_t:   list[float] = [0.0]
    _HOVER_INTERVAL  = 1.0 / 60.0

    def _on_mouse_move(obj, event):
        import time
        now = time.perf_counter()
        if now < _hover_next_t[0]:
            return
        _hover_next_t[0] = now + _HOVER_INTERVAL
        try:
            x, y = plotter.iren.get_event_position()
        except AttributeError:
            x, y = plotter.iren.GetEventPosition()
        _hover_picker.Pick(x, y, 0, _hover_renderer)
        pid  = _hover_picker.GetPointId()
        dset = _hover_picker.GetDataSet()

        # --- 3-D world-space cursor coordinates (always shown) ---
        if dset is not None and pid >= 0:
            wx, wy, wz = _hover_picker.GetPickPosition()
        else:
            # Project mouse onto the camera focal plane for a meaningful 3-D pos
            cam = _hover_renderer.GetActiveCamera()
            fp  = cam.GetFocalPoint()
            _hover_renderer.SetWorldPoint(fp[0], fp[1], fp[2], 1.0)
            _hover_renderer.WorldToDisplay()
            dz = _hover_renderer.GetDisplayPoint()[2]
            _hover_renderer.SetDisplayPoint(x, y, dz)
            _hover_renderer.DisplayToWorld()
            wp = _hover_renderer.GetWorldPoint()
            if wp[3] != 0.0:
                wx, wy, wz = wp[0] / wp[3], wp[1] / wp[3], wp[2] / wp[3]
            else:
                wx, wy, wz = 0.0, 0.0, 0.0
        _coord_actor.SetInput(f"X: {wx:.2f}  Y: {wy:.2f}  Z: {wz:.2f} mm")

        # --- Hover inspection ---
        new_text = ""
        if dset is not None and pid >= 0:
            pd      = dset.GetPointData()
            vid_arr = pd.GetArray("vertex_id")
            fid_arr = pd.GetArray("face_id")
            g = _cur_graph[0]
            if vid_arr is not None and pid < vid_arr.GetNumberOfTuples():
                vid = int(vid_arr.GetValue(pid))
                x_, y_, z_ = g.vertices.get(vid, (0.0, 0.0, 0.0))
                new_text = (f"Vertex {vid}\n  x = {x_:.3f} mm\n"
                            f"  y = {y_:.3f} mm\n  z = {z_:.3f} mm")
            elif fid_arr is not None and pid < fid_arr.GetNumberOfTuples():
                fid    = int(fid_arr.GetValue(pid))
                nedges = len(g.lamina_face_edges.get(fid, []))
                new_text = f"Lamina Face {fid}\n  {nedges} boundary edges"
        if new_text != _hover_last[0]:
            _hover_last[0] = new_text
            _hover_actor.SetInput(new_text)

    try:
        plotter.iren.AddObserver("MouseMoveEvent", _on_mouse_move)
    except AttributeError:
        _raw = (getattr(plotter.iren, "_iren", None)
                or getattr(plotter.iren, "interactor", None))
        if _raw is not None:
            _raw.AddObserver("MouseMoveEvent", _on_mouse_move)

    # ------------------------------------------------------------------ #
    #  Drag-and-drop (Win32 WndProc; handles both files and folders)     #
    # ------------------------------------------------------------------ #
    def _get_iren():
        iren = plotter.iren
        return getattr(iren, "interactor", None) or getattr(iren, "_iren", None) or iren

    _dnd_ready = [False]

    def _on_render_start(obj, event):
        if _dnd_ready[0]:
            return

        def _on_dropped(raw_path: str) -> None:
            if os.path.isdir(raw_path):
                files = _scan_folder(raw_path)
                if not files:
                    print(f"No .graph files found in: {raw_path}")
                    return
                _next_file[0]     = files[0]
                _next_playlist[0] = files
            else:
                _next_file[0]     = raw_path
                _next_playlist[0] = None
            _get_iren().TerminateApp()

        if _setup_win32_dnd(plotter.render_window, _on_dropped):
            _dnd_ready[0] = True

    plotter.render_window.AddObserver("StartEvent", _on_render_start)

    # ------------------------------------------------------------------ #
    #  File list panel + Up/Down navigation (playlist mode only)         #
    # ------------------------------------------------------------------ #
    if playlist and len(playlist) > 1:
        _pl_idx = [playlist_idx]
        _panel  = _FileListPanel(plotter.renderer, vtk, playlist,
                                 playlist_idx, _WIN_W, _WIN_H)

        def _navigate(delta: int) -> None:
            new_idx = _pl_idx[0] + delta
            if not (0 <= new_idx < len(playlist)):
                return
            new_path = playlist[new_idx]
            try:
                new_graph = _load_graph(new_path)
            except SystemExit as exc:
                print(exc)
                return
            _pl_idx[0] = new_idx
            _panel.select(new_idx)
            _build_graph(new_graph, new_path)
            plotter.render()

        plotter.add_key_event("Up",   lambda: _navigate(-1))
        plotter.add_key_event("Down", lambda: _navigate( 1))

    # ------------------------------------------------------------------ #
    #  View preset keyboard shortcuts                                     #
    # ------------------------------------------------------------------ #
    _VIEW_KEYS: dict[str, tuple] = {
        "1": ((0.0,  0.0,  1.0), (0.0, 1.0, 0.0)),   # Top
        "2": ((0.0, -1.0,  0.0), (0.0, 0.0, 1.0)),   # Front
        "3": ((1.0,  0.0,  0.0), (0.0, 0.0, 1.0)),   # Side
        "4": ((1.0, -1.0,  1.0), (0.0, 0.0, 1.0)),   # Iso
        "5": ((0.0,  1.0,  0.0), (0.0, 0.0, 1.0)),   # Back
        "0": ((1.0, -1.0,  1.0), (0.0, 0.0, 1.0)),   # Home
    }

    def _make_view_cb(pos_dir: tuple, up: tuple):
        def _cb():
            plotter.view_vector(pos_dir, viewup=up)
            plotter.reset_camera()
        return _cb

    for _key, (_pos_dir, _up) in _VIEW_KEYS.items():
        plotter.add_key_event(_key, _make_view_cb(_pos_dir, _up))

    def _reset_view():
        plotter.view_vector((1.0, -1.0, 1.0), viewup=(0.0, 0.0, 1.0))
        plotter.reset_camera()

    plotter.add_key_event("r", _reset_view)

    def _open_file():
        path = _pick_file_dialog()
        if path:
            _next_file[0]     = path
            _next_playlist[0] = None
            _get_iren().TerminateApp()

    def _open_folder():
        folder = _pick_folder_dialog()
        if not folder:
            return
        files = _scan_folder(folder)
        if not files:
            print(f"No .graph files found in: {folder}")
            return
        _next_file[0]     = files[0]
        _next_playlist[0] = files
        _get_iren().TerminateApp()

    plotter.add_key_event("o", _open_file)
    plotter.add_key_event("O", _open_file)
    plotter.add_key_event("f", _open_folder)
    plotter.add_key_event("F", _open_folder)

    # ------------------------------------------------------------------ #
    #  Screenshot                                                         #
    # ------------------------------------------------------------------ #
    import datetime

    def _save_screenshot():
        ts    = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        fname = f"sheet_metal_{ts}.png"
        plotter.screenshot(fname, transparent_background=False)
        print(f"Screenshot saved: {fname}")

    plotter.add_key_event("s", _save_screenshot)

    # ------------------------------------------------------------------ #
    #  Axes + console summary                                             #
    # ------------------------------------------------------------------ #
    plotter.add_axes(line_width=2)

    print("\n--- Controls ---")
    print("  Keys : 1=Top  2=Front  3=Side  4=Iso  5=Back  0=Home  R=Reset  S=Screenshot")
    print("         O=Open file  F=Open folder  +/-=Zoom")
    if playlist and len(playlist) > 1:
        print(f"         Up/Down = navigate playlist ({len(playlist)} files)")
    print("  Mouse: Left-drag=Orbit  Middle-drag=Pan  Scroll=Zoom")
    print("  Hover over geometry to inspect elements (vertex coords, face info)")
    print("  Top-left checkboxes: toggle layer visibility")
    print("  Drag a file or folder onto the window to reload / open playlist")

    plotter.show(auto_close=False)
    try:
        plotter.close()
    except Exception:
        pass
    return _next_file[0], _next_playlist[0]


# ---------------------------------------------------------------------------
# Entry point helpers
# ---------------------------------------------------------------------------

_REQUIRED_GRAPH_KEYS = [
    "VertexNodes", "EdgeToVertexArchs", "EdgeNodes",
    "LaminaFaceToEdgeArchs", "BendToLaminaFaceArchs",
    "BendNodes", "FaceToLaminaFaceArchs",
]


def _load_graph(path: str) -> SheetMetalGraph:
    # Try encodings in order — Windows often saves files as CP-1252, not UTF-8
    raw = None
    for enc in ("utf-8-sig", "utf-8", "cp1252", "latin-1"):
        try:
            with open(path, encoding=enc) as f:
                raw = json.load(f)
            break
        except FileNotFoundError:
            raise SystemExit(f"ERROR: File not found:\n{path}")
        except PermissionError:
            raise SystemExit(f"ERROR: Permission denied reading:\n{path}")
        except UnicodeDecodeError:
            continue
        except json.JSONDecodeError as exc:
            raise SystemExit(f"ERROR: Invalid JSON in:\n{path}\n\n{exc}")
        except OSError as exc:
            raise SystemExit(f"ERROR: Cannot open '{path}':\n{exc}")

    if raw is None:
        raise SystemExit(
            f"ERROR: Could not decode '{os.path.basename(path)}'.\n"
            "Tried utf-8, cp1252, latin-1 — file may be corrupt."
        )

    if "Graph" not in raw:
        raise SystemExit(
            f"ERROR: JSON root must contain a 'Graph' key.\n"
            f"Keys found: {list(raw.keys())}"
        )
    graph_dict = raw["Graph"]
    missing = [k for k in _REQUIRED_GRAPH_KEYS if k not in graph_dict]
    if missing:
        raise SystemExit(
            f"ERROR: 'Graph' section is missing required keys:\n{missing}"
        )
    try:
        return SheetMetalGraph.from_dict(graph_dict)
    except Exception as exc:
        raise SystemExit(
            f"ERROR: Failed to parse graph data in:\n{path}\n\n{exc}"
        ) from exc


def _show_error(message: str) -> None:
    """Show a blocking error dialog — works even in .pyw (no console) mode."""
    try:
        import tkinter as tk
        from tkinter import messagebox
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        messagebox.showerror("Graph Load Error", message)
        root.destroy()
    except Exception:
        pass


def _pick_file_dialog() -> Optional[str]:
    """Native OS file-open dialog.  Returns selected path or None."""
    try:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        root.lift()
        p = filedialog.askopenfilename(
            title="Open Graph File",
            filetypes=[("Graph files", "*.graph *.txt *.json"),
                       ("All files", "*.*")],
        )
        root.destroy()
        return p.strip() or None
    except Exception:
        return None


def _pick_folder_dialog() -> Optional[str]:
    """Native OS folder-browse dialog.  Returns selected folder or None."""
    _dbg("_pick_folder_dialog: opening dialog")
    try:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        root.lift()
        p = filedialog.askdirectory(title="Select Folder with .graph Files")
        root.destroy()
        result = p.strip() or None
        _dbg(f"_pick_folder_dialog: returning {result!r}")
        return result
    except Exception as _e:
        _dbg(f"_pick_folder_dialog: EXCEPTION {_e}")
        return None


def _welcome_screen() -> tuple[Optional[str], Optional[list[str]]]:
    """Splash window shown when no file is provided at startup.
    Returns (first_path, playlist) — playlist is None for single files."""
    import vtk as _vtk

    pv.global_theme.background = "#1C1C2E"
    pl = pv.Plotter(
        title="Sheet Metal Graph Visualizer",
        window_size=[640, 400],
    )
    pl.enable_trackball_style()

    # Title — centred via vtkTextActor so justification is exact
    _ta = _vtk.vtkTextActor()
    _ta.SetInput("Sheet Metal Graph Visualizer")
    _ta.GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
    _ta.SetPosition(0.5, 0.82)
    _tp2 = _ta.GetTextProperty()
    _tp2.SetFontSize(18)
    _tp2.SetColor(1.0, 1.0, 1.0)
    _tp2.SetJustificationToCentered()
    _tp2.SetBold(True)
    pl.renderer.AddActor2D(_ta)

    # Instructions — centred, placed in the lower half of the window
    _ia = _vtk.vtkTextActor()
    _ia.SetInput(
        "Drop a .graph file or folder onto this window\n\n"
        "O  =  open file        F  =  open folder\n\n"
        "Q  =  quit"
    )
    _ia.GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
    _ia.SetPosition(0.5, 0.28)
    _ip = _ia.GetTextProperty()
    _ip.SetFontSize(12)
    _ip.SetColor(0.67, 0.67, 0.67)
    _ip.SetJustificationToCentered()
    pl.renderer.AddActor2D(_ia)

    result_path:     list[Optional[str]]        = [None]
    result_playlist: list[Optional[list[str]]]  = [None]

    def _get_iren():
        iren = pl.iren
        return getattr(iren, "interactor", None) or getattr(iren, "_iren", None) or iren

    _dnd_ready = [False]

    def _on_render_start(obj, event):
        if _dnd_ready[0]:
            return

        def _on_dropped(raw_path: str) -> None:
            _dbg(f"_on_dropped (welcome): {raw_path!r}  isdir={os.path.isdir(raw_path)}")
            if os.path.isdir(raw_path):
                files = _scan_folder(raw_path)
                _dbg(f"_on_dropped (welcome): scan found {len(files)} files")
                if not files:
                    return
                result_path[0]     = files[0]
                result_playlist[0] = files
            else:
                result_path[0] = raw_path
            _dbg("_on_dropped (welcome): calling TerminateApp")
            _get_iren().TerminateApp()

        ok = _setup_win32_dnd(pl.render_window, _on_dropped)
        _dbg(f"_on_render_start (welcome): _setup_win32_dnd returned {ok}")
        if ok:
            _dnd_ready[0] = True

    pl.render_window.AddObserver("StartEvent", _on_render_start)

    def _open_file():
        path = _pick_file_dialog()
        if path:
            result_path[0] = path
            _get_iren().TerminateApp()

    def _open_folder():
        _dbg("_open_folder (welcome): called")
        folder = _pick_folder_dialog()
        if not folder:
            _dbg("_open_folder (welcome): no folder selected")
            return
        files = _scan_folder(folder)
        _dbg(f"_open_folder (welcome): folder={folder!r}  files={len(files)}")
        if not files:
            return
        result_path[0]     = files[0]
        result_playlist[0] = files
        _dbg(f"_open_folder (welcome): calling TerminateApp, first={files[0]!r}")
        _get_iren().TerminateApp()

    pl.add_key_event("o", _open_file)
    pl.add_key_event("O", _open_file)
    pl.add_key_event("f", _open_folder)
    pl.add_key_event("F", _open_folder)

    pl.show(auto_close=False)
    try:
        pl.close()
    except Exception:
        pass

    pv.global_theme.background = "white"
    return result_path[0], result_playlist[0]


def main() -> None:
    import sys
    import datetime

    # Clear the debug log at start of each run for a clean read
    _dbg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "debug.log")
    try:
        with open(_dbg_path, "w", encoding="utf-8") as _f:
            _f.write(f"=== run {datetime.datetime.now()} ===\n")
    except Exception:
        pass

    _dbg(f"main: argv={sys.argv}")

    show_labels = "--labels" in sys.argv
    cli_args    = [a for a in sys.argv[1:] if not a.startswith("-")]

    if cli_args:
        raw = cli_args[0]
        if os.path.isdir(raw):
            playlist: Optional[list[str]] = _scan_folder(raw)
            path: Optional[str]           = playlist[0] if playlist else None
        else:
            playlist = None
            path     = raw
    else:
        path, playlist = _welcome_screen()
        _dbg(f"main: welcome returned path={path!r}  playlist_len={len(playlist) if playlist else 0}")

    while path:
        try:
            graph = _load_graph(path)
        except (SystemExit, Exception) as exc:
            msg = str(exc).replace("SystemExit: ", "")
            _dbg(f"load error: {msg}")
            print(msg)
            _show_error(msg)
            path, playlist = _welcome_screen()
            continue

        pl_idx = 0
        if playlist:
            try:
                pl_idx = playlist.index(path)
            except ValueError:
                pass

        _dbg(f"main: calling visualize  path={path!r}  playlist_len={len(playlist) if playlist else 0}")
        path, playlist = visualize(
            graph, title=path, show_labels=show_labels,
            playlist=playlist, playlist_idx=pl_idx,
        )
        _dbg(f"main: visualize returned path={path!r}")


if __name__ == "__main__":
    main()
