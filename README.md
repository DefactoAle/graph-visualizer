# Sheet Metal Graph Visualizer

Interactive 3D viewer for Salvagnini SMSerializer sheet metal graphs.
Reads the JSON topology file produced by SMSerializer and renders the full
sheet metal structure — panels, bend edges, vertices, and bend annotations —
in a GPU-accelerated VTK window.

![Sheet metal part rendered in isometric view](docs/screenshot.png)

---

## Features

- **GPU-accelerated rendering** via PyVista / VTK (OpenGL) — smooth at any
  graph size
- **Layer toggles** — show/hide Faces, Edges, Bend Edges, Vertices, Labels,
  and Bend Annotations independently with the checkboxes in the top-left corner
- **Hover inspection** — move the mouse over any vertex or face panel to see
  its ID, coordinates (mm), or edge count in a tooltip at the bottom-right
- **Screenshot** — press `S` to save a timestamped PNG to the working directory
- **View presets** — one key to jump to Top / Front / Side / Iso / Back
- **CAD-style navigation** — left-drag orbits, middle-drag pans, scroll zooms
- **Validated input** — helpful error messages for missing files, bad JSON, or
  incomplete graph data

---

## Requirements

| Package | Version | License |
|---------|---------|---------|
| Python  | 3.9+    | PSF     |
| NumPy   | any     | BSD-3   |
| PyVista | 0.44+   | MIT     |
| VTK     | 9.0+    | BSD-3   |

All dependencies are permissively licensed (BSD / MIT) and are free for
commercial use.

### Install

```bash
pip install "pyvista[all]" numpy
```

---

## Usage

```bash
# Basic
python visualize_graph.py path/to/grafo.txt

# With vertex ID labels
python visualize_graph.py path/to/grafo.txt --labels

# Default (looks for grafo.txt in the current directory)
python visualize_graph.py
```

---

## Controls

| Action | Result |
|--------|--------|
| **Left-drag** | Orbit |
| **Middle-drag** | Pan |
| **Scroll** | Zoom |
| **1** | Top view |
| **2** | Front view |
| **3** | Side view |
| **4** | Isometric view (default) |
| **5** | Back view |
| **0** | Home (reset to iso) |
| **R** | Reset camera to fit |
| **S** | Save screenshot as PNG |
| **Hover** | Inspect vertex / face under cursor |
| **Top-left checkboxes** | Toggle layer visibility |

---

## Input format

The visualizer reads the JSON file emitted by Salvagnini SMSerializer.
The root object must have a `"Graph"` key containing:

```
Graph
├── VertexNodes              [[id, x, y, z], ...]
├── EdgeNodes                [[type, weight], ...]
├── EdgeToVertexArchs        [[edge_ids], [vertex_ids]]
├── LaminaFaceToEdgeArchs    [[face_ids], [edge_ids]]
├── BendNodes                [[id, ..., angle_rad, ...], ...]
├── BendToLaminaFaceArchs    [[bend_seq], [face_ids]]
└── FaceToLaminaFaceArchs    [[coarse_face_ids], [lamina_face_ids]]
```

All coordinates are in **millimetres**. Bend angles are in **radians**.

---

## Architecture

```
visualize_graph.py
├── Layer 1 – Data model & parser
│   └── SheetMetalGraph  dataclass + from_dict() factory
├── Layer 2 – Vertex clustering
│   └── classify_vertices()  (z-percentile bucketing)
├── Layer 2.5 – Topology
│   └── _extract_face_polygon()  (edge-tracing)
├── Layer 3 – Renderer
│   └── visualize()  (PyVista / VTK)
└── Entry point
    └── _load_graph() + argparse CLI
```
