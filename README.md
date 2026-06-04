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
- **Drag-and-drop** — drop a `.graph` file or folder onto the window to load it
- **Folder scanning** — select a folder to scan for all `.graph` files and browse
  them with Up/Down arrow keys (playlist mode)
- **File/folder dialogs** — press `O` to open a file or `F` to open a folder
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

### GUI (double-click app launcher)

On Windows, double-click **`launch.pyw`** to open the visualizer with a welcome
screen. You can then:
- Drag a `.graph` file or folder onto the window
- Press **O** to open a file dialog
- Press **F** to open a folder dialog (scans for all `.graph` files)

When a folder is open, the right panel shows a playlist. Navigate with **Up/Down**.

### Command-line

```bash
# Basic
python visualize_graph.py path/to/grafo.txt

# With vertex ID labels
python visualize_graph.py path/to/grafo.txt --labels

# Open folder (scans for .graph files)
python visualize_graph.py path/to/folder

# Default (shows welcome screen)
python visualize_graph.py
```

---

## Controls

### Mouse & Scroll
| Action | Result |
|--------|--------|
| **Left-drag** | Orbit |
| **Middle-drag** | Pan |
| **Scroll** | Zoom |
| **Hover** | Inspect vertex / face under cursor |

### Keyboard
| Key | Result |
|-----|--------|
| **1–5, 0** | View presets (Top, Front, Side, Iso, Back, Home) |
| **R** | Reset camera to fit |
| **S** | Save screenshot as PNG |
| **O** | Open file dialog |
| **F** | Open folder dialog |
| **Up / Down** | Navigate playlist (when folder is open) |
| **+/−** | Zoom in/out |

### UI
| Element | Result |
|---------|--------|
| **Top-left checkboxes** | Toggle layer visibility |
| **Right panel** (in playlist mode) | File list with arrow navigation |

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
