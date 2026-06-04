"""
Double-click launcher for the Sheet Metal Graph Visualizer.

Using .pyw instead of .py runs the script with pythonw.exe on Windows,
which suppresses the console window so the app feels like a native GUI app.
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from visualize_graph import main

main()
