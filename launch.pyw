"""
Double-click launcher for the Sheet Metal Graph Visualizer.

Using .pyw instead of .py runs the script with pythonw.exe on Windows,
which suppresses the console window so the app feels like a native GUI app.

Any unhandled exception is written to  error.log  in the same folder so
crashes are diagnosable even without a console.
"""
import sys
import os
import traceback
import datetime

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _here)

_log_path = os.path.join(_here, "error.log")


def _write_log(msg: str) -> None:
    try:
        with open(_log_path, "a", encoding="utf-8") as _f:
            _f.write(msg)
    except Exception:
        pass


try:
    from visualize_graph import main
    main()
except Exception:
    _write_log(
        f"\n=== CRASH  {datetime.datetime.now()} ===\n"
        + traceback.format_exc()
    )
