# Building a standalone ampliscan.exe (Windows)

This produces a single double-click executable that opens the GUI, so users who
don't have Python can run ampliscan. Do this **after** placing your icon (see
`src/ampliscan/assets/README.md`).

## 1. Install build tools

```bash
pip install -e .[fast,excel]
pip install pyinstaller pillow
```

## 2. Generate the icon (once)

```bash
python scripts/make_icon.py path/to/your_logo.png
```

## 3. Build

```bash
pyinstaller --noconfirm --windowed --onefile ^
  --name ampliscan ^
  --icon src/ampliscan/assets/ampliscan_icon.ico ^
  --add-data "src/ampliscan/assets;ampliscan/assets" ^
  --collect-submodules ampliscan ^
  gui_launch.py
```

Where `gui_launch.py` is a one-line entry script:

```python
from ampliscan.gui import launch
launch()
```

Notes:

- `--windowed` hides the console so it looks like a native app.
- `--add-data "SRC;DEST"` uses `;` on Windows and `:` on macOS/Linux.
- The `.exe` appears in `dist/ampliscan.exe`. Ship that single file.
- On macOS, use `--icon ...icns` (convert the PNG to `.icns`) and the same
  `--add-data` with a `:` separator; the result is `dist/ampliscan.app`.

## 4. Result

`dist/ampliscan.exe` — a self-contained GUI application with your DNA icon in
the title bar, taskbar, and file icon. No Python required on the target machine.
