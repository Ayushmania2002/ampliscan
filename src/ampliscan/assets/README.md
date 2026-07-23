# ampliscan assets

Drop your generated brand images here.

## Required for the app icon

Save your **square logo** (the DNA-helix-into-bars mark) and run the converter,
which produces both files the app looks for:

```bash
python scripts/make_icon.py path/to/your_logo.png
```

This writes:

- `ampliscan_icon.ico` — Windows title bar, taskbar, and the PyInstaller `.exe`
- `ampliscan_icon.png` — cross-platform window icon (macOS/Linux)

The GUI loads these automatically at startup (`_set_window_icon` in `gui.py`)
and silently skips them if they are absent, so the app always launches.

## Optional

- `cover.png` — the wide banner (used in the GitHub README / social preview).

These files are bundled into the wheel, so they ship with `pip install ampliscan`.
