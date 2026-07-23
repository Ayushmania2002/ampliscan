"""Convert the ampliscan logo PNG into a multi-size Windows .ico.

Usage:
    python scripts/make_icon.py path/to/your_icon.png

Writes src/ampliscan/assets/ampliscan_icon.ico (used for the GUI window,
the Windows taskbar, and the PyInstaller executable) and also copies a
square PNG to src/ampliscan/assets/ampliscan_icon.png (cross-platform).

Requires Pillow:  pip install pillow
"""
import sys
from pathlib import Path


def main(src_png: str) -> int:
    try:
        from PIL import Image
    except ImportError:
        print("Pillow is required:  pip install pillow")
        return 1

    src = Path(src_png)
    if not src.exists():
        print(f"No such file: {src}")
        return 1

    assets = Path(__file__).resolve().parent.parent / "src" / "ampliscan" / "assets"
    assets.mkdir(parents=True, exist_ok=True)

    img = Image.open(src).convert("RGBA")

    # pad to a square canvas so the icon isn't distorted
    w, h = img.size
    side = max(w, h)
    square = Image.new("RGBA", (side, side), (0, 0, 0, 0))
    square.paste(img, ((side - w) // 2, (side - h) // 2))

    png_out = assets / "ampliscan_icon.png"
    square.resize((256, 256)).save(png_out)

    ico_out = assets / "ampliscan_icon.ico"
    sizes = [(16, 16), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)]
    square.save(ico_out, format="ICO", sizes=sizes)

    print(f"Wrote {png_out}")
    print(f"Wrote {ico_out}")
    print("The GUI will now show this icon; rebuild the exe to embed it.")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(2)
    raise SystemExit(main(sys.argv[1]))
