"""Tkinter desktop GUI for ampliscan.

Uses only the Python standard library (tkinter), so it ships with the package
and needs no extra dependencies -- ``pip install ampliscan`` is enough. Launch
with ``ampliscan gui`` or ``ampliscan-gui``.

The heavy demux runs in a background thread so the window stays responsive; the
worker talks back to the UI through a thread-safe queue polled on the Tk event
loop.
"""
from __future__ import annotations

import queue
import threading
import traceback
from pathlib import Path
from typing import Optional

from ampliscan.config import Panel, MatchPolicy


def _require_tk():
    try:
        import tkinter as tk  # noqa: F401
        from tkinter import ttk  # noqa: F401
        return True
    except Exception:  # pragma: no cover - only on headless/no-Tk builds
        return False


# Small viridis-ish ramp for the heatmap cells.
_VIRIDIS = [
    (68, 1, 84), (59, 82, 139), (33, 145, 140), (94, 201, 98), (253, 231, 37),
]


def _viridis(t: float) -> str:
    t = max(0.0, min(1.0, t))
    s = t * (len(_VIRIDIS) - 1)
    i = int(s)
    f = s - i
    a = _VIRIDIS[i]
    b = _VIRIDIS[min(i + 1, len(_VIRIDIS) - 1)]
    c = tuple(round(a[k] + (b[k] - a[k]) * f) for k in range(3))
    return f"#{c[0]:02x}{c[1]:02x}{c[2]:02x}"


def _fmt_count(v: int) -> str:
    """Compact read count that fits a small heatmap cell (891396 -> '891k')."""
    if v >= 1_000_000:
        return f"{v / 1e6:.1f}M"
    if v >= 10_000:
        return f"{v / 1e3:.0f}k"
    return str(v)


def _parse_barcodes(text: str) -> dict:
    out = {}
    for line in text.strip().splitlines():
        parts = line.replace(",", " ").replace("\t", " ").split()
        if len(parts) >= 2 and parts[0] and parts[1]:
            out[parts[0]] = parts[1].upper()
    return out


def _asset(name: str):
    """Absolute path to a bundled asset (works installed or from source)."""
    return Path(__file__).resolve().parent / "assets" / name


def _set_window_icon(root) -> None:
    """Apply the ampliscan icon to the window/taskbar if the asset is present.

    Looks for src/ampliscan/assets/ampliscan_icon.ico (Windows title bar +
    taskbar + PyInstaller exe) and .png (cross-platform). Silently no-ops if
    the files aren't there yet, so the app always launches.
    """
    ico = _asset("ampliscan_icon.ico")
    png = _asset("ampliscan_icon.png")
    try:
        if ico.exists():
            root.iconbitmap(default=str(ico))
    except Exception:  # noqa: BLE001 - non-fatal; not all platforms take .ico
        pass
    try:
        if png.exists():
            import tkinter as tk
            img = tk.PhotoImage(file=str(png))
            root.iconphoto(True, img)
            root._ampliscan_icon_ref = img  # keep a reference so Tk doesn't GC it
    except Exception:  # noqa: BLE001
        pass


class AmpliscanGUI:
    # Modern flat palette (single white surface + one accent).
    PALETTE = {
        "bg": "#ffffff",
        "subtle": "#f3f4f6",
        "accent": "#2563eb",
        "accent_hover": "#1d4ed8",
        "text": "#111827",
        "muted": "#6b7280",
        "border": "#e5e7eb",
        "success": "#059669",
    }

    def __init__(self, root):
        import tkinter as tk
        from tkinter import ttk

        self.tk = tk
        self.ttk = ttk
        self.root = root
        root.title("ampliscan — amplicon demultiplexer")
        root.geometry("1200x760")
        root.minsize(980, 640)

        self.q: "queue.Queue" = queue.Queue()
        self.worker: Optional[threading.Thread] = None

        self._setup_style()
        self._build_widgets()
        self._poll_queue()

    # --- theming ------------------------------------------------------------

    def _setup_style(self):
        style = self.ttk.Style()
        try:
            style.theme_use("clam")  # clam honours custom colours; native themes don't
        except Exception:  # pragma: no cover
            pass
        P = self.PALETTE
        base = ("Segoe UI", 10)

        style.configure(".", background=P["bg"], foreground=P["text"],
                        font=base, borderwidth=0, focuscolor=P["bg"])
        style.configure("TFrame", background=P["bg"])
        style.configure("TLabel", background=P["bg"], foreground=P["text"])
        style.configure("Muted.TLabel", background=P["bg"], foreground=P["muted"],
                        font=("Segoe UI", 9))
        style.configure("Section.TLabel", background=P["bg"], foreground=P["muted"],
                        font=("Segoe UI", 9, "bold"))
        style.configure("Field.TLabel", background=P["bg"], foreground=P["muted"],
                        font=("Segoe UI", 9))
        style.configure("Title.TLabel", background=P["bg"], foreground=P["text"],
                        font=("Segoe UI", 21, "bold"))
        style.configure("Status.TLabel", background=P["bg"], foreground=P["accent"],
                        font=("Segoe UI", 10))

        # buttons
        style.configure("TButton", background=P["subtle"], foreground=P["text"],
                        relief="flat", padding=(10, 6), borderwidth=0)
        style.map("TButton", background=[("active", P["border"]), ("pressed", P["border"])])
        style.configure("Small.TButton", padding=(7, 4), font=("Segoe UI", 9))
        style.configure("Accent.TButton", background=P["accent"], foreground="#ffffff",
                        font=("Segoe UI", 11, "bold"), padding=(12, 11),
                        relief="flat", borderwidth=0)
        style.map("Accent.TButton",
                  background=[("active", P["accent_hover"]), ("disabled", "#9ca3af")],
                  foreground=[("disabled", "#f3f4f6")])

        # entries / spinboxes
        for s in ("TEntry", "TSpinbox"):
            style.configure(s, fieldbackground="#ffffff", background="#ffffff",
                            foreground=P["text"], bordercolor=P["border"],
                            lightcolor=P["border"], darkcolor=P["border"],
                            insertcolor=P["text"], relief="solid", borderwidth=1,
                            padding=5, arrowsize=13)
            style.map(s, bordercolor=[("focus", P["accent"])],
                      lightcolor=[("focus", P["accent"])],
                      darkcolor=[("focus", P["accent"])])

        style.configure("TCheckbutton", background=P["bg"], foreground=P["text"],
                        focuscolor=P["bg"])
        style.map("TCheckbutton", background=[("active", P["bg"])])

        style.configure("Accent.Horizontal.TProgressbar", background=P["accent"],
                        troughcolor=P["subtle"], borderwidth=0, thickness=8)

        # card container
        style.configure("Card.TLabelframe", background=P["bg"], bordercolor=P["border"],
                        relief="solid", borderwidth=1, padding=14)
        style.configure("Card.TLabelframe.Label", background=P["bg"],
                        foreground=P["muted"], font=("Segoe UI", 10, "bold"))
        style.configure("TSeparator", background=P["border"])

        self.root.configure(background=P["bg"])

    def _style_text(self, widget):
        """Apply the modern flat look to a tk.Text widget."""
        P = self.PALETTE
        widget.configure(
            relief="solid", borderwidth=1, highlightthickness=1,
            highlightbackground=P["border"], highlightcolor=P["accent"],
            background="#ffffff", foreground=P["text"], insertbackground=P["text"],
            padx=8, pady=6, font=("Consolas", 10), selectbackground="#dbe4ff",
        )

    # --- layout -------------------------------------------------------------

    def _make_scrollable(self, parent, width):
        """Return an inner frame that scrolls vertically inside ``parent``.

        Fixes the problem where the lower options in a tall input column get
        clipped: the whole column now scrolls with the mouse wheel / scrollbar.
        """
        tk, ttk = self.tk, self.ttk
        P = self.PALETTE
        container = ttk.Frame(parent)
        container.pack(side="left", fill="y", padx=(0, 16))

        canvas = tk.Canvas(container, width=width, background=P["bg"],
                           borderwidth=0, highlightthickness=0)
        vsb = ttk.Scrollbar(container, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vsb.set)
        vsb.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)

        inner = ttk.Frame(canvas)
        win = canvas.create_window((0, 0), window=inner, anchor="nw")

        def _on_inner_configure(_):
            canvas.configure(scrollregion=canvas.bbox("all"))
        inner.bind("<Configure>", _on_inner_configure)
        # keep the inner frame the same width as the canvas
        canvas.bind("<Configure>", lambda e: canvas.itemconfigure(win, width=e.width))

        # mouse-wheel scrolling, active only while the pointer is over this column
        def _on_wheel(event):
            delta = -1 * (event.delta // 120) if event.delta else 0
            canvas.yview_scroll(delta, "units")
        canvas.bind("<Enter>", lambda e: canvas.bind_all("<MouseWheel>", _on_wheel))
        canvas.bind("<Leave>", lambda e: canvas.unbind_all("<MouseWheel>"))
        # Linux uses Button-4/5 for the wheel
        canvas.bind("<Enter>", lambda e: (canvas.bind_all("<MouseWheel>", _on_wheel),
                                          canvas.bind_all("<Button-4>", lambda ev: canvas.yview_scroll(-1, "units")),
                                          canvas.bind_all("<Button-5>", lambda ev: canvas.yview_scroll(1, "units"))))
        return inner

    def _build_widgets(self):
        tk, ttk = self.tk, self.ttk
        P = self.PALETTE

        # --- header ---
        header = ttk.Frame(self.root, padding=(20, 16, 20, 12))
        header.pack(fill="x")
        ttk.Label(header, text="🧬  ampliscan", style="Title.TLabel").pack(anchor="w")
        ttk.Label(header, text="Anchor-based amplicon demultiplexer — sort reads "
                  "into per-sample bins by their barcodes.",
                  style="Muted.TLabel").pack(anchor="w", pady=(2, 0))
        ttk.Separator(self.root).pack(fill="x")

        outer = ttk.Frame(self.root, padding=18)
        outer.pack(fill="both", expand=True)

        # two columns: inputs (left, scrollable), results (right)
        left = self._make_scrollable(outer, width=440)
        right = ttk.Frame(outer)
        right.pack(side="left", fill="both", expand=True)

        # --- Panel card ---
        pf = ttk.Labelframe(left, text="  Panel  ", style="Card.TLabelframe")
        pf.pack(fill="x")

        arow = ttk.Frame(pf); arow.pack(fill="x")
        ttk.Label(arow, text="5′ anchor", style="Field.TLabel").grid(row=0, column=0, sticky="w")
        ttk.Label(arow, text="3′ anchor", style="Field.TLabel").grid(row=0, column=1, sticky="w", padx=(10, 0))
        ttk.Label(arow, text="BC length", style="Field.TLabel").grid(row=0, column=2, sticky="w", padx=(10, 0))
        self.a5 = tk.StringVar(value="GCTT")
        self.a3 = tk.StringVar(value="ACAG")
        self.bclen = tk.StringVar(value="4")
        ttk.Entry(arow, textvariable=self.a5, width=9).grid(row=1, column=0, sticky="w", pady=(2, 0))
        ttk.Entry(arow, textvariable=self.a3, width=9).grid(row=1, column=1, sticky="w", padx=(10, 0), pady=(2, 0))
        ttk.Entry(arow, textvariable=self.bclen, width=7).grid(row=1, column=2, sticky="w", padx=(10, 0), pady=(2, 0))

        ttk.Label(pf, text="Forward barcodes   (name  sequence, one per line)",
                  style="Field.TLabel").pack(anchor="w", pady=(12, 2))
        self.fbc = tk.Text(pf, height=6, width=32)
        self._style_text(self.fbc)
        self.fbc.pack(fill="x")
        self.fbc.insert("1.0", "F01 GCGT\nF02 GTAG\nF03 ACGC\nF04 CTCG\nF05 GCTC")

        ttk.Label(pf, text="Reverse barcodes", style="Field.TLabel").pack(anchor="w", pady=(10, 2))
        self.rbc = tk.Text(pf, height=6, width=32)
        self._style_text(self.rbc)
        self.rbc.pack(fill="x")
        self.rbc.insert("1.0", "R01 ACGC\nR02 CTAC\nR03 GCGT\nR04 CGAG\nR05 GAGC")

        pbtn = ttk.Frame(pf); pbtn.pack(fill="x", pady=(10, 0))
        ttk.Button(pbtn, text="Load YAML…", style="Small.TButton", command=self._load_yaml).pack(side="left")
        ttk.Button(pbtn, text="From Excel…", style="Small.TButton", command=self._load_excel).pack(side="left", padx=(6, 0))
        ttk.Button(pbtn, text="Save YAML…", style="Small.TButton", command=self._save_yaml).pack(side="left", padx=(6, 0))

        # --- Reads card ---
        rf = ttk.Labelframe(left, text="  Reads & options  ", style="Card.TLabelframe")
        rf.pack(fill="x", pady=(14, 0))

        self.r1_path = tk.StringVar()
        self.r2_path = tk.StringVar()
        self.out_path = tk.StringVar()

        self._file_row(rf, "R1 file", self.r1_path, self._pick_r1)
        self._file_row(rf, "R2 (optional)", self.r2_path, self._pick_r2)
        self._file_row(rf, "Output folder", self.out_path, self._pick_out, folder=True)

        mrow = ttk.Frame(rf); mrow.pack(fill="x", pady=(12, 0))
        ttk.Label(mrow, text="Anchor mm", style="Field.TLabel").grid(row=0, column=0, sticky="w")
        ttk.Label(mrow, text="Barcode mm", style="Field.TLabel").grid(row=0, column=1, sticky="w", padx=(12, 0))
        self.amm = tk.StringVar(value="0")
        self.bmm = tk.StringVar(value="1")
        self.indels = tk.BooleanVar(value=False)
        ttk.Spinbox(mrow, from_=0, to=3, textvariable=self.amm, width=5).grid(row=1, column=0, sticky="w", pady=(2, 0))
        ttk.Spinbox(mrow, from_=0, to=3, textvariable=self.bmm, width=5).grid(row=1, column=1, sticky="w", padx=(12, 0), pady=(2, 0))
        ttk.Checkbutton(mrow, text="allow indels", variable=self.indels).grid(row=1, column=2, sticky="w", padx=(16, 0))

        self.compress = tk.BooleanVar(value=True)
        ttk.Checkbutton(rf, text="gzip output (.fq.gz)", variable=self.compress).pack(anchor="w", pady=(10, 0))

        self.also_fasta = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            rf, text="also write FASTA (.fa) alongside FASTQ",
            variable=self.also_fasta,
        ).pack(anchor="w", pady=(2, 0))

        self.split_strand = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            rf, text="separate + / − strands  (F01_R01_fwd vs _rev)",
            variable=self.split_strand,
        ).pack(anchor="w", pady=(2, 0))

        # CPU cores
        import os as _os
        ncpu = _os.cpu_count() or 1
        crow = ttk.Frame(rf); crow.pack(fill="x", pady=(10, 0))
        ttk.Label(crow, text="CPU cores", style="Field.TLabel").pack(side="left")
        self.cores = tk.StringVar(value="1")
        ttk.Spinbox(crow, from_=1, to=ncpu, textvariable=self.cores, width=5).pack(side="left", padx=(8, 6))
        ttk.Label(crow, text=f"of {ncpu}  (1 is usually fastest for gzipped input)",
                  style="Muted.TLabel").pack(side="left")

        self.run_btn = ttk.Button(rf, text="▶   Run demux", style="Accent.TButton", command=self._run)
        self.run_btn.pack(fill="x", pady=(16, 4))

        # --- Results (right) ---
        self.status = tk.StringVar(value="Ready — configure the panel, pick your reads, and run.")
        ttk.Label(right, textvariable=self.status, wraplength=560,
                  style="Status.TLabel").pack(anchor="w")

        self.prog = ttk.Progressbar(right, mode="determinate", maximum=100,
                                    style="Accent.Horizontal.TProgressbar")
        self.prog.pack(fill="x", pady=(10, 12))

        statf = ttk.Frame(right); statf.pack(fill="x")
        self.stat_total = tk.StringVar(value="—")
        self.stat_assigned = tk.StringVar(value="—")
        self.stat_unassigned = tk.StringVar(value="—")
        self._stat_value_labels = []
        for i, (lbl, var, colour) in enumerate([
            ("total reads", self.stat_total, P["text"]),
            ("assigned", self.stat_assigned, P["success"]),
            ("unassigned", self.stat_unassigned, P["muted"]),
        ]):
            cell = tk.Frame(statf, background=P["subtle"], highlightthickness=0)
            cell.grid(row=0, column=i, sticky="ew", padx=(0 if i == 0 else 10, 0), ipady=10)
            val = tk.Label(cell, textvariable=var, background=P["subtle"],
                           foreground=colour, font=("Segoe UI", 17, "bold"))
            val.pack(pady=(6, 0))
            self._stat_value_labels.append(val)
            tk.Label(cell, text=lbl, background=P["subtle"], foreground=P["muted"],
                     font=("Segoe UI", 9)).pack(pady=(0, 4))
            statf.columnconfigure(i, weight=1, uniform="stat")

        ttk.Label(right, text="READS PER F × R BIN", style="Section.TLabel").pack(anchor="w", pady=(18, 6))
        self._last_bins = None
        # small requested width so the canvas doesn't force the frame wider than
        # the window; fill/expand lets it grow into the available space.
        self.canvas = tk.Canvas(right, width=360, height=250, background="#ffffff",
                                highlightthickness=1, highlightbackground=P["border"])
        self.canvas.pack(fill="both", expand=True)
        # redraw the heatmap whenever the canvas is resized (or first mapped)
        self.canvas.bind("<Configure>",
                         lambda e: self._last_bins and self._draw_heatmap(self._last_bins))

        ttk.Label(right, text="UNASSIGNED REASONS", style="Section.TLabel").pack(anchor="w", pady=(14, 6))
        self.reasons = tk.Text(right, height=5)
        self._style_text(self.reasons)
        self.reasons.configure(background=P["subtle"], highlightthickness=0, borderwidth=0)
        self.reasons.pack(fill="x")
        self.reasons.configure(state="disabled")

        self.open_btn = ttk.Button(right, text="📂  Open output folder",
                                   command=self._open_output, state="disabled")
        self.open_btn.pack(anchor="e", pady=(12, 0))

    def _file_row(self, parent, label, var, cmd, folder=False):
        ttk, tk = self.ttk, self.tk
        row = ttk.Frame(parent); row.pack(fill="x", pady=(6, 0))
        ttk.Label(row, text=label, width=14, style="Field.TLabel").pack(side="left")
        ttk.Entry(row, textvariable=var).pack(side="left", fill="x", expand=True)
        ttk.Button(row, text="Browse", style="Small.TButton", command=cmd).pack(side="left", padx=(6, 0))

    # --- file dialogs -------------------------------------------------------

    def _pick_r1(self):
        self._pick_file(self.r1_path)

    def _pick_r2(self):
        self._pick_file(self.r2_path)

    def _pick_file(self, var):
        from tkinter import filedialog
        p = filedialog.askopenfilename(
            title="Select FASTQ/FASTA",
            filetypes=[("Sequence files", "*.fq *.fastq *.fa *.fasta *.gz"), ("All files", "*.*")],
        )
        if p:
            var.set(p)

    def _pick_out(self):
        from tkinter import filedialog
        p = filedialog.askdirectory(title="Select output folder")
        if p:
            self.out_path.set(p)

    # --- panel load/save ----------------------------------------------------

    def _panel_from_form(self) -> Panel:
        return Panel(
            name="gui_panel",
            forward_5p_anchor=self.a5.get().strip(),
            forward_3p_anchor=self.a3.get().strip(),
            barcode_length=int(self.bclen.get()),
            forward_barcodes=_parse_barcodes(self.fbc.get("1.0", "end")),
            reverse_barcodes=_parse_barcodes(self.rbc.get("1.0", "end")),
        )

    def _fill_form(self, panel: Panel):
        self.a5.set(panel.forward_5p_anchor)
        self.a3.set(panel.forward_3p_anchor)
        self.bclen.set(str(panel.barcode_length))
        self.fbc.delete("1.0", "end")
        self.fbc.insert("1.0", "\n".join(f"{k} {v}" for k, v in panel.forward_barcodes.items()))
        self.rbc.delete("1.0", "end")
        self.rbc.insert("1.0", "\n".join(f"{k} {v}" for k, v in panel.reverse_barcodes.items()))

    def _load_yaml(self):
        from tkinter import filedialog, messagebox
        p = filedialog.askopenfilename(filetypes=[("YAML", "*.yaml *.yml"), ("All", "*.*")])
        if not p:
            return
        try:
            self._fill_form(Panel.from_yaml(p))
            self.status.set(f"Loaded panel: {Path(p).name}")
        except Exception as e:  # noqa: BLE001
            messagebox.showerror("Load failed", str(e))

    def _load_excel(self):
        from tkinter import filedialog, messagebox, simpledialog
        p = filedialog.askopenfilename(filetypes=[("Excel", "*.xlsx *.xls"), ("All", "*.*")])
        if not p:
            return
        try:
            from ampliscan import panel_from_excel
            blen = int(self.bclen.get() or "4")
            self._fill_form(panel_from_excel(p, barcode_length=blen))
            self.status.set(f"Built panel from Excel: {Path(p).name}")
        except Exception as e:  # noqa: BLE001
            messagebox.showerror("Excel parse failed",
                                 f"{e}\n\n(Needs the 'excel' extra: pip install ampliscan[excel])")

    def _save_yaml(self):
        from tkinter import filedialog, messagebox
        try:
            panel = self._panel_from_form()
        except Exception as e:  # noqa: BLE001
            messagebox.showerror("Invalid panel", str(e)); return
        p = filedialog.asksaveasfilename(defaultextension=".yaml",
                                         filetypes=[("YAML", "*.yaml")])
        if p:
            panel.to_yaml(p)
            self.status.set(f"Saved panel: {Path(p).name}")

    # --- run ----------------------------------------------------------------

    def _run(self):
        from tkinter import messagebox
        if self.worker and self.worker.is_alive():
            return
        try:
            panel = self._panel_from_form()
        except Exception as e:  # noqa: BLE001
            messagebox.showerror("Invalid panel", str(e)); return
        if not self.r1_path.get():
            messagebox.showwarning("Missing R1", "Choose an R1 file."); return
        if not self.out_path.get():
            messagebox.showwarning("Missing output", "Choose an output folder."); return

        policy = MatchPolicy(
            anchor_mismatch=int(self.amm.get()),
            barcode_mismatch=int(self.bmm.get()),
            allow_indels=bool(self.indels.get()),
        )
        fmts = ["fastq"] + (["fasta"] if self.also_fasta.get() else [])
        try:
            workers = max(1, int(self.cores.get()))
        except (ValueError, TypeError):
            workers = 1
        args = dict(
            r1=self.r1_path.get(),
            r2=self.r2_path.get() or None,
            out=self.out_path.get(),
            panel=panel, policy=policy,
            fmt=fmts, compress=bool(self.compress.get()),
            split_by_strand=bool(self.split_strand.get()),
            workers=workers,
        )
        self.run_btn.configure(state="disabled")
        self.open_btn.configure(state="disabled")
        self.prog.configure(mode="indeterminate"); self.prog.start(12)
        self.status.set("Running demux… (large files take a while)")
        self.worker = threading.Thread(target=self._worker, kwargs=args, daemon=True)
        self.worker.start()

    def _worker(self, r1, r2, out, panel, policy, fmt, compress, split_by_strand, workers):
        try:
            from ampliscan import demux_to_bins
            last = [0]

            def progress(n):
                # throttle UI messages
                if n - last[0] >= 100000:
                    last[0] = n
                    self.q.put(("progress", n))

            stats = demux_to_bins(
                r1, panel, out, r2_path=r2, policy=policy,
                workers=workers, output_format=fmt, compress=compress,
                split_by_strand=split_by_strand, progress=progress,
            )
            self.q.put(("done", stats))
        except Exception:  # noqa: BLE001
            self.q.put(("error", traceback.format_exc()))

    # --- queue polling / UI updates ----------------------------------------

    def _poll_queue(self):
        try:
            while True:
                kind, payload = self.q.get_nowait()
                if kind == "progress":
                    self.status.set(f"Running demux… {payload:,} reads processed")
                elif kind == "done":
                    self._show_results(payload)
                elif kind == "error":
                    from tkinter import messagebox
                    self.prog.stop(); self.prog.configure(mode="determinate", value=0)
                    self.run_btn.configure(state="normal")
                    self.status.set("Failed.")
                    messagebox.showerror("Demux failed", payload)
        except queue.Empty:
            pass
        self.root.after(120, self._poll_queue)

    def _show_results(self, stats):
        self.prog.stop()
        self.prog.configure(mode="determinate", value=100)
        self.run_btn.configure(state="normal")
        self.open_btn.configure(state="normal")
        pct = stats.assigned_fraction * 100
        self.status.set(f"Done — {stats.total:,} reads, {pct:.1f}% assigned.")
        self.stat_total.set(f"{stats.total:,}")
        self.stat_assigned.set(f"{pct:.1f}%")
        self.stat_unassigned.set(f"{stats.unassigned:,}")

        self._draw_heatmap(stats.bin_counts)

        self.reasons.configure(state="normal")
        self.reasons.delete("1.0", "end")
        if stats.reason_counts:
            for k, v in stats.reason_counts.most_common():
                self.reasons.insert("end", f"  {k:28s} {v:,}\n")
        else:
            self.reasons.insert("end", "  none — everything assigned\n")
        self.reasons.configure(state="disabled")

    def _draw_heatmap(self, bin_counts):
        # Remember the data so we can redraw responsively on canvas resize.
        self._last_bins = bin_counts
        c = self.canvas
        c.delete("all")
        if not bin_counts:
            return
        f_names, r_names = set(), set()
        for k in bin_counts:
            f, _, r = k.partition("_")
            f_names.add(f); r_names.add(r)
        f_names = sorted(f_names)
        r_names = sorted(r_names)
        mx = max(bin_counts.values()) or 1

        c.update_idletasks()
        W = c.winfo_width()
        H = c.winfo_height()
        if W < 60 or H < 60:   # not laid out yet; a <Configure> redraw will follow
            return
        pad_l, pad_t = 44, 24
        # cap cell size so a few columns don't stretch into giant blocks
        gw = min(96, max(28, (W - pad_l - 12) // max(1, len(r_names))))
        gh = min(70, max(24, (H - pad_t - 12) // max(1, len(f_names))))

        muted = self.PALETTE["muted"]
        for ci, r in enumerate(r_names):
            x = pad_l + ci * gw
            c.create_text(x + gw / 2, pad_t / 2, text=r, fill=muted,
                          font=("Segoe UI", 8, "bold"))
        for ri, f in enumerate(f_names):
            y = pad_t + ri * gh
            c.create_text(pad_l / 2, y + gh / 2, text=f, fill=muted,
                          font=("Segoe UI", 8, "bold"))
            for ci, r in enumerate(r_names):
                v = bin_counts.get(f"{f}_{r}", 0)
                x = pad_l + ci * gw
                colour = _viridis(v / mx)
                c.create_rectangle(x + 1, y + 1, x + gw - 3, y + gh - 3,
                                   fill=colour, outline="")
                # bright (high) cells get dark text, dark cells get light text
                fg = "#111827" if v / mx >= 0.55 else "#ffffff"
                c.create_text(x + gw / 2, y + gh / 2, text=_fmt_count(v),
                              fill=fg, font=("Segoe UI", 8))

    def _open_output(self):
        import os
        import subprocess
        import sys
        path = self.out_path.get()
        if not path or not Path(path).exists():
            return
        try:
            if sys.platform == "win32":
                os.startfile(path)  # noqa: S606
            elif sys.platform == "darwin":
                subprocess.run(["open", path], check=False)
            else:
                subprocess.run(["xdg-open", path], check=False)
        except Exception:  # noqa: BLE001
            pass


def launch() -> int:
    """Entry point for ``ampliscan gui`` / ``ampliscan-gui``."""
    if not _require_tk():
        print(
            "The GUI needs tkinter, which is missing from this Python build.\n"
            "  - Windows/macOS: reinstall Python from python.org (tkinter is included)\n"
            "  - Debian/Ubuntu: sudo apt install python3-tk\n"
            "You can still use the terminal interface: run  ampliscan demux --help",
        )
        return 1
    import tkinter as tk
    root = tk.Tk()
    _set_window_icon(root)
    AmpliscanGUI(root)  # applies its own modern theme in _setup_style
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(launch())
