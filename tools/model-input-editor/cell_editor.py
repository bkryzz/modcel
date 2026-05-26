#!/usr/bin/env python3
"""Editor for cell phenotypes and interaction matrix."""

from __future__ import annotations

import tkinter as tk
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from tkinter import messagebox, ttk


APP_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = APP_DIR / "default"
OUTPUT_FILE = APP_DIR / "input_cell"
PLAY_SCRIPT = APP_DIR / "play.py"
ASSIGNMENT_RE = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)\s*$")
COMMENT_MARKERS = ("!", "#")

KIND_OPTIONS = ["solid", "fluid", "ameboid"]
COLOR_OPTIONS = [
    "white",
    "bisque",
    "plum",
    "coral",
    "rosybrown",
    "lightcyan",
    "gold",
    "lightblue",
    "crimson",
    "chartreuse",
    "azure",
    "darkblue",
    "gray",
]
CELL_COLUMNS = [
    ("cell\ntype", "type", "text"),
    ("kind", "kind", "choice"),
    ("color", "color", "choice"),
    ("size\n(relative)", "size", "number"),
    ("mobility", "mobility", "number"),
    ("duplication\nprobability", "duplication", "number"),
    ("necrotic\nthreshold", "necrotic", "number"),
    ("adhesion\n(self)", "adhesion_self", "number"),
    ("adhesion\n(hetero)", "adhesion_hetero", "number"),
]

PRESETS = {
    2: {
        "type": ["hepato", "sinus", "kupfer", "star", "myofib"],
        "kind": ["solid", "fluid", "ameboid", "ameboid", "solid"],
        "color": ["rosybrown", "lightblue", "bisque", "coral", "azure"],
        "size": ["1.0", "1.0", "0.1", "0.1", "1.0"],
        "mobility": ["0.1", "0.", "0.03", "0.", "0.5"],
    },
    3: {
        "type": ["astro", "macro", "glioma", "capillary"],
        "kind": ["solid", "solid", "solid", "fluid"],
        "color": ["bisque", "plum", "gray", "crimson"],
        "size": ["1.0", "1.0", "0.6", "0.5"],
        "mobility": ["0.1", "0.1", "0.1", "0.1"],
    },
}

BASE_CELL_COUNTS = {2: 5, 3: 4}
DEFAULT_CELL_VALUES = {
    "duplication": "1.",
    "necrotic": "0.",
    "adhesion_self": "-1.",
    "adhesion_hetero": "1.",
}


@dataclass(frozen=True)
class TemplateValues:
    nr_cel: int
    imeta: int


def read_template_values(path: Path = DEFAULT_INPUT) -> TemplateValues:
    if not path.exists():
        raise ValueError("No default parameter file found")

    values: dict[str, str] = {}
    for raw in path.read_text().splitlines():
        match = ASSIGNMENT_RE.match(raw)
        if not match:
            continue
        name, value = match.groups()
        values[name.upper()] = split_comment(value).strip()

    missing = [name for name in ("NR_CEL", "IMETA") if name not in values]
    if missing:
        raise ValueError(f"Missing required variable(s) in default: {', '.join(missing)}")

    try:
        nr_cel = int(values["NR_CEL"])
        imeta = int(values["IMETA"])
    except ValueError as exc:
        raise ValueError("NR_CEL and IMETA must be integer values") from exc

    if nr_cel < 1:
        raise ValueError("NR_CEL must be at least 1")
    if imeta in BASE_CELL_COUNTS and nr_cel < BASE_CELL_COUNTS[imeta]:
        raise ValueError(f"IMETA = {imeta} requires NR_CEL to be at least {BASE_CELL_COUNTS[imeta]}")

    return TemplateValues(nr_cel=nr_cel, imeta=imeta)


def split_comment(value: str) -> str:
    quote: str | None = None
    for index, char in enumerate(value):
        if char in ("'", '"'):
            quote = None if quote == char else char if quote is None else quote
        if quote is None and char in COMMENT_MARKERS:
            return value[:index]
    return value


def preset_value(imeta: int, field: str, index: int) -> str:
    preset = PRESETS.get(imeta, {})
    values = preset.get(field, [])
    if index < len(values):
        return values[index]
    if index < BASE_CELL_COUNTS.get(imeta, 0):
        return DEFAULT_CELL_VALUES.get(field, "")
    return ""


class CellEditor(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("Define Cell Phenotypes")
        self.geometry("1180x760")
        self.minsize(900, 520)

        self.template: TemplateValues | None = None
        self.cell_vars: list[dict[str, tk.StringVar]] = []
        self.matrix_vars: list[list[tk.StringVar]] = []
        self.play_process: subprocess.Popen[str] | None = None

        self._build_shell()
        self.load_default()

    def _build_shell(self) -> None:
        top = ttk.Frame(self, padding=10)
        top.pack(fill=tk.X)

        ttk.Button(top, text="Reload default", command=self.load_default).pack(side=tk.LEFT)
        ttk.Button(top, text="Save input_cell", command=self.save).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Button(top, text="Continue", command=self.continue_to_play).pack(side=tk.LEFT, padx=(8, 0))

        self.info = ttk.Label(top, text="Loading default parameter file.")
        self.info.pack(side=tk.LEFT, padx=14)

        container = ttk.Frame(self)
        container.pack(fill=tk.BOTH, expand=True, padx=10, pady=(0, 10))

        self.canvas = tk.Canvas(container, highlightthickness=0)
        vbar = ttk.Scrollbar(container, orient=tk.VERTICAL, command=self.canvas.yview)
        hbar = ttk.Scrollbar(container, orient=tk.HORIZONTAL, command=self.canvas.xview)
        self.form = ttk.Frame(self.canvas, padding=(0, 0, 10, 10))

        self.form.bind("<Configure>", lambda event: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        self.canvas.create_window((0, 0), window=self.form, anchor="nw")
        self.canvas.configure(yscrollcommand=vbar.set, xscrollcommand=hbar.set)

        self.canvas.grid(row=0, column=0, sticky="nsew")
        vbar.grid(row=0, column=1, sticky="ns")
        hbar.grid(row=1, column=0, sticky="ew")
        container.columnconfigure(0, weight=1)
        container.rowconfigure(0, weight=1)

        self.status = ttk.Label(self, text="", padding=(10, 0, 10, 10))
        self.status.pack(fill=tk.X)

    def load_default(self) -> None:
        try:
            self.template = read_template_values()
        except ValueError as exc:
            self.template = None
            self.cell_vars.clear()
            self.matrix_vars.clear()
            self.info.configure(text=str(exc))
            self.status.configure(text="")
            self.render()
            messagebox.showerror("Cannot create cell input", str(exc))
            return

        self.build_variables()
        self.info.configure(text=f"NR_CEL = {self.template.nr_cel}, IMETA = {self.template.imeta}")
        self.render()

    def build_variables(self) -> None:
        assert self.template is not None
        self.cell_vars = []
        for i in range(self.template.nr_cel):
            row: dict[str, tk.StringVar] = {}
            for _, key, _ in CELL_COLUMNS:
                row[key] = tk.StringVar(value=preset_value(self.template.imeta, key, i))
            self.cell_vars.append(row)

        self.matrix_vars = []
        base_count = BASE_CELL_COUNTS.get(self.template.imeta, 0)
        for i in range(self.template.nr_cel):
            row = []
            for j in range(self.template.nr_cel):
                if i < base_count and j < base_count:
                    value = "-1.0" if i == j else "1.0"
                else:
                    value = ""
                row.append(tk.StringVar(value=value))
            self.matrix_vars.append(row)

    def render(self) -> None:
        for child in self.form.winfo_children():
            child.destroy()

        if self.template is None:
            ttk.Label(self.form, text="No default parameter file found", padding=20).grid(row=0, column=0, sticky="w")
            return

        row = 0
        ttk.Label(self.form, text="Cell Characteristics", font=("", 14, "bold")).grid(
            row=row, column=0, columnspan=len(CELL_COLUMNS), sticky="w", pady=(0, 8)
        )
        row += 1

        for col, (label, _, _) in enumerate(CELL_COLUMNS):
            ttk.Label(self.form, text=label, foreground="#164a8b", justify=tk.CENTER).grid(
                row=row, column=col, sticky="ew", padx=3, pady=(0, 4)
            )
        row += 1

        for cell_index, fields in enumerate(self.cell_vars, start=1):
            for col, (_, key, kind) in enumerate(CELL_COLUMNS):
                if key == "kind":
                    widget = ttk.Combobox(self.form, textvariable=fields[key], values=KIND_OPTIONS, width=11)
                elif key == "color":
                    widget = ttk.Combobox(self.form, textvariable=fields[key], values=COLOR_OPTIONS, width=12)
                else:
                    widget = ttk.Entry(self.form, textvariable=fields[key], width=12 if kind == "number" else 14)
                widget.grid(row=row, column=col, sticky="ew", padx=3, pady=2)
            row += 1

        row += 2
        ttk.Label(self.form, text="Cell Interaction Matrix", font=("", 14, "bold")).grid(
            row=row, column=0, columnspan=self.template.nr_cel + 1, sticky="w", pady=(8, 8)
        )
        row += 1

        for col in range(self.template.nr_cel):
            ttk.Label(self.form, text=str(col + 1), foreground="#164a8b").grid(row=row, column=col + 1, padx=3)
        row += 1

        for i, matrix_row in enumerate(self.matrix_vars):
            ttk.Label(self.form, text=str(i + 1), foreground="#164a8b").grid(row=row, column=0, sticky="e", padx=(0, 6))
            for j, variable in enumerate(matrix_row):
                ttk.Entry(self.form, textvariable=variable, width=8).grid(row=row, column=j + 1, padx=3, pady=2)
            row += 1

        self.status.configure(text=f"Ready to save {OUTPUT_FILE}")

    def save(self) -> bool:
        if self.template is None:
            messagebox.showinfo("Nothing to save", "No default parameter file found.")
            return False

        try:
            self.validate_values()
            OUTPUT_FILE.write_text(self.build_output() + "\n")
        except ValueError as exc:
            messagebox.showerror("Invalid value", str(exc))
            return False
        except OSError as exc:
            messagebox.showerror("Could not save file", str(exc))
            return False

        self.status.configure(text=f"Saved {OUTPUT_FILE}")
        return True

    def continue_to_play(self) -> None:
        if not self.save():
            return
        if not PLAY_SCRIPT.exists():
            messagebox.showerror("Cannot continue", f"Could not find:\n{PLAY_SCRIPT}")
            return

        try:
            self.play_process = subprocess.Popen([sys.executable, str(PLAY_SCRIPT)])
        except OSError as exc:
            messagebox.showerror("Cannot continue", str(exc))
            return

        self.status.configure(text=f"Saved {OUTPUT_FILE} and opened {PLAY_SCRIPT.name}")
        self.after(1000, self.watch_play)

    def watch_play(self) -> None:
        if self.play_process is None:
            return
        if self.play_process.poll() is None:
            self.after(1000, self.watch_play)
            return
        self.destroy()

    def validate_values(self) -> None:
        for i, fields in enumerate(self.cell_vars, start=1):
            for label, key, kind in CELL_COLUMNS:
                value = fields[key].get().strip()
                if not value:
                    raise ValueError(f"Cell {i}: {label.replace(chr(10), ' ')} cannot be empty.")
                if kind == "number":
                    self.validate_number(value, f"Cell {i}: {label.replace(chr(10), ' ')}")

        for i, row in enumerate(self.matrix_vars, start=1):
            for j, variable in enumerate(row, start=1):
                self.validate_number(variable.get().strip(), f"Matrix value ({i}, {j})")

    @staticmethod
    def validate_number(value: str, label: str) -> None:
        try:
            float(value.replace("d", "e").replace("D", "E"))
        except ValueError as exc:
            raise ValueError(f"{label} must be numeric.") from exc

    def build_output(self) -> str:
        assert self.template is not None
        lines = [
            "# number of cell types",
            str(self.template.nr_cel),
            "# define cell characteristics",
            "#    type      kind     color     size  mobility  duplic  necros  self     potf",
        ]

        for fields in self.cell_vars:
            values = [fields[key].get().strip() for _, key, _ in CELL_COLUMNS]
            lines.append("      ".join(values))

        lines.append("# define cell interaction matrix")
        for row in self.matrix_vars:
            lines.append("      ".join(variable.get().strip() for variable in row))

        return "\n".join(lines)


if __name__ == "__main__":
    app = CellEditor()
    app.mainloop()
