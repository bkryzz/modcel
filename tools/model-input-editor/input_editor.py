#!/usr/bin/env python3
"""Dynamic editor for simple VAR = value input files."""

from __future__ import annotations

import re
import subprocess
import sys
import tkinter as tk
from dataclasses import dataclass
from pathlib import Path
from tkinter import filedialog, messagebox, ttk


APP_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = APP_DIR / "default"
DEFAULT_OUTPUT = APP_DIR / "input"
CELL_EDITOR = APP_DIR / "cell_editor.py"
ASSIGNMENT_RE = re.compile(r"^(\s*)([A-Za-z_][A-Za-z0-9_]*)(\s*=\s*)(.*?)(\s*)$")
SECTION_RE = re.compile(r"^\s*#\s*([A-Za-z0-9_][A-Za-z0-9_. /]*?)\s*-{2,}\s*$")
COMMENT_MARKERS = ("!", "#")


@dataclass
class ParsedLine:
    raw: str
    number: int = 0
    indent: str = ""
    name: str | None = None
    equals: str = " = "
    value: str = ""
    trailing: str = ""
    kind: str = "text"
    section: str = "General"
    section_name: str | None = None

    @property
    def is_variable(self) -> bool:
        return self.name is not None

    @property
    def is_section(self) -> bool:
        return self.section_name is not None


def split_comment(value: str) -> tuple[str, str]:
    quote: str | None = None
    for index, char in enumerate(value):
        if char in ("'", '"'):
            quote = None if quote == char else char if quote is None else quote
        if quote is None and char in COMMENT_MARKERS:
            return value[:index].rstrip(), " " + value[index:].lstrip()
    return value.rstrip(), ""


def infer_kind(value: str) -> str:
    text = value.strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in ("'", '"'):
        return "character"
    if re.fullmatch(r"[+-]?\d+", text):
        return "integer"
    if re.fullmatch(r"[+-]?((\d+\.\d*)|(\.\d+)|(\d+))([eEdD][+-]?\d+)?", text):
        return "real"
    if text.lower() in {".true.", ".false.", "true", "false", "t", "f"}:
        return "logical"
    return "text"


def parse_input(path: Path) -> list[ParsedLine]:
    lines: list[ParsedLine] = []
    current_section = "General"
    for number, raw in enumerate(path.read_text().splitlines()):
        section_match = SECTION_RE.match(raw)
        if section_match:
            current_section = section_match.group(1).strip()
            lines.append(ParsedLine(raw=raw, number=number, section=current_section, section_name=current_section))
            continue

        match = ASSIGNMENT_RE.match(raw)
        if not match:
            lines.append(ParsedLine(raw=raw, number=number, section=current_section))
            continue

        indent, name, equals, remainder, _ = match.groups()
        value, comment = split_comment(remainder)
        kind = infer_kind(value)
        lines.append(
            ParsedLine(
                raw=raw,
                number=number,
                indent=indent,
                name=name,
                equals=equals,
                value=value.strip(),
                trailing=comment,
                kind=kind,
                section=current_section,
            )
        )
    return lines


class InputEditor(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("Input File Editor")
        self.geometry("860x620")
        self.minsize(720, 420)

        self.path: Path | None = None
        self.lines: list[ParsedLine] = []
        self.entries: dict[int, tk.StringVar] = {}
        self.cell_editor_process: subprocess.Popen[str] | None = None

        self._build_shell()
        self.open_default()

    def _build_shell(self) -> None:
        top = ttk.Frame(self, padding=10)
        top.pack(fill=tk.X)

        ttk.Button(top, text="Reload default", command=self.open_default).pack(side=tk.LEFT)
        ttk.Button(top, text="Save input", command=self.save_input).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Button(top, text="Continue", command=self.continue_to_cell_editor).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Button(top, text="Save as...", command=self.save_as).pack(side=tk.LEFT, padx=(8, 0))

        self.file_label = ttk.Label(top, text="No input file loaded")
        self.file_label.pack(side=tk.LEFT, padx=14)

        search_frame = ttk.Frame(self, padding=(10, 0, 10, 8))
        search_frame.pack(fill=tk.X)
        ttk.Label(search_frame, text="Filter").pack(side=tk.LEFT)
        self.filter_text = tk.StringVar()
        self.filter_text.trace_add("write", lambda *_: self.render_form())
        ttk.Entry(search_frame, textvariable=self.filter_text).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(8, 0))

        container = ttk.Frame(self)
        container.pack(fill=tk.BOTH, expand=True, padx=10, pady=(0, 10))

        self.canvas = tk.Canvas(container, highlightthickness=0)
        scrollbar = ttk.Scrollbar(container, orient=tk.VERTICAL, command=self.canvas.yview)
        self.form = ttk.Frame(self.canvas, padding=(0, 0, 8, 0))

        self.form.bind("<Configure>", lambda event: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        self.canvas_window = self.canvas.create_window((0, 0), window=self.form, anchor="nw")
        self.canvas.configure(yscrollcommand=scrollbar.set)
        self.canvas.bind("<Configure>", lambda event: self.canvas.itemconfigure(self.canvas_window, width=event.width))

        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        self.status = ttk.Label(self, text="Loading default input template.", padding=(10, 0, 10, 10))
        self.status.pack(fill=tk.X)

    def open_default(self) -> None:
        if not DEFAULT_INPUT.exists():
            self.lines = []
            self.entries.clear()
            self.file_label.configure(text="No default parameter file found")
            self.render_form()
            messagebox.showerror("No default parameter file found", f"Expected file:\n{DEFAULT_INPUT}")
            return

        self.load_file(DEFAULT_INPUT)

    def open_file(self) -> None:
        filename = filedialog.askopenfilename(title="Open input file")
        if not filename:
            return

        self.load_file(Path(filename))

    def load_file(self, path: Path) -> None:
        self.path = path
        try:
            self.lines = parse_input(self.path)
        except OSError as exc:
            messagebox.showerror("Could not open file", str(exc))
            return

        self.file_label.configure(text=str(self.path))
        self.entries.clear()
        for line in self.lines:
            if line.is_variable:
                self.entries[line.number] = tk.StringVar(value=line.value)
        self.render_form()

    def save_input(self) -> None:
        self.save_to(DEFAULT_OUTPUT)

    def continue_to_cell_editor(self) -> None:
        if not self.save_to(DEFAULT_OUTPUT):
            return
        if not CELL_EDITOR.exists():
            messagebox.showerror("Cannot continue", f"Could not find:\n{CELL_EDITOR}")
            return

        try:
            self.cell_editor_process = subprocess.Popen([sys.executable, str(CELL_EDITOR)])
        except OSError as exc:
            messagebox.showerror("Cannot continue", str(exc))
            return

        self.status.configure(text=f"Saved {DEFAULT_OUTPUT} and opened {CELL_EDITOR.name}")
        self.after(1000, self.watch_cell_editor)

    def watch_cell_editor(self) -> None:
        if self.cell_editor_process is None:
            return
        if self.cell_editor_process.poll() is None:
            self.after(1000, self.watch_cell_editor)
            return
        self.destroy()

    def render_form(self) -> None:
        for child in self.form.winfo_children():
            child.destroy()

        needle = self.filter_text.get().strip().lower()
        variables = [line for line in self.lines if line.is_variable and line.name]
        if needle:
            variables = [
                line
                for line in variables
                if needle in line.name.lower() or needle in line.section.lower()
            ]

        row = 0
        visible_section: str | None = None
        for line in variables:
            assert line.name is not None
            if line.section != visible_section:
                visible_section = line.section
                header = ttk.Label(
                    self.form,
                    text=visible_section,
                    font=("", 12, "bold"),
                    padding=(0, 12, 0, 4),
                )
                header.grid(row=row, column=0, columnspan=4, sticky="ew")
                row += 1

            ttk.Label(self.form, text=line.name, width=24).grid(row=row, column=0, sticky="w", pady=3)
            ttk.Label(self.form, text=line.kind, width=10).grid(row=row, column=1, sticky="w", padx=(8, 8))
            entry = ttk.Entry(self.form, textvariable=self.entries[line.number])
            entry.grid(row=row, column=2, sticky="ew", pady=3)
            if line.trailing:
                ttk.Label(self.form, text=line.trailing.strip(), foreground="#666").grid(row=row, column=3, sticky="w", padx=(8, 0))
            row += 1

        self.form.columnconfigure(2, weight=1)
        total = len([line for line in self.lines if line.is_variable])
        sections = len({line.section for line in variables})
        self.status.configure(text=f"Showing {len(variables)} of {total} variables in {sections} sections.")

    def save_as(self) -> None:
        if not self.lines:
            messagebox.showinfo("Nothing to save", "Open an input file first.")
            return

        filename = filedialog.asksaveasfilename(title="Save input file as")
        if not filename:
            return

        self.save_to(Path(filename))

    def save_to(self, path: Path) -> bool:
        if not self.lines:
            messagebox.showinfo("Nothing to save", "No default input file is loaded.")
            return False

        try:
            output = self.build_output()
            path.write_text(output + "\n")
        except ValueError as exc:
            messagebox.showerror("Invalid value", str(exc))
            return False
        except OSError as exc:
            messagebox.showerror("Could not save file", str(exc))
            return False

        self.status.configure(text=f"Saved {path}")
        return True

    def build_output(self) -> str:
        rendered: list[str] = []
        for line in self.lines:
            if not line.is_variable or not line.name:
                rendered.append(line.raw)
                continue

            value = self.entries[line.number].get().strip()
            self.validate_value(line.name, value, line.kind)
            rendered.append(f"{line.indent}{line.name}{line.equals}{value}{line.trailing}")
        return "\n".join(rendered)

    @staticmethod
    def validate_value(name: str, value: str, kind: str) -> None:
        if kind == "integer" and not re.fullmatch(r"[+-]?\d+", value):
            raise ValueError(f"{name} must be an integer.")
        if kind == "real" and not re.fullmatch(r"[+-]?((\d+\.\d*)|(\.\d+)|(\d+))([eEdD][+-]?\d+)?", value):
            raise ValueError(f"{name} must be a real number.")
        if kind == "logical" and value.lower() not in {".true.", ".false.", "true", "false", "t", "f"}:
            raise ValueError(f"{name} must be a logical value.")


if __name__ == "__main__":
    app = InputEditor()
    app.mainloop()
