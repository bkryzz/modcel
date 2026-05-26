from __future__ import annotations

from pathlib import Path
import tkinter as tk
from tkinter import ttk

from metabolang.fortran_export import export_fortran_model
from metabolang.kinetics import actions, get_rule, kinetics_for, pretty_expression
from metabolang.model import Gate, Relation, build_odes
from metabolang.parser import parse_species_names


class ScrollableFrame(ttk.Frame):
    def __init__(self, parent: tk.Widget) -> None:
        super().__init__(parent)
        self.canvas = tk.Canvas(self, highlightthickness=0)
        self.scrollbar = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.content = ttk.Frame(self.canvas)
        self.window_id = self.canvas.create_window((0, 0), window=self.content, anchor="nw")

        self.canvas.configure(yscrollcommand=self.scrollbar.set)
        self.canvas.grid(row=0, column=0, sticky="nsew")
        self.scrollbar.grid(row=0, column=1, sticky="ns")

        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        self.content.bind("<Configure>", self._update_scroll_region)
        self.canvas.bind("<Configure>", self._match_content_width)
        self.canvas.bind("<Enter>", self._bind_mousewheel)
        self.canvas.bind("<Leave>", self._unbind_mousewheel)

    def _update_scroll_region(self, *_: object) -> None:
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))

    def _match_content_width(self, event: tk.Event) -> None:
        self.canvas.itemconfigure(self.window_id, width=event.width)

    def _bind_mousewheel(self, *_: object) -> None:
        self.canvas.bind_all("<MouseWheel>", self._on_mousewheel)

    def _unbind_mousewheel(self, *_: object) -> None:
        self.canvas.unbind_all("<MouseWheel>")

    def _on_mousewheel(self, event: tk.Event) -> None:
        self.canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")


class RelationEditor(ttk.Frame):
    def __init__(self, parent: tk.Widget, app: "MetabolicInterpreterApp") -> None:
        super().__init__(parent, padding=(0, 0, 0, 14))
        self.app = app
        self.gates: list[GateEditor] = []

        self.source_var = tk.StringVar()
        self.action_var = tk.StringVar(value=actions()[0])
        self.target_var = tk.StringVar()
        self.kinetic_var = tk.StringVar()
        self.partner_var = tk.StringVar()
        self.preview_var = tk.StringVar(value="")
        self.status_var = tk.StringVar(value="")
        self.param_vars: dict[str, tk.StringVar] = {}
        self.committed_index: int | None = None

        self.source = ttk.Combobox(self, textvariable=self.source_var, state="readonly", width=12)
        self.action = ttk.Combobox(self, textvariable=self.action_var, state="readonly", width=13)
        self.target = ttk.Combobox(self, textvariable=self.target_var, state="readonly", width=12)
        self.kinetic = ttk.Combobox(self, textvariable=self.kinetic_var, state="readonly", width=18)
        self.partner = ttk.Combobox(self, textvariable=self.partner_var, state="readonly", width=12)
        self.params_frame = ttk.Frame(self)
        self.gate_button = ttk.Button(self, text="GATE", command=self.add_gate)
        self.commit_button = ttk.Button(self, text="COMMIT", command=self.commit)
        self.undo_button = ttk.Button(self, text="UNDO", command=self.undo_commit, state="disabled")
        self.preview = ttk.Label(self, textvariable=self.preview_var)
        self.status = ttk.Label(self, textvariable=self.status_var, foreground="#b42318")

        self.source.grid(row=0, column=0, padx=(0, 8), sticky="ew")
        self.action.grid(row=0, column=1, padx=(0, 8), sticky="ew")
        self.target.grid(row=0, column=2, padx=(0, 8), sticky="ew")
        self.kinetic.grid(row=0, column=3, padx=(0, 8), sticky="ew")
        self.partner.grid(row=0, column=4, padx=(0, 8), sticky="ew")
        self.params_frame.grid(row=0, column=5, padx=(0, 8), sticky="w")
        self.gate_button.grid(row=0, column=6, padx=(0, 8))
        self.commit_button.grid(row=0, column=7, padx=(0, 8))
        self.undo_button.grid(row=0, column=8)
        self.preview.grid(row=1, column=0, columnspan=9, sticky="w", pady=(5, 0))
        self.status.grid(row=2, column=0, columnspan=9, sticky="w")

        self.gates_frame = ttk.Frame(self)
        self.gates_frame.grid(row=3, column=0, columnspan=9, sticky="ew", pady=(4, 0))

        self.action["values"] = actions()
        self.refresh_species()
        self._refresh_kinetics()
        self._refresh_params()
        self._bind_updates()

    def refresh_species(self) -> None:
        species = self.app.species
        for combo in (self.source, self.target, self.partner):
            combo["values"] = species
        if species:
            self.source_var.set(self.source_var.get() or species[0])
            self.target_var.set(self.target_var.get() or species[min(1, len(species) - 1)])
        for gate in self.gates:
            gate.refresh_species()
        self.refresh_preview()

    def add_gate(self) -> None:
        gate = GateEditor(self.gates_frame, self)
        gate.grid(row=len(self.gates), column=0, sticky="ew", pady=(2, 0))
        self.gates.append(gate)
        self.refresh_preview()

    def commit(self) -> None:
        relation = self.to_relation()
        if relation is None:
            return
        self.committed_index = self.app.commit_relation(relation, self.committed_index)
        self.status_var.set("Committed.")
        self._set_enabled(False)

    def undo_commit(self) -> None:
        if self.committed_index is None:
            return
        self.app.clear_relation(self.committed_index)
        self.status_var.set("Reopened for editing.")
        self._set_enabled(True)

    def to_relation(self) -> Relation | None:
        error = self._validation_error()
        if error:
            self.status_var.set(error)
            return None

        gates = []
        for gate_editor in self.gates:
            gate = gate_editor.to_gate()
            if gate is None:
                self.status_var.set("Complete all gate fields before committing.")
                return None
            gates.append(gate)

        return Relation(
            source=self.source_var.get(),
            action=self.action_var.get(),
            target=self.target_var.get(),
            kinetic=self.kinetic_var.get(),
            partner=self.partner_var.get(),
            params=self._params(),
            gates=gates,
        )

    def refresh_preview(self, *_: object) -> None:
        relation = self.to_relation_preview()
        if relation is None:
            self.preview_var.set("")
            return

        expression = relation.flux_expression()
        stoich = ", ".join(
            f"d{name}/dt {self._format_preview_coefficient(coeff)} v"
            for name, coeff in relation.stoichiometry().items()
        )
        self.preview_var.set(f"v = {pretty_expression(expression)}    |    {stoich}")
        self.status_var.set("")

    def to_relation_preview(self) -> Relation | None:
        try:
            relation = Relation(
                source=self.source_var.get(),
                action=self.action_var.get(),
                target=self.target_var.get(),
                kinetic=self.kinetic_var.get(),
                partner=self.partner_var.get(),
                params=self._params(default_symbols=True),
                gates=[
                    gate.to_gate(default_symbols=True)
                    for gate in self.gates
                    if gate.to_gate(default_symbols=True) is not None
                ],
            )
            relation.base_expression()
            return relation
        except (ValueError, KeyError):
            return None

    def _bind_updates(self) -> None:
        self.action_var.trace_add("write", self._on_action_changed)
        self.kinetic_var.trace_add("write", self._on_kinetic_changed)
        for var in (self.source_var, self.target_var, self.partner_var):
            var.trace_add("write", self.refresh_preview)

    def _on_action_changed(self, *_: object) -> None:
        self._refresh_kinetics()
        self._refresh_params()
        self.refresh_preview()

    def _on_kinetic_changed(self, *_: object) -> None:
        self._refresh_kinetics()
        self._refresh_params()
        self.refresh_preview()

    def _refresh_kinetics(self) -> None:
        choices = kinetics_for(self.action_var.get())
        self.kinetic["values"] = choices
        if self.kinetic_var.get() not in choices:
            self.kinetic_var.set(choices[0])

        try:
            rule = get_rule(self.action_var.get(), self.kinetic_var.get())
        except ValueError:
            return
        self.partner.configure(state="readonly" if rule.requires_partner else "disabled")
        if not rule.requires_partner:
            self.partner_var.set("")
        elif self.partner_var.get() not in self.app.species and self.app.species:
            self.partner_var.set(self.app.species[-1])

    def _refresh_params(self) -> None:
        try:
            rule = get_rule(self.action_var.get(), self.kinetic_var.get())
        except ValueError:
            return

        existing = self._params(default_symbols=True)
        for child in self.params_frame.winfo_children():
            child.destroy()
        self.param_vars = {}
        for index, name in enumerate(rule.params):
            var = tk.StringVar(value="" if existing.get(name) == name else existing.get(name, ""))
            var.trace_add("write", self.refresh_preview)
            self.param_vars[name] = var
            ttk.Label(self.params_frame, text=name).grid(row=0, column=index * 2, sticky="e")
            ttk.Entry(self.params_frame, textvariable=var, width=6).grid(
                row=0,
                column=index * 2 + 1,
                padx=(2, 6),
            )

    def _params(self, default_symbols: bool = False) -> dict[str, str]:
        params = {}
        for name, var in self.param_vars.items():
            value = var.get().strip()
            params[name] = value or (name if default_symbols else "")
        return params

    def _validation_error(self) -> str:
        if not self.source_var.get() or not self.target_var.get():
            return "Choose source and target species."
        try:
            rule = get_rule(self.action_var.get(), self.kinetic_var.get())
        except ValueError as error:
            return str(error)
        if rule.requires_partner and not self.partner_var.get():
            return "Choose the third partner for this kinetic law."
        missing = [name for name, value in self._params().items() if not value]
        if missing:
            return f"Define parameter values: {', '.join(missing)}."
        return ""

    def _set_enabled(self, enabled: bool) -> None:
        state = "readonly" if enabled else "disabled"
        for combo in (self.source, self.action, self.target, self.kinetic, self.partner):
            combo.configure(state=state if combo is not self.partner or combo["values"] else "disabled")
        for child in self.params_frame.winfo_children():
            if isinstance(child, ttk.Entry):
                child.configure(state="normal" if enabled else "disabled")
        self.gate_button.configure(state="normal" if enabled else "disabled")
        self.commit_button.configure(state="normal" if enabled else "disabled")
        self.undo_button.configure(state="disabled" if enabled else "normal")
        for gate in self.gates:
            gate.set_enabled(enabled)
        if enabled:
            self._refresh_kinetics()

    def _format_preview_coefficient(self, coefficient: int) -> str:
        if coefficient == 1:
            return "+"
        if coefficient == -1:
            return "-"
        return f"{'+' if coefficient > 0 else ''}{coefficient}"


class GateEditor(ttk.Frame):
    def __init__(self, parent: tk.Widget, relation: RelationEditor) -> None:
        super().__init__(parent)
        self.relation = relation
        self.source_var = tk.StringVar()
        self.action_var = tk.StringVar(value="ACTIVATES")
        self.kinetic_var = tk.StringVar()
        self.param_vars: dict[str, tk.StringVar] = {}

        ttk.Label(self, text="multiplied by").grid(row=0, column=0, padx=(28, 8))
        self.source = ttk.Combobox(self, textvariable=self.source_var, state="readonly", width=12)
        self.action = ttk.Combobox(self, textvariable=self.action_var, state="readonly", width=13)
        self.target = ttk.Label(self, text="previous flux", width=14)
        self.kinetic = ttk.Combobox(self, textvariable=self.kinetic_var, state="readonly", width=18)
        self.params_frame = ttk.Frame(self)

        self.source.grid(row=0, column=1, padx=(0, 8))
        self.action.grid(row=0, column=2, padx=(0, 8))
        self.target.grid(row=0, column=3, padx=(0, 8))
        self.kinetic.grid(row=0, column=4, padx=(0, 8))
        self.params_frame.grid(row=0, column=5, sticky="w")

        self.action["values"] = actions()
        self.refresh_species()
        self._refresh_kinetics()
        self._refresh_params()
        self._bind_updates()

    def refresh_species(self) -> None:
        self.source["values"] = self.relation.app.species
        if self.relation.app.species and not self.source_var.get():
            self.source_var.set(self.relation.app.species[0])

    def to_gate(self, default_symbols: bool = False) -> Gate | None:
        if not self.source_var.get() or not self.action_var.get() or not self.kinetic_var.get():
            return None
        params = {
            name: var.get().strip() or (name if default_symbols else "")
            for name, var in self.param_vars.items()
        }
        if not default_symbols and any(not value for value in params.values()):
            return None
        return Gate(
            action=self.action_var.get(),
            kinetic=self.kinetic_var.get(),
            source=self.source_var.get(),
            params=params,
        )

    def set_enabled(self, enabled: bool) -> None:
        state = "readonly" if enabled else "disabled"
        for combo in (self.source, self.action, self.kinetic):
            combo.configure(state=state)
        for child in self.params_frame.winfo_children():
            if isinstance(child, ttk.Entry):
                child.configure(state="normal" if enabled else "disabled")

    def _bind_updates(self) -> None:
        self.source_var.trace_add("write", self.relation.refresh_preview)
        self.action_var.trace_add("write", self._on_action_changed)
        self.kinetic_var.trace_add("write", self._on_kinetic_changed)

    def _on_action_changed(self, *_: object) -> None:
        self._refresh_kinetics()
        self._refresh_params()
        self.relation.refresh_preview()

    def _on_kinetic_changed(self, *_: object) -> None:
        self._refresh_params()
        self.relation.refresh_preview()

    def _refresh_kinetics(self) -> None:
        choices = kinetics_for(self.action_var.get(), gate=True)
        self.kinetic["values"] = choices
        if self.kinetic_var.get() not in choices:
            self.kinetic_var.set(choices[0])
        self._refresh_params()

    def _refresh_params(self) -> None:
        try:
            rule = get_rule(self.action_var.get(), self.kinetic_var.get())
        except ValueError:
            return
        existing = {
            name: var.get().strip()
            for name, var in self.param_vars.items()
        }
        for child in self.params_frame.winfo_children():
            child.destroy()
        self.param_vars = {}
        for index, name in enumerate(rule.params):
            var = tk.StringVar(value=existing.get(name, ""))
            var.trace_add("write", self.relation.refresh_preview)
            self.param_vars[name] = var
            ttk.Label(self.params_frame, text=name).grid(row=0, column=index * 2)
            ttk.Entry(self.params_frame, textvariable=var, width=6).grid(
                row=0,
                column=index * 2 + 1,
                padx=(2, 6),
            )


class MetabolicInterpreterApp(tk.Tk):
    SPECIES_PLACEHOLDER = "A, B, C, ATP, ADP, Pi"

    def __init__(self) -> None:
        super().__init__()
        self.title("Metabolic Equation Interpreter")
        self.geometry("1180x720")
        self.minsize(980, 520)

        self.species_var = tk.StringVar()
        self.species_entry: tk.Entry | None = None
        self.species_placeholder_active = False
        self.species: list[str] = []
        self.editors: list[RelationEditor] = []
        self.relations: list[Relation | None] = []
        self.odes_text: tk.Text | None = None
        self.submit_status_var = tk.StringVar(value="")

        self._build_layout()
        self._show_species_placeholder()

    def confirm_species(self) -> None:
        if self.species_placeholder_active:
            self.species = []
        else:
            self.species = parse_species_names(self.species_var.get())
        self.species_status_var.set(f"{len(self.species)} species confirmed.")
        for editor in self.editors:
            editor.refresh_species()
        if not self.editors:
            self.add_relation_editor()
        self.refresh_odes()

    def add_relation_editor(self) -> None:
        editor = RelationEditor(self.rows_frame, self)
        editor.grid(row=len(self.editors), column=0, sticky="ew")
        self.editors.append(editor)

    def commit_relation(self, relation: Relation, index: int | None = None) -> int:
        if index is None:
            self.relations.append(relation)
            index = len(self.relations) - 1
            self.add_relation_editor()
        else:
            self.relations[index] = relation
        self.refresh_odes()
        return index

    def clear_relation(self, index: int) -> None:
        self.relations[index] = None
        self.refresh_odes()

    def refresh_odes(self) -> None:
        if not self.species:
            self._set_odes_text("")
            return
        committed_relations = [relation for relation in self.relations if relation is not None]
        odes = build_odes(self.species, committed_relations)
        lines = [f"d{name}/dt = {expression}" for name, expression in odes.items()]
        self._set_odes_text("\n".join(lines))

    def submit_model(self) -> None:
        if not self.species:
            self.submit_status_var.set("Confirm species before submitting.")
            return

        committed_relations = [relation for relation in self.relations if relation is not None]
        if not committed_relations:
            self.submit_status_var.set("Commit at least one equation before submitting.")
            return

        output_dir = Path(__file__).resolve().parent / "generated"
        fortran_path, parameter_path = export_fortran_model(
            self.species,
            committed_relations,
            output_dir,
        )
        self.submit_status_var.set(
            f"Generated {fortran_path.name} and {parameter_path.name} in {output_dir}"
        )

    def _build_layout(self) -> None:
        container = ttk.Frame(self, padding=16)
        container.grid(row=0, column=0, sticky="nsew")
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        panes = tk.PanedWindow(
            container,
            orient="vertical",
            sashwidth=8,
            sashrelief="raised",
            bg="#c7ced8",
            bd=0,
        )
        panes.grid(row=0, column=0, sticky="nsew")
        container.columnconfigure(0, weight=1)
        container.rowconfigure(0, weight=1)

        input_scroll = ScrollableFrame(panes)
        equations_frame = ttk.Frame(panes, padding=(0, 12, 0, 0))
        panes.add(input_scroll, stretch="always", minsize=220)
        panes.add(equations_frame, stretch="always", minsize=140)

        top = ttk.Frame(input_scroll.content, padding=(8, 8, 16, 8))
        top.grid(row=0, column=0, sticky="nsew")
        input_scroll.content.columnconfigure(0, weight=1)
        input_scroll.content.rowconfigure(0, weight=1)

        ttk.Label(
            top,
            text="Metabolic Equation Interpreter",
            font=("Helvetica", 22, "bold"),
        ).grid(row=0, column=0, sticky="w", pady=(0, 18))

        ttk.Label(
            top,
            text="Define your metabolites with short acronyms separated by a comma",
        ).grid(row=1, column=0, sticky="w", pady=(0, 4))

        species_line = ttk.Frame(top)
        species_line.grid(row=2, column=0, sticky="ew", pady=(0, 16))
        species_line.columnconfigure(0, weight=1)

        self.species_entry = tk.Entry(
            species_line,
            textvariable=self.species_var,
            font=("Menlo", 13),
            relief="solid",
            borderwidth=1,
        )
        self.species_entry.grid(
            row=0,
            column=0,
            sticky="ew",
            padx=(0, 8),
        )
        self.species_entry.bind("<FocusIn>", self._hide_species_placeholder)
        self.species_entry.bind("<FocusOut>", self._show_species_placeholder_if_empty)
        ttk.Button(species_line, text="CONFIRM", command=self.confirm_species).grid(row=0, column=1)
        self.species_status_var = tk.StringVar(value="Confirm species to start.")
        ttk.Label(species_line, textvariable=self.species_status_var).grid(
            row=1,
            column=0,
            columnspan=2,
            sticky="w",
            pady=(4, 0),
        )

        header = ttk.Frame(top)
        header.grid(row=3, column=0, sticky="ew", pady=(0, 6))
        labels = (
            ("source", 12),
            ("action", 13),
            ("target", 12),
            ("kinetics", 18),
            ("partner", 12),
            ("parameters", 20),
        )
        for index, (label, width) in enumerate(labels):
            ttk.Label(header, text=label, width=width).grid(
                row=0,
                column=index,
                sticky="w",
                padx=(0, 8),
            )

        self.rows_frame = ttk.Frame(top)
        self.rows_frame.grid(row=4, column=0, sticky="nsew")
        self.rows_frame.columnconfigure(0, weight=1)

        top.columnconfigure(0, weight=1)
        top.rowconfigure(4, weight=1)

        equations_header = ttk.Frame(equations_frame)
        equations_header.grid(row=0, column=0, sticky="ew", pady=(0, 6))
        equations_header.columnconfigure(1, weight=1)
        ttk.Label(
            equations_header,
            text="Generated equations",
            font=("Helvetica", 14, "bold"),
        ).grid(row=0, column=0, sticky="w")
        ttk.Button(equations_header, text="SUBMIT", command=self.submit_model).grid(
            row=0,
            column=2,
            sticky="e",
        )
        equations_frame.columnconfigure(0, weight=1)
        equations_frame.rowconfigure(1, weight=1)

        text_frame = ttk.Frame(equations_frame)
        text_frame.grid(row=1, column=0, sticky="nsew")
        text_frame.columnconfigure(0, weight=1)
        text_frame.rowconfigure(0, weight=1)

        self.odes_text = tk.Text(
            text_frame,
            height=10,
            font=("Menlo", 12),
            wrap="word",
            relief="solid",
            borderwidth=1,
        )
        equations_scrollbar = ttk.Scrollbar(
            text_frame,
            orient="vertical",
            command=self.odes_text.yview,
        )
        self.odes_text.configure(yscrollcommand=equations_scrollbar.set, state="disabled")
        self.odes_text.grid(row=0, column=0, sticky="nsew")
        equations_scrollbar.grid(row=0, column=1, sticky="ns")
        ttk.Label(equations_frame, textvariable=self.submit_status_var).grid(
            row=2,
            column=0,
            sticky="w",
            pady=(6, 0),
        )

    def _show_species_placeholder(self) -> None:
        if self.species_entry is None:
            return
        self.species_placeholder_active = True
        self.species_var.set(self.SPECIES_PLACEHOLDER)
        self.species_entry.configure(fg="#8a94a6")

    def _hide_species_placeholder(self, *_: object) -> None:
        if self.species_entry is None or not self.species_placeholder_active:
            return
        self.species_placeholder_active = False
        self.species_var.set("")
        self.species_entry.configure(fg="#111827")

    def _show_species_placeholder_if_empty(self, *_: object) -> None:
        if self.species_entry is None:
            return
        if not self.species_var.get().strip():
            self._show_species_placeholder()

    def _set_odes_text(self, text: str) -> None:
        if self.odes_text is None:
            return
        self.odes_text.configure(state="normal")
        self.odes_text.delete("1.0", "end")
        self.odes_text.insert("1.0", text)
        self.odes_text.configure(state="disabled")


def main() -> None:
    app = MetabolicInterpreterApp()
    app.mainloop()


if __name__ == "__main__":
    main()
