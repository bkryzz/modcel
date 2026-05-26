# Metabolic Equation Interpreter

Prototype GUI for defining metabolic species, assembling kinetic relations, and
exporting the resulting ODE system as a Fortran90 `computeRates` routine.

This tool is an independent input-authoring prototype. It does not contain
project-specific simulation code or biological parameter sets.

## Run

```bash
python3 app.py
```

## Workflow

1. Define species with short acronyms separated by commas.
2. Confirm species to populate dropdown menus.
3. Define reaction/interaction rows using source, action, target, kinetics,
   optional partner, and parameters.
4. Commit rows to the generated ODE display.
5. Use `SUBMIT` to write Fortran output into `generated/`.

## Test

```bash
python3 -m unittest
```

## Notes

The files in `docs/specs/` record the design notes used to build this
prototype. The examples in this public repository are deliberately artificial.

