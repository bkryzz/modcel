# Specification Review

Source: `model.tex`

## Coherent Core

The general workflow is internally consistent:

1. Define species in a top field.
2. Confirm species to populate dropdown menus.
3. Add relation rows with:
   - `FIELD1`: source species
   - `FIELD2`: verbal action
   - `FIELD3`: target species
   - `FIELD4`: kinetic law
   - `FIELD5`: optional third partner
   - `GATE`: add multiplicative regulation
   - `COMMIT`: accept the current relation
4. Convert committed flux terms into species ODEs.
5. GATE rows should multiply an existing flux and should not be counted as
   independent contributions.

## Ambiguities And Corrections Used In This Prototype

- `FIELD5` is described as optional but is essential for the `PRODUCES / MASS ACTION`
  example, where `A + B -> C`. This prototype treats `FIELD1` and `FIELD3` as
  reactants and `FIELD5` as the product for that specific case.
- `CONSUMES / LINEAR` states `v = kAB`, but the prompt says `v = kB` and the
  contributions keep `A` unchanged. This prototype uses `v = kB`.
- `CONSUMES / MASS ACTION` states `v = kAB`, but the prompt says `v = kB`.
  This prototype uses the equation text, `v = kAB`.
- `ACTIVATES / ACTIVATION (MM)` has flux depending on `A`, but the prompt says
  `B/(K+B)`. This prototype uses the equation text, `A/(K+A)`.
- `INHIBITION (HILL)` defines `K^n/(K^n + A^n)`, but the prompt text says
  `A^n/(K^n+A^n)`. This prototype uses the inhibitory equation,
  `K^n/(K^n + A^n)`.
- `INHIBITION (COMP)` asks for `V,a` but says to read three fields. This
  prototype reads two fields.
- Inhibition entries say they have no independent contribution, but show
  derivative equations using `v * f`. This prototype treats all GATE entries as
  multipliers of the parent flux and applies only the parent flux stoichiometry.

## First-Version Scope

This version implements the dropdown-based workflow, parameter prompts, GATE
multipliers, committed flux collection, and generated ODE preview. It keeps the
older symbolic parser available as a separate module for later near-human text
input.
