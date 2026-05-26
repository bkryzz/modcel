"""Parser for a tiny metabolic reaction language.

The first language version is deliberately small:

    A + B -> C
    2 ATP + Glucose <-> 2 ADP + G6P

Supported operators:
    +      separates terms
    ->     irreversible reaction
    <->    reversible reaction

Terms may have an optional leading positive numeric coefficient.
"""

from __future__ import annotations

from dataclasses import dataclass
import re


ARROWS = ("<->", "->")


@dataclass(frozen=True)
class Term:
    """One species term in a reaction side."""

    species: str
    coefficient: float = 1.0

    def display(self) -> str:
        if self.coefficient == 1:
            return self.species
        if self.coefficient.is_integer():
            return f"{int(self.coefficient)} {self.species}"
        return f"{self.coefficient:g} {self.species}"


@dataclass(frozen=True)
class Reaction:
    """Parsed reaction."""

    reactants: tuple[Term, ...]
    products: tuple[Term, ...]
    arrow: str

    @property
    def reversible(self) -> bool:
        return self.arrow == "<->"

    def display(self) -> str:
        arrow = "⇌" if self.reversible else "→"
        return f"{format_side(self.reactants)} {arrow} {format_side(self.products)}"

    def latex(self) -> str:
        arrow = r"\rightleftharpoons" if self.reversible else r"\longrightarrow"
        return f"{format_side_latex(self.reactants)} {arrow} {format_side_latex(self.products)}"


class ReactionSyntaxError(ValueError):
    """Raised when reaction input cannot be parsed."""


def parse_species_list(text: str) -> set[str]:
    """Parse comma/newline separated species names."""

    return set(parse_species_names(text))


def parse_species_names(text: str) -> list[str]:
    """Parse comma/newline separated species names, preserving first-seen order."""

    names = re.split(r"[,\n]+", text)
    seen = set()
    parsed = []
    for name in names:
        stripped = name.strip()
        if stripped and stripped not in seen:
            seen.add(stripped)
            parsed.append(stripped)
    return parsed


def parse_reaction(text: str, allowed_species: set[str] | None = None) -> Reaction:
    """Parse a reaction string into structured terms."""

    stripped = text.strip()
    if not stripped:
        raise ReactionSyntaxError("Enter a reaction.")

    arrow = _find_arrow(stripped)
    left, right = stripped.split(arrow, maxsplit=1)
    reactants = _parse_side(left, allowed_species)
    products = _parse_side(right, allowed_species)

    if not reactants:
        raise ReactionSyntaxError("Add at least one reactant before the arrow.")
    if not products:
        raise ReactionSyntaxError("Add at least one product after the arrow.")

    return Reaction(tuple(reactants), tuple(products), arrow)


def format_side(terms: tuple[Term, ...]) -> str:
    return " + ".join(term.display() for term in terms)


def format_side_latex(terms: tuple[Term, ...]) -> str:
    return " + ".join(_term_latex(term) for term in terms)


def _find_arrow(text: str) -> str:
    arrow_pattern = re.compile(r"<->|->")
    matches = arrow_pattern.findall(text)
    if not matches:
        raise ReactionSyntaxError("Use -> or <-> to separate reactants and products.")
    if len(matches) != 1:
        raise ReactionSyntaxError("Use exactly one reaction arrow.")
    return matches[0]


def _parse_side(text: str, allowed_species: set[str] | None) -> list[Term]:
    pieces = [piece.strip() for piece in text.split("+")]
    terms = []
    for piece in pieces:
        if not piece:
            raise ReactionSyntaxError("Empty term near '+'.")
        terms.append(_parse_term(piece, allowed_species))
    return terms


def _parse_term(text: str, allowed_species: set[str] | None) -> Term:
    match = re.fullmatch(r"(?:(\d+(?:\.\d+)?)\s+)?([A-Za-z][A-Za-z0-9_]*)", text)
    if not match:
        raise ReactionSyntaxError(f"Cannot parse term: {text!r}.")

    coefficient_text, species = match.groups()
    if allowed_species is not None and species not in allowed_species:
        raise ReactionSyntaxError(f"Unknown species: {species}.")

    coefficient = float(coefficient_text) if coefficient_text else 1.0
    if coefficient <= 0:
        raise ReactionSyntaxError("Coefficients must be positive.")

    return Term(species=species, coefficient=coefficient)


def _term_latex(term: Term) -> str:
    species = re.sub(r"(\d+)", r"_{\1}", term.species)
    if term.coefficient == 1:
        return species
    if term.coefficient.is_integer():
        return f"{int(term.coefficient)}\\,{species}"
    return f"{term.coefficient:g}\\,{species}"
