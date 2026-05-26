"""Model assembly for committed metabolic relation rows."""

from __future__ import annotations

from dataclasses import dataclass, field

from metabolang.kinetics import make_expression, make_stoichiometry


@dataclass
class Gate:
    action: str
    kinetic: str
    source: str
    params: dict[str, str]

    def expression(self) -> str:
        return make_expression(self.action, self.kinetic, self.source, "previous flux", "", self.params)


@dataclass
class Relation:
    source: str
    action: str
    target: str
    kinetic: str
    partner: str = ""
    params: dict[str, str] = field(default_factory=dict)
    gates: list[Gate] = field(default_factory=list)

    def base_expression(self) -> str:
        return make_expression(
            self.action,
            self.kinetic,
            self.source,
            self.target,
            self.partner,
            self.params,
        )

    def flux_expression(self) -> str:
        expression = self.base_expression()
        for gate in self.gates:
            expression = f"({expression})*({gate.expression()})"
        return expression

    def stoichiometry(self) -> dict[str, int]:
        return make_stoichiometry(
            self.action,
            self.kinetic,
            self.source,
            self.target,
            self.partner,
        )


def build_odes(species: list[str], relations: list[Relation]) -> dict[str, str]:
    """Build dX/dt expressions by regrouping committed flux terms."""

    terms: dict[str, list[str]] = {name: [] for name in species}
    for index, relation in enumerate(relations, start=1):
        flux_name = f"v{index}"
        flux = relation.flux_expression()
        for species_name, coefficient in relation.stoichiometry().items():
            if species_name not in terms or coefficient == 0:
                continue
            sign = "+" if coefficient > 0 else "-"
            magnitude = abs(coefficient)
            body = flux if magnitude == 1 else f"{magnitude}*({flux})"
            terms[species_name].append(f"{sign} {body}")

    return {
        species_name: " ".join(parts).lstrip("+ ").strip() or "0"
        for species_name, parts in terms.items()
    }
