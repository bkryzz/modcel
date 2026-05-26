"""Kinetic catalog for the dropdown-based metabolic interpreter."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class Action(str, Enum):
    PRODUCES = "PRODUCES"
    CONSUMES = "CONSUMES"
    ACTIVATES = "ACTIVATES"


class Kinetic(str, Enum):
    LINEAR = "LINEAR"
    MASS_ACTION = "MASS ACTION"
    ACTIVATION_MM = "ACTIVATION (MM)"
    ACTIVATION_HILL = "ACTIVATION (HILL)"
    INHIBITION_HILL = "INHIBITION (HILL)"
    INHIBITION_COMP = "INHIBITION (COMP)"
    REVERSIBLE = "REVERSIBLE"


@dataclass(frozen=True)
class KineticRule:
    action: Action
    kinetic: Kinetic
    params: tuple[str, ...]
    gate_only: bool = False
    requires_partner: bool = False


RULES: tuple[KineticRule, ...] = (
    KineticRule(Action.PRODUCES, Kinetic.LINEAR, ("k",)),
    KineticRule(Action.PRODUCES, Kinetic.MASS_ACTION, ("k",), requires_partner=True),
    KineticRule(Action.PRODUCES, Kinetic.ACTIVATION_MM, ("V", "K")),
    KineticRule(Action.PRODUCES, Kinetic.ACTIVATION_HILL, ("V", "K", "n")),
    KineticRule(Action.PRODUCES, Kinetic.REVERSIBLE, ("k1", "k2")),
    KineticRule(Action.CONSUMES, Kinetic.LINEAR, ("k",)),
    KineticRule(Action.CONSUMES, Kinetic.MASS_ACTION, ("k",)),
    KineticRule(Action.CONSUMES, Kinetic.ACTIVATION_MM, ("V", "K")),
    KineticRule(Action.ACTIVATES, Kinetic.LINEAR, ("k",)),
    KineticRule(Action.ACTIVATES, Kinetic.ACTIVATION_MM, ("V", "K")),
    KineticRule(Action.ACTIVATES, Kinetic.ACTIVATION_HILL, ("V", "K", "n")),
    KineticRule(Action.ACTIVATES, Kinetic.INHIBITION_HILL, ("V", "K", "n"), gate_only=True),
    KineticRule(Action.ACTIVATES, Kinetic.INHIBITION_COMP, ("V", "a"), gate_only=True),
)


def actions() -> list[str]:
    return [action.value for action in Action]


def kinetics_for(action: str, gate: bool = False) -> list[str]:
    selected = Action(action)
    return [
        rule.kinetic.value
        for rule in RULES
        if rule.action == selected and (gate or not rule.gate_only)
    ]


def get_rule(action: str, kinetic: str) -> KineticRule:
    selected_action = Action(action)
    selected_kinetic = Kinetic(kinetic)
    for rule in RULES:
        if rule.action == selected_action and rule.kinetic == selected_kinetic:
            return rule
    raise ValueError(f"Unsupported combination: {action} / {kinetic}")


def make_expression(
    action: str,
    kinetic: str,
    source: str,
    target: str,
    partner: str,
    params: dict[str, str],
) -> str:
    """Return a math-like expression for the flux or gate multiplier."""

    rule = get_rule(action, kinetic)
    p = lambda name: params.get(name, name)
    a = source
    b = target

    if rule.kinetic == Kinetic.LINEAR:
        if rule.action == Action.CONSUMES:
            return f"{p('k')}*{b}"
        return f"{p('k')}*{a}"

    if rule.kinetic == Kinetic.MASS_ACTION:
        return f"{p('k')}*{a}*{b}"

    if rule.kinetic == Kinetic.ACTIVATION_MM:
        return f"{p('V')}*{a}/({p('K')}+{a})"

    if rule.kinetic == Kinetic.ACTIVATION_HILL:
        return f"{p('V')}*{a}^{p('n')}/({p('K')}^{p('n')}+{a}^{p('n')})"

    if rule.kinetic == Kinetic.INHIBITION_HILL:
        return f"{p('V')}*{p('K')}^{p('n')}/({p('K')}^{p('n')}+{a}^{p('n')})"

    if rule.kinetic == Kinetic.INHIBITION_COMP:
        return f"{p('V')}/(1+{p('a')}*{a})"

    if rule.kinetic == Kinetic.REVERSIBLE:
        return f"{p('k1')}*{a}-{p('k2')}*{b}"

    raise ValueError(f"Unsupported kinetic: {kinetic}")


def make_stoichiometry(
    action: str,
    kinetic: str,
    source: str,
    target: str,
    partner: str,
) -> dict[str, int]:
    """Return species coefficients for one base flux."""

    selected_action = Action(action)
    selected_kinetic = Kinetic(kinetic)

    if selected_action == Action.PRODUCES and selected_kinetic == Kinetic.MASS_ACTION:
        return {source: -1, target: -1, partner: 1}

    if selected_action == Action.PRODUCES:
        return {source: -1, target: 1}

    if selected_action == Action.CONSUMES and selected_kinetic == Kinetic.MASS_ACTION:
        return {source: -1, target: -1}

    if selected_action == Action.CONSUMES:
        return {target: -1}

    if selected_action == Action.ACTIVATES:
        return {target: 1}

    raise ValueError(f"Unsupported action: {action}")


def pretty_expression(expression: str) -> str:
    return expression.replace("*", " ")
