import unittest

from metabolang.kinetics import get_rule, kinetics_for
from metabolang.model import Gate, Relation, build_odes
from metabolang.parser import parse_species_names


class KineticCatalogTests(unittest.TestCase):
    def test_allowed_kinetics_exclude_gate_only_entries_for_base_rows(self) -> None:
        choices = kinetics_for("ACTIVATES")
        self.assertIn("ACTIVATION (HILL)", choices)
        self.assertNotIn("INHIBITION (HILL)", choices)

    def test_produces_mass_action_requires_partner(self) -> None:
        rule = get_rule("PRODUCES", "MASS ACTION")
        self.assertTrue(rule.requires_partner)


class ModelAssemblyTests(unittest.TestCase):
    def test_linear_production_builds_expected_odes(self) -> None:
        relation = Relation(
            source="A",
            action="PRODUCES",
            target="B",
            kinetic="LINEAR",
            params={"k": "0.2"},
        )

        self.assertEqual(relation.flux_expression(), "0.2*A")
        self.assertEqual(
            build_odes(["A", "B"], [relation]),
            {
                "A": "- 0.2*A",
                "B": "0.2*A",
            },
        )

    def test_gate_multiplier_modifies_only_final_flux(self) -> None:
        relation = Relation(
            source="A",
            action="PRODUCES",
            target="B",
            kinetic="LINEAR",
            params={"k": "0.2"},
            gates=[
                Gate(
                    action="ACTIVATES",
                    kinetic="INHIBITION (COMP)",
                    source="C",
                    params={"V": "1", "a": "2"},
                )
            ],
        )

        odes = build_odes(["A", "B", "C"], [relation])
        self.assertEqual(odes["C"], "0")
        self.assertIn("(0.2*A)*(1/(1+2*C))", odes["A"])
        self.assertIn("(0.2*A)*(1/(1+2*C))", odes["B"])

    def test_produces_mass_action_uses_partner_as_product(self) -> None:
        relation = Relation(
            source="A",
            action="PRODUCES",
            target="B",
            kinetic="MASS ACTION",
            partner="C",
            params={"k": "3"},
        )

        self.assertEqual(relation.stoichiometry(), {"A": -1, "B": -1, "C": 1})
        self.assertEqual(relation.flux_expression(), "3*A*B")


class SpeciesParsingTests(unittest.TestCase):
    def test_species_names_preserve_first_seen_order(self) -> None:
        self.assertEqual(
            parse_species_names("ATP, ADP, Pi, ATP"),
            ["ATP", "ADP", "Pi"],
        )


if __name__ == "__main__":
    unittest.main()
