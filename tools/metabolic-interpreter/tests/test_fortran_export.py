import tempfile
from pathlib import Path
import unittest

from metabolang.fortran_export import export_fortran_model
from metabolang.model import Relation


class FortranExportTests(unittest.TestCase):
    def test_export_writes_parameter_file_and_compute_rates_subroutine(self) -> None:
        relation = Relation(
            source="A",
            action="PRODUCES",
            target="B",
            kinetic="ACTIVATION (MM)",
            params={"V": "2", "K": "0.5"},
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            fortran_path, parameter_path = export_fortran_model(
                ["A", "B"],
                [relation],
                Path(temp_dir),
            )

            parameter_text = parameter_path.read_text(encoding="utf-8")
            fortran_text = fortran_path.read_text(encoding="utf-8")

        self.assertIn("V_R01 = 2", parameter_text)
        self.assertIn("KM_R01 = 0.5", parameter_text)
        self.assertIn("subroutine computeRates(time,CONSTS,RATES,STATES,k)", fortran_text)
        self.assertIn("A = STATES(1)", fortran_text)
        self.assertIn("B = STATES(2)", fortran_text)
        self.assertIn("V_R01 = CONSTS(1)", fortran_text)
        self.assertIn("KM_R01 = CONSTS(2)", fortran_text)
        self.assertIn("RATES(1) = - V_R01*A/(KM_R01+A)", fortran_text)
        self.assertIn("RATES(2) = V_R01*A/(KM_R01+A)", fortran_text)

    def test_export_converts_powers_to_fortran_syntax(self) -> None:
        relation = Relation(
            source="A",
            action="ACTIVATES",
            target="B",
            kinetic="ACTIVATION (HILL)",
            params={"V": "1", "K": "3", "n": "2"},
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            fortran_path, _parameter_path = export_fortran_model(
                ["A", "B"],
                [relation],
                Path(temp_dir),
            )

            fortran_text = fortran_path.read_text(encoding="utf-8")

        self.assertIn("A**N_R01", fortran_text)
        self.assertIn("KM_R01**N_R01", fortran_text)


if __name__ == "__main__":
    unittest.main()
