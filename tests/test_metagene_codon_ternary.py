from __future__ import annotations

import collections
import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src" / "miao_orf" / "ribo_metagene_qc.py"
SPEC = importlib.util.spec_from_file_location("test_ribo_metagene_qc", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class MetageneCodonTernaryTests(unittest.TestCase):
    def make_stats(self):
        return MODULE.TemplateStats(
            sum_counts=[250, 20, 30],
            sum_props=[0.0, 0.0, 0.0],
            patterns=collections.Counter({
                (3, 0, 0): 70,
                (0, 0, 3): 10,
                (2, 1, 0): 20,
            }),
            n_codons=100,
            n_cds_kept=5,
        )

    def test_grid_key_uses_deterministic_largest_remainders(self) -> None:
        self.assertEqual(MODULE.codon_ternary_grid_key((1, 1, 1)), (7, 7, 6))
        self.assertEqual(MODULE.codon_ternary_grid_key((1, 0, 1)), (10, 0, 10))
        self.assertEqual(MODULE.codon_ternary_grid_key((9, 1, 0)), (18, 2, 0))

    def test_summary_is_normalized_and_reports_vertex_mass(self) -> None:
        rows = MODULE.summarize_codon_ternary(self.make_stats())
        self.assertEqual(len(rows), 3)
        self.assertEqual(sum(int(row["codon_count"]) for row in rows), 100)
        self.assertAlmostEqual(sum(float(row["sample_percentage"]) for row in rows), 100.0)
        self.assertEqual(
            MODULE.codon_ternary_vertex_percentages(rows),
            (70.0, 0.0, 10.0),
        )
        mean = MODULE.codon_ternary_mean_percentages(self.make_stats())
        self.assertAlmostEqual(mean[0], 83.33333333333333)
        self.assertAlmostEqual(mean[1], 6.666666666666667)
        self.assertAlmostEqual(mean[2], 10.0)

    def test_table_and_figures_are_written(self) -> None:
        try:
            MODULE.require_dependencies()
        except SystemExit as exc:
            self.skipTest(str(exc))
        rows = MODULE.summarize_codon_ternary(self.make_stats())
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "Example"
            table = Path(f"{prefix}.codon_frame_ternary.tsv")
            MODULE.write_codon_ternary(str(table), rows)
            MODULE.plot_codon_ternary(
                str(prefix),
                rows,
                3.25,
                "success",
                MODULE.codon_ternary_mean_percentages(self.make_stats()),
            )
            with table.open("r", encoding="utf-8", newline="") as handle:
                parsed = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(parsed), 3)
            self.assertEqual(int(parsed[0]["total_codons"]), 100)
            self.assertTrue(Path(f"{prefix}.codon_frame_ternary.png").stat().st_size > 0)
            self.assertTrue(Path(f"{prefix}.codon_frame_ternary.pdf").stat().st_size > 0)


if __name__ == "__main__":
    unittest.main()
