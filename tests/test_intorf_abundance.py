#!/usr/bin/env python3
"""Tests for standalone intORF abundance quantification and orchestration."""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


ABUNDANCE = load_module(
    "quantify_intorf_abundance_test",
    ROOT / "src/miao_orf/quantify_intorf_abundance.py",
)
RUNNER = load_module("miao_orf_runner_test", ROOT / "miao_orf.py")


class IntorfAbundanceTest(unittest.TestCase):
    def test_model_allocated_pFPKM_and_interval(self) -> None:
        row = {
            "gorf_id": "g1",
            "n_core_codons": "100",
            "analyzed_core_reads": "50",
            "lambda_hat": "0.4",
            "lambda_profile_ci95_low": "0.2",
            "lambda_profile_ci95_high": "0.6",
            "model_expected_intorf_core_reads": "20",
            "model_expected_host_cds_core_reads": "30",
        }
        result = ABUNDANCE.quantify_row(row, 1_000_000)
        self.assertEqual(result["effective_core_nt"], 300)
        self.assertEqual(result["usable_library_psites"], 1_000_000)
        self.assertEqual(result["abundance_status"], "quantified")
        self.assertTrue(math.isclose(result["intorf_psite_RPM"], 20.0))
        self.assertTrue(math.isclose(result["intorf_pFPKM"], 200.0 / 3.0))
        self.assertTrue(math.isclose(result["observed_core_pFPKM"], 500.0 / 3.0))
        self.assertTrue(math.isclose(result["host_component_pFPKM"], 100.0))
        self.assertTrue(math.isclose(result["intorf_to_host_ratio"], 2.0 / 3.0))
        self.assertTrue(math.isclose(result["intorf_pFPKM_ci95_low"], 100.0 / 3.0))
        self.assertTrue(math.isclose(result["intorf_pFPKM_ci95_high"], 100.0))

    def test_unfitted_and_invalid_core_rows_are_preserved(self) -> None:
        unfitted = ABUNDANCE.quantify_row(
            {
                "gorf_id": "g2",
                "n_core_codons": "10",
                "analyzed_core_reads": "3",
                "model_expected_intorf_core_reads": "",
            },
            100,
        )
        invalid = ABUNDANCE.quantify_row(
            {"gorf_id": "g3", "n_core_codons": "0"},
            100,
        )
        self.assertEqual(unfitted["abundance_status"], "not_model_quantifiable")
        self.assertEqual(unfitted["intorf_pFPKM"], "")
        self.assertEqual(invalid["abundance_status"], "invalid_core_length")

    def test_overall_runner_places_postprocessing_between_dm_and_visualize(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            parser = RUNNER.build_parser()
            args = parser.parse_args(
                [
                    "--sample", "sample1",
                    "--out-root", temporary,
                    "--from-stage", "dm",
                    "--to-stage", "visualize",
                    "--torf", str(Path(temporary) / "reference.torf.tsv"),
                    "--psite-bam", str(Path(temporary) / "sample.psite.bam"),
                    "--dm-background", str(Path(temporary) / "background.tsv"),
                    "--export-orf-psites",
                ]
            )
            selected = RUNNER.validate_cli(args, parser)
            self.assertEqual(
                selected, ["dm", "abundance", "context", "codon", "visualize"]
            )
            paths = RUNNER.build_paths(args, selected)
            specs = RUNNER.build_stage_specs(args, paths)
            self.assertIn("quantify_intorf_abundance.py", specs["abundance"].command[1])
            self.assertEqual(specs["abundance"].required_inputs[0], paths.dm_results)
            self.assertEqual(specs["abundance"].required_inputs[1], paths.psite_bam)
            self.assertTrue(
                str(specs["abundance"].expected_outputs[0]).endswith(
                    ".intorf_abundance.tsv"
                )
            )
            self.assertIn("annotate_gene_cds_context.py", specs["context"].command[1])
            self.assertEqual(specs["context"].required_inputs, [paths.dm_results, paths.torf])
            self.assertEqual(specs["context"].expected_outputs[0], paths.gene_context_results)
            self.assertIn("export_orf_psites.py", specs["codon"].command[1])
            self.assertIn(str(paths.gene_context_results), specs["codon"].command)
            self.assertIn(str(paths.dm_run_manifest), specs["codon"].command)
            self.assertIn("--gene-level-pure-intorf-only", specs["codon"].command)
            self.assertTrue(
                str(specs["codon"].expected_outputs[0]).endswith(
                    ".orf_psite_codons.tsv.gz"
                )
            )
            self.assertEqual(specs["visualize"].required_inputs, [paths.gene_context_results])
            self.assertIn(str(paths.gene_context_results), specs["visualize"].command)
            self.assertIn("--gene-level-pure-intorf-only", specs["visualize"].command)

    def test_orf_psite_export_is_disabled_by_default(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            parser = RUNNER.build_parser()
            args = parser.parse_args(
                [
                    "--sample", "sample1",
                    "--out-root", temporary,
                    "--from-stage", "dm",
                    "--to-stage", "visualize",
                    "--torf", str(Path(temporary) / "reference.torf.tsv"),
                    "--psite-bam", str(Path(temporary) / "sample.psite.bam"),
                    "--dm-background", str(Path(temporary) / "background.tsv"),
                ]
            )
            selected = RUNNER.validate_cli(args, parser)
            self.assertEqual(selected, ["dm", "abundance", "context", "visualize"])

    def test_explicit_orf_ids_enable_export_and_override_default_filter(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            parser = RUNNER.build_parser()
            args = parser.parse_args(
                [
                    "--sample", "sample1",
                    "--out-root", temporary,
                    "--only-stage", "codon",
                    "--torf", str(Path(temporary) / "reference.torf.tsv"),
                    "--psite-bam", str(Path(temporary) / "sample.psite.bam"),
                    "--dm-run-manifest", str(Path(temporary) / "run_manifest.json"),
                    "--gene-context-results", str(Path(temporary) / "results.tsv"),
                    "--orf-psite-id", "g1",
                    "--orf-psite-id", "t2,t3",
                ]
            )
            selected = RUNNER.validate_cli(args, parser)
            self.assertEqual(selected, ["codon"])
            paths = RUNNER.build_paths(args, selected)
            command = RUNNER.build_stage_specs(args, paths)["codon"].command
            self.assertEqual(command.count("--orf-id"), 2)
            self.assertNotIn("--gene-level-pure-intorf-only", command)

    @unittest.skipUnless(importlib.util.find_spec("pysam"), "pysam is not installed")
    def test_standalone_main_counts_indexed_psite_bam(self) -> None:
        import pysam

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bam_path = root / "sample.psite.bam"
            with pysam.AlignmentFile(
                str(bam_path), "wb", header={"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
            ) as bam:
                for index, position in enumerate((10, 20, 30)):
                    read = pysam.AlignedSegment()
                    read.query_name = f"read{index}"
                    read.query_sequence = "A"
                    read.flag = 0
                    read.reference_id = 0
                    read.reference_start = position
                    read.mapping_quality = 60
                    read.cigar = ((0, 1),)
                    read.query_qualities = pysam.qualitystring_to_array("I")
                    bam.write(read)
            pysam.index(str(bam_path))

            dm_path = root / "sample.intorf_dm_results.tsv"
            columns = sorted(ABUNDANCE.REQUIRED_DM_COLUMNS)
            with dm_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
                writer.writeheader()
                writer.writerow(
                    {
                        "gorf_id": "g1",
                        "n_core_codons": "10",
                        "analyzed_core_reads": "5",
                        "lambda_hat": "0.4",
                        "lambda_profile_ci95_low": "0.2",
                        "lambda_profile_ci95_high": "0.6",
                        "model_expected_intorf_core_reads": "2",
                        "model_expected_host_cds_core_reads": "3",
                    }
                )

            prefix = root / "out" / "sample"
            self.assertEqual(
                ABUNDANCE.main(
                    [
                        "--sample", "sample",
                        "--dm-results", str(dm_path),
                        "--psite-bam", str(bam_path),
                        "--out-prefix", str(prefix),
                    ]
                ),
                0,
            )
            with Path(f"{prefix}.intorf_abundance.tsv").open(
                "r", encoding="utf-8", newline=""
            ) as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["usable_library_psites"], "3")
            self.assertTrue(math.isclose(float(row["intorf_pFPKM"]), 1.0e9 * 2 / 90))
            manifest = json.loads(
                Path(f"{prefix}.intorf_abundance_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["usable_library_psites"], 3)
            self.assertEqual(manifest["library_count_method"], "BAM_index_mapped_alignments")


if __name__ == "__main__":
    unittest.main()
