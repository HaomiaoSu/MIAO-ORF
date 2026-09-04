#!/usr/bin/env python3
"""Regression tests for DM model-expected component read summaries."""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


CALLER_PATH = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "miao_orf"
    / "ribo_intorf_dm_caller.py"
)
SPEC = importlib.util.spec_from_file_location("dm_caller_under_test", CALLER_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Cannot load caller module from {CALLER_PATH}")
DM = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = DM
SPEC.loader.exec_module(DM)
DM.require_worker_dependencies()


class ModelExpectedReadsTest(unittest.TestCase):
    def test_fixed_importance_reports_actual_draw_count(self) -> None:
        counts = DM.np.asarray([[2, 1, 7], [1, 1, 8]], dtype=int)
        host = DM.np.asarray([0.1, 0.2, 0.7])
        trans = DM.np.asarray([0.8, 0.1, 0.1])
        result = DM.importance_sampling_bootstrap_p(
            counts,
            host,
            trans,
            20.0,
            observed_stat=0.5,
            observed_lambda=0.25,
            reps=7,
            rng=DM.np.random.default_rng(1),
            eta_values=[0.0, 0.25, 1.0],
            grid_size=11,
        )
        self.assertEqual(result["importance_total_reps_used"], 7)
        self.assertEqual(
            sum(int(value) for value in result["importance_component_draws"].split(",")),
            7,
        )

    def test_exact_zero_lrt_reports_no_simulated_draws(self) -> None:
        counts = DM.np.asarray([[2, 1, 7]], dtype=int)
        result = DM.importance_sampling_bootstrap_p(
            counts,
            DM.np.asarray([0.1, 0.2, 0.7]),
            DM.np.asarray([0.8, 0.1, 0.1]),
            20.0,
            observed_stat=0.0,
            observed_lambda=0.0,
            reps=7,
            rng=DM.np.random.default_rng(1),
            eta_values=[0.0, 1.0],
            grid_size=11,
        )
        self.assertEqual(result["importance_total_reps_used"], 0)
        self.assertEqual(result["importance_component_draws"], "exact_lrt_zero")

    def test_component_and_phase_totals(self) -> None:
        host = DM.np.asarray([0.1, 0.2, 0.7])
        trans = DM.np.asarray([0.8, 0.1, 0.1])
        result = DM.model_expected_component_reads(
            100,
            0.25,
            host,
            trans,
        )
        self.assertAlmostEqual(result["model_expected_host_cds_core_reads"], 75.0)
        self.assertAlmostEqual(result["model_expected_intorf_core_reads"], 25.0)
        self.assertAlmostEqual(
            sum(
                result[f"model_expected_host_cds_phase{i}_core_reads"]
                for i in range(3)
            ),
            75.0,
        )
        self.assertAlmostEqual(
            sum(
                result[f"model_expected_intorf_phase{i}_core_reads"]
                for i in range(3)
            ),
            25.0,
        )
        self.assertAlmostEqual(result["model_expected_intorf_phase0_core_reads"], 20.0)
        fitted_phase_counts = 100.0 * DM.pi_mix(0.25, host, trans)
        for phase in range(3):
            reconstructed = (
                result[f"model_expected_host_cds_phase{phase}_core_reads"]
                + result[f"model_expected_intorf_phase{phase}_core_reads"]
            )
            self.assertAlmostEqual(reconstructed, fitted_phase_counts[phase])

    def test_boundary_mixtures(self) -> None:
        host = DM.np.asarray([0.2, 0.3, 0.5])
        trans = DM.np.asarray([0.7, 0.2, 0.1])
        host_only = DM.model_expected_component_reads(17, 0.0, host, trans)
        intorf_only = DM.model_expected_component_reads(17, 1.0, host, trans)
        self.assertAlmostEqual(host_only["model_expected_host_cds_core_reads"], 17.0)
        self.assertAlmostEqual(host_only["model_expected_intorf_core_reads"], 0.0)
        self.assertAlmostEqual(intorf_only["model_expected_host_cds_core_reads"], 0.0)
        self.assertAlmostEqual(intorf_only["model_expected_intorf_core_reads"], 17.0)

    def test_templates_are_normalized_for_rounding(self) -> None:
        result = DM.model_expected_component_reads(
            10,
            0.4,
            DM.np.asarray([2.0, 3.0, 5.0]),
            DM.np.asarray([7.0, 2.0, 1.0]),
        )
        self.assertAlmostEqual(
            sum(
                result[f"model_expected_host_cds_phase{i}_core_reads"]
                for i in range(3)
            ),
            6.0,
        )
        self.assertAlmostEqual(
            sum(
                result[f"model_expected_intorf_phase{i}_core_reads"]
                for i in range(3)
            ),
            4.0,
        )

    def test_invalid_inputs_fail_loudly(self) -> None:
        host = DM.np.asarray([0.2, 0.3, 0.5])
        trans = DM.np.asarray([0.7, 0.2, 0.1])
        for bad_lambda in (-0.01, 1.01, float("nan")):
            with self.subTest(lambda_hat=bad_lambda):
                with self.assertRaises(ValueError):
                    DM.model_expected_component_reads(10, bad_lambda, host, trans)
        with self.assertRaises(ValueError):
            DM.model_expected_component_reads(-1, 0.5, host, trans)
        invalid_templates = (
            DM.np.asarray([0.5, 0.5]),
            DM.np.asarray([0.5, -0.1, 0.6]),
            DM.np.asarray([0.0, 0.0, 0.0]),
            DM.np.asarray([0.5, float("nan"), 0.5]),
        )
        for bad_template in invalid_templates:
            with self.subTest(template=bad_template):
                with self.assertRaises(ValueError):
                    DM.model_expected_component_reads(10, 0.5, bad_template, trans)

    def test_output_schema_contains_expected_columns(self) -> None:
        columns = DM.result_columns([2, 3, 4])
        expected = {
            "orf_region_reads",
            "analyzed_core_reads",
            "lambda_profile_ci95_low",
            "lambda_profile_ci95_high",
            "model_expected_host_cds_fraction",
            "model_expected_intorf_fraction",
            "dm_evidence_score",
            "dm_fdr_score",
            "qc_status",
            "qc_flags",
            "filter_reason",
            "background_template_sensitivity_status",
            "model_expected_host_cds_core_reads",
            "model_expected_intorf_core_reads",
            *{
                f"model_expected_{component}_phase{phase}_core_reads"
                for component in ("host_cds", "intorf")
                for phase in range(3)
            },
            *{
                f"observed_phase{phase}_core_{suffix}"
                for phase in range(3)
                for suffix in ("reads", "fraction")
            },
        }
        self.assertTrue(expected.issubset(columns))
        self.assertEqual(len(columns), len(set(columns)))
        expected_v12 = [
            "model_expected_host_cds_core_reads",
            "model_expected_intorf_core_reads",
            "model_expected_host_cds_phase0_core_reads",
            "model_expected_host_cds_phase1_core_reads",
            "model_expected_host_cds_phase2_core_reads",
            "model_expected_intorf_phase0_core_reads",
            "model_expected_intorf_phase1_core_reads",
            "model_expected_intorf_phase2_core_reads",
        ]
        start = columns.index(expected_v12[0])
        self.assertEqual(columns[start:start + len(expected_v12)], expected_v12)

    def test_observed_read_summary_distinguishes_region_and_core(self) -> None:
        mat_all = DM.np.asarray(
            [[1, 0, 0], [2, 3, 5], [1, 1, 2], [0, 0, 4]], dtype=int
        )
        mat_core = mat_all[1:3]
        active = DM.np.asarray([True, True])
        result = DM.summarize_candidate_reads(mat_all, mat_core, active)
        self.assertEqual(result["orf_region_reads"], 19)
        self.assertEqual(result["core_reads"], 14)
        self.assertEqual(result["analyzed_core_reads"], 14)
        self.assertEqual(result["observed_phase0_core_reads"], 3)
        self.assertEqual(result["observed_phase1_core_reads"], 4)
        self.assertEqual(result["observed_phase2_core_reads"], 7)
        self.assertAlmostEqual(
            sum(result[f"observed_phase{i}_core_fraction"] for i in range(3)),
            1.0,
        )

    def test_lambda_profile_interval_contains_fit_and_is_bounded(self) -> None:
        counts = DM.np.asarray(
            [[2, 1, 7], [1, 1, 8], [2, 0, 8], [1, 2, 7], [2, 1, 7]],
            dtype=int,
        )
        host = DM.np.asarray([0.1, 0.2, 0.7])
        trans = DM.np.asarray([0.8, 0.1, 0.1])
        fitted, _, ll1, _, _ = DM.fit_lambda_lrt(counts, host, trans, 20.0)
        low, high = DM.lambda_profile_likelihood_interval(
            counts, host, trans, 20.0, fitted, ll1
        )
        self.assertLessEqual(0.0, low)
        self.assertLessEqual(low, fitted)
        self.assertLessEqual(fitted, high)
        self.assertLessEqual(high, 1.0)

    def test_qc_contract_does_not_treat_nonsignificance_as_qc_failure(self) -> None:
        row = {
            "classification": "host_only_supported",
            "classification_pre_fdr": "tested",
            "p_final": 0.5,
            "q_BH": 0.9,
            "lambda_hat": 0.1,
            "template_separation_ok": 1,
            "mixture_geometry_consistent": 1,
            "credible_active_core_breadth_ok": 1,
            "credible_target_residual_breadth_ok": 1,
        }
        DM.add_output_qc_contract(row)
        self.assertEqual(row["qc_status"], "pass")
        self.assertEqual(row["qc_flags"], "")
        self.assertEqual(row["filter_reason"], "host_only_supported")
        self.assertGreater(row["dm_evidence_score"], 0.0)

    def test_run_manifest_freezes_parameters_and_schema(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            bam = root / "sample.psite.bam"
            torf = root / "reference.torf.tsv"
            background = root / "sample.dm_background.tsv"
            output = root / "sample.intorf_dm_results.tsv"
            for path in (bam, torf, background, output):
                path.write_text("test\n", encoding="utf-8")
            args = SimpleNamespace(
                psite_bam=str(bam),
                dm_background=str(background),
                pi_method="codon_equal",
                block_sizes_parsed=[2, 3, 4],
                seed=123,
            )
            manifest_path = root / "sample.run_manifest.json"
            DM.write_run_manifest(
                str(manifest_path),
                args,
                [str(torf)],
                [str(output)],
                {"background_source": "qc_file", "background_path": str(background), "A0": 20.0},
                1,
            )
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(manifest["version"], "1.0.0")
        self.assertEqual(manifest["output_schema_version"], "1.0")
        self.assertEqual(manifest["parameters"]["seed"], 123)
        self.assertIn("model_expected_intorf_core_reads", manifest["output_schema"])
        self.assertEqual(manifest["inputs"]["torf_files"][0]["sha256_status"], "computed")

    def test_writer_emits_values_and_blanks_for_unfitted_rows(self) -> None:
        fitted = {
            "torf_id": "fitted",
            **DM.model_expected_component_reads(
                20,
                0.25,
                DM.np.asarray([0.2, 0.3, 0.5]),
                DM.np.asarray([0.7, 0.2, 0.1]),
            ),
        }
        unfitted = {"torf_id": "insufficient", "classification": "insufficient_data"}
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "results.tsv"
            DM.write_results(str(output), [fitted, unfitted], [2, 3, 4])
            with output.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(rows[0]["model_expected_host_cds_core_reads"], "15")
        self.assertEqual(rows[0]["model_expected_intorf_core_reads"], "5")
        self.assertEqual(rows[1]["model_expected_host_cds_core_reads"], "")
        self.assertEqual(rows[1]["model_expected_intorf_core_reads"], "")

    def test_analyze_candidate_populates_expected_counts(self) -> None:
        counts = DM.np.asarray(
            [[2, 1, 7], [1, 1, 8], [2, 0, 8], [1, 2, 7], [2, 1, 7]],
            dtype=int,
        )
        row = DM.ORFRow(
            raw={},
            torf_id="integration",
            gene_id="gene",
            gene_name="GENE",
            transcript_id="tx",
            chrom="chr1",
            strand="+",
            blocks=[(100, 130)],
            t_start=1,
            t_end=31,
            aa_len=10,
            peptide_len=10,
            cds_first_t=0,
            orf_biotype="intORF_altframe",
            is_annotated_cds=False,
            has_mane=True,
            prepared_mat_core=counts,
            prepared_core_indices=DM.np.arange(len(counts), dtype=int),
        )
        args = SimpleNamespace(
            importance_mode_used="not_run",
            template_separation_min=0.01,
            min_intorf_aa=1,
            exclude_start_codons=0,
            exclude_stop_codons=0,
            min_core_codons=1,
            min_active_core_codons=1,
            min_core_reads=1,
            min_credible_active_core_codons=1,
            min_credible_active_core_frac=0.0,
            min_credible_target_residual_frac=0.0,
            lambda_rel_eps=1e-8,
            lambda_abs_diff_max=1.0,
            lambda_rel_diff_max=10.0,
            lambda_grid_size=21,
            bootstrap_gate_p=0.0,
            review_gate_window=0.0,
        )
        result = DM.analyze_candidate(
            None,
            row,
            DM.np.asarray([0.8, 0.1, 0.1]),
            20.0,
            args,
            DM.np.random.default_rng(1),
        )
        self.assertIn("model_expected_host_cds_core_reads", result)
        self.assertIn("model_expected_intorf_core_reads", result)
        self.assertIn("lambda_profile_ci95_low", result)
        self.assertLessEqual(result["lambda_profile_ci95_low"], result["lambda_hat"])
        self.assertLessEqual(result["lambda_hat"], result["lambda_profile_ci95_high"])
        self.assertAlmostEqual(
            result["model_expected_host_cds_core_reads"]
            + result["model_expected_intorf_core_reads"],
            result["core_reads"],
        )


if __name__ == "__main__":
    unittest.main()
