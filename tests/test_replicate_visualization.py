#!/usr/bin/env python3
"""Tests for replicate-integration visualization data preparation."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src" / "miao_orf" / "visualize_replicate_integration.py"
SPEC = importlib.util.spec_from_file_location("replicate_visualization_under_test", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class ReplicateVisualizationTest(unittest.TestCase):
    def setUp(self) -> None:
        self.replicates = ["Rep1", "Rep2", "RepUnavailable"]
        self.long_rows = [
            {
                "candidate_key": "g1|intORF_altframe", "replicate_id": "Rep1",
                "q_BH": "0.01", "lambda_hat": "0.8", "core_reads": "20",
                "primary_credible_call": "1",
            },
            {
                "candidate_key": "g1|intORF_altframe", "replicate_id": "Rep2",
                "q_BH": "0.02", "lambda_hat": "0.7", "core_reads": "18",
                "primary_credible_call": "1",
            },
            {
                "candidate_key": "g2|intORF_altframe", "replicate_id": "Rep1",
                "q_BH": "0.60", "lambda_hat": "0.2", "core_reads": "5",
                "primary_credible_call": "0",
            },
            {
                "candidate_key": "g2|intORF_altframe", "replicate_id": "Rep2",
                "q_BH": "0.03", "lambda_hat": "0.6", "core_reads": "12",
                "primary_credible_call": "1",
            },
        ]
        self.consensus_rows = [
            {
                "candidate_key": "g1|intORF_altframe", "n_primary_credible": "2",
                "gorf_id": "g|chr1|+|blocks|h111111abcdef", "overlap_type": "F1",
                "gene_id": "ENSG000001.1", "gene_name": "GENE1",
                "max_core_reads": "20",
                "replicate::Rep1::present": "1",
                "replicate::Rep1::primary_credible_call": "1",
                "replicate::Rep2::present": "1",
                "replicate::Rep2::primary_credible_call": "1",
                "replicate::RepUnavailable::present": "0",
                "replicate::RepUnavailable::primary_credible_call": "0",
            },
            {
                "candidate_key": "g2|intORF_altframe", "n_primary_credible": "1",
                "gorf_id": "g|chr2|-|blocks|h222222abcdef", "overlap_type": "F2",
                "gene_id": "ENSG000002.1", "gene_name": "GENE2",
                "max_core_reads": "12",
                "replicate::Rep1::present": "1",
                "replicate::Rep1::primary_credible_call": "0",
                "replicate::Rep2::present": "1",
                "replicate::Rep2::primary_credible_call": "1",
                "replicate::RepUnavailable::present": "0",
                "replicate::RepUnavailable::primary_credible_call": "0",
            },
        ]

    def test_call_counts_keep_observed_significant_and_primary_distinct(self) -> None:
        counts = MODULE.call_counts(self.long_rows, self.replicates, 0.05)
        self.assertEqual(counts["Rep1"], {"total": 2, "significant": 1, "primary": 1})
        self.assertEqual(counts["Rep2"], {"total": 2, "significant": 2, "primary": 2})
        self.assertEqual(counts["RepUnavailable"], {"total": 0, "significant": 0, "primary": 0})

    def test_gene_context_filter_removes_candidates_from_long_and_consensus(self) -> None:
        context_rows = [dict(row) for row in self.consensus_rows]
        context_rows[0]["gene_level_pure_intorf_eligible"] = "1"
        context_rows[1]["gene_level_pure_intorf_eligible"] = "0"
        filtered_long, filtered_consensus, excluded = MODULE.filter_gene_level_pure_intorfs(
            self.long_rows, context_rows
        )
        self.assertEqual([row["candidate_key"] for row in filtered_consensus], ["g1|intORF_altframe"])
        self.assertEqual(len(filtered_long), 2)
        self.assertEqual(excluded, 1)
        counts = MODULE.call_counts(filtered_long, self.replicates, 0.05)
        self.assertEqual(counts["Rep1"]["primary"], 1)
        self.assertEqual(counts["Rep2"]["primary"], 1)

    def test_support_combinations_use_primary_calls_only(self) -> None:
        counts = MODULE.support_combinations(self.consensus_rows, self.replicates)
        self.assertEqual(counts[("Rep1", "Rep2")], 1)
        self.assertEqual(counts[("Rep2",)], 1)

    def test_pairwise_lambda_uses_only_shared_candidates(self) -> None:
        pairs = MODULE.paired_lambda_values(self.long_rows, self.replicates)
        self.assertEqual(len(pairs), 1)
        pair = pairs[0]
        self.assertEqual((pair.left, pair.right), ("Rep1", "Rep2"))
        np.testing.assert_allclose(pair.x, [0.8, 0.2])
        np.testing.assert_allclose(pair.y, [0.7, 0.6])
        np.testing.assert_array_equal(pair.left_primary, [True, False])
        np.testing.assert_array_equal(pair.right_primary, [True, True])

    def test_lambda_scatter_matrices_are_dynamic_for_other_replicate_counts(self) -> None:
        for count in (1, 2, 3, 5):
            with self.subTest(replicate_count=count):
                replicates = [f"R{index + 1}" for index in range(count)]
                rows = []
                for replicate_index, replicate in enumerate(replicates):
                    for candidate in range(6):
                        rows.append({
                            "candidate_key": f"candidate-{candidate}",
                            "replicate_id": replicate,
                            "lambda_hat": str((candidate + replicate_index + 1) / (count + 7)),
                            "primary_credible_call": "1" if candidate < 3 else "0",
                        })
                pairs = MODULE.paired_lambda_values(rows, replicates)
                figures = (
                    MODULE.plot_lambda_concordance(pairs, rows, replicates),
                    MODULE.plot_all_eligible_lambda_concordance(pairs, rows, replicates),
                )
                try:
                    for figure in figures:
                        self.assertEqual(len(figure.axes), count * count)
                        for row in range(count):
                            for column in range(count):
                                axis = figure.axes[row * count + column]
                                self.assertEqual(axis.axison, column <= row)
                finally:
                    for figure in figures:
                        MODULE.plt.close(figure)

    def test_lambda_heatmap_matrix_uses_primary_union_and_absolute_lambda(self) -> None:
        rows, matrix = MODULE.lambda_candidate_rows(
            self.consensus_rows, self.long_rows, self.replicates
        )
        self.assertEqual(
            [row["candidate_key"] for row in rows],
            ["g1|intORF_altframe", "g2|intORF_altframe"],
        )
        np.testing.assert_allclose(matrix[:, :2], [[0.8, 0.7], [0.2, 0.6]])
        self.assertTrue(np.all(np.isnan(matrix[:, 2])))

    def test_pfpkm_alignment_and_correlations_use_log2p1_primary_union(self) -> None:
        abundance = {
            "Rep1": {
                "g1|intORF_altframe": {
                    "abundance_status": "quantified", "intorf_pFPKM": "3",
                },
                "g2|intORF_altframe": {
                    "abundance_status": "quantified", "intorf_pFPKM": "1",
                },
            },
            "Rep2": {
                "g1|intORF_altframe": {
                    "abundance_status": "quantified", "intorf_pFPKM": "7",
                },
                "g2|intORF_altframe": {
                    "abundance_status": "quantified", "intorf_pFPKM": "3",
                },
            },
        }
        rows, matrix = MODULE.pfpkm_candidate_rows(
            self.consensus_rows, self.replicates, abundance
        )
        self.assertEqual(
            [row["candidate_key"] for row in rows],
            ["g1|intORF_altframe", "g2|intORF_altframe"],
        )
        np.testing.assert_allclose(matrix[:, :2], [[3.0, 7.0], [1.0, 3.0]])
        self.assertTrue(np.all(np.isnan(matrix[:, 2])))
        pairs = MODULE.paired_pfpkm_values(rows, matrix, self.replicates)
        pair = next(item for item in pairs if (item.left, item.right) == ("Rep1", "Rep2"))
        np.testing.assert_allclose(pair.x, [2.0, 1.0])
        np.testing.assert_allclose(pair.y, [3.0, 2.0])
        self.assertAlmostEqual(pair.pearson, 1.0)
        self.assertAlmostEqual(pair.spearman, 1.0)

        figure = MODULE.plot_pfpkm_correlation(pairs, matrix, self.replicates)
        try:
            self.assertEqual(len(figure.axes), 9)
            self.assertFalse(figure.axes[1].axison)
            self.assertFalse(figure.axes[2].axison)
            self.assertFalse(figure.axes[5].axison)
            self.assertGreater(len(figure.axes[0].patches), 0)
        finally:
            MODULE.plt.close(figure)

    def test_pfpkm_heatmap_contains_absolute_abundance_only(self) -> None:
        matrix = np.asarray([[3.0, 7.0, np.nan], [1.0, 3.0, np.nan]])
        figure = MODULE.plot_pfpkm_heatmap(
            self.consensus_rows, matrix, self.replicates, 2
        )
        try:
            self.assertEqual(len(figure.axes), 2)  # heatmap plus one colorbar
            self.assertEqual(figure.axes[0].get_title(), "Absolute abundance")
            all_text = " ".join(
                text.get_text() for axis in figure.axes for text in axis.texts
            )
            self.assertNotIn("relative abundance", all_text.lower())
            self.assertNotIn("z-score", all_text.lower())
        finally:
            MODULE.plt.close(figure)

    def test_pfpkm_scatter_matrix_is_dynamic_for_other_replicate_counts(self) -> None:
        for count in (1, 2, 3, 5):
            with self.subTest(replicate_count=count):
                replicates = [f"R{index + 1}" for index in range(count)]
                matrix = np.asarray([
                    [float(candidate + replicate + 1) for replicate in range(count)]
                    for candidate in range(6)
                ])
                pairs = MODULE.paired_pfpkm_values([], matrix, replicates)
                figure = MODULE.plot_pfpkm_correlation(pairs, matrix, replicates)
                try:
                    self.assertEqual(len(figure.axes), count * count)
                    for row in range(count):
                        for column in range(count):
                            axis = figure.axes[row * count + column]
                            self.assertEqual(axis.axison, column <= row)
                finally:
                    MODULE.plt.close(figure)

    def test_primary_union_excludes_candidates_not_primary_in_either_replicate(self) -> None:
        rows = [*self.long_rows,
            {
                "candidate_key": "g3|intORF_altframe", "replicate_id": "Rep1",
                "lambda_hat": "0.4", "primary_credible_call": "0",
            },
            {
                "candidate_key": "g3|intORF_altframe", "replicate_id": "Rep2",
                "lambda_hat": "0.5", "primary_credible_call": "0",
            },
        ]
        pair = MODULE.paired_lambda_values(rows, self.replicates)[0]
        union = pair.left_primary | pair.right_primary
        np.testing.assert_allclose(pair.x[union], [0.8, 0.2])
        np.testing.assert_allclose(pair.y[union], [0.7, 0.6])
        self.assertEqual(int(np.sum(pair.left_primary & pair.right_primary)), 1)

    def test_heatmap_marks_unavailable_replicate_separately(self) -> None:
        matrix, labels, _shown, unavailable = MODULE.reproducibility_matrix(
            self.consensus_rows, self.replicates, 10
        )
        self.assertEqual(labels, ["GENE1", "GENE2"])
        self.assertEqual(unavailable, {"RepUnavailable"})
        np.testing.assert_array_equal(matrix, [[2, 2, -1], [1, 2, -1]])

    def test_global_heatmap_includes_all_primary_candidates_and_groups_patterns(self) -> None:
        rows = [dict(row) for row in self.consensus_rows]
        rows.extend([
            {
                "candidate_key": "g3|intORF_altframe", "n_primary_credible": "1",
                "max_core_reads": "8", "gene_name": "GENE3",
                "replicate::Rep1::present": "1", "replicate::Rep1::primary_credible_call": "1",
                "replicate::Rep2::present": "1", "replicate::Rep2::primary_credible_call": "0",
                "replicate::RepUnavailable::present": "0",
                "replicate::RepUnavailable::primary_credible_call": "0",
            },
            {
                "candidate_key": "never-primary", "n_primary_credible": "0",
                "max_core_reads": "100", "gene_name": "IGNORED",
                "replicate::Rep1::present": "1", "replicate::Rep1::primary_credible_call": "0",
                "replicate::Rep2::present": "1", "replicate::Rep2::primary_credible_call": "0",
                "replicate::RepUnavailable::present": "0",
                "replicate::RepUnavailable::primary_credible_call": "0",
            },
        ])
        matrix, labels, shown, unavailable = MODULE.reproducibility_matrix(
            rows, self.replicates, 100, group_by_pattern=True
        )
        self.assertEqual(labels, ["GENE1", "GENE3", "GENE2"])
        self.assertEqual([row["candidate_key"] for row in shown], [
            "g1|intORF_altframe", "g3|intORF_altframe", "g2|intORF_altframe",
        ])
        self.assertEqual(unavailable, {"RepUnavailable"})
        np.testing.assert_array_equal(matrix, [[2, 2, -1], [2, 1, -1], [1, 2, -1]])

    def test_global_heatmap_height_is_configurable(self) -> None:
        matrix = np.asarray([[2, 1, -1], [1, 2, -1]])
        figure = MODULE.plot_global_reproducibility_heatmap(
            matrix, self.consensus_rows, self.replicates, 5.5
        )
        try:
            self.assertAlmostEqual(figure.get_size_inches()[1], 5.5)
        finally:
            MODULE.plt.close(figure)

    def test_quality_labels_include_psite_depth_frame0_and_a0(self) -> None:
        metadata = {
            "Rep1": {
                "replicate_id": "Rep1", "psite_alignments": "22973977",
                "frame0_prop": "0.9332772207", "A0": "3.284990995",
                "selected_length_offsets": "28:12,29:12,30:12",
            }
        }
        labels = MODULE.replicate_quality_labels(["Rep1", "Rep2"], metadata)
        self.assertEqual(
            labels[0], "Rep1\n23.0M P-sites · L 28/29/30\nF0 93.3% · A0 3.28"
        )
        self.assertEqual(labels[1], "Rep2")

    def test_batch_summary_metadata_columns_are_accepted(self) -> None:
        import tempfile

        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "batch_summary.tsv"
            path.write_text(
                "replicate_id\tselected_length_offsets\tpsite_alignments\t"
                "qc_frame0_prop\tqc_A0\n"
                "Rep1\t28:12,29:12\t1000\t0.9\t2.5\n",
                encoding="utf-8",
            )
            metadata = MODULE.read_sample_metadata(path)
        self.assertEqual(metadata["Rep1"]["frame0_prop"], "0.9")
        self.assertEqual(metadata["Rep1"]["A0"], "2.5")

    def test_heatmaps_render_quality_labels_when_metadata_is_available(self) -> None:
        matrix = np.asarray([[2, 1, -1], [1, 2, -1]])
        metadata = {
            "Rep1": {
                "replicate_id": "Rep1", "psite_alignments": "22973977",
                "frame0_prop": "0.9332772207", "A0": "3.284990995",
                "selected_length_offsets": "28:12,29:12,30:12",
            }
        }
        figures = [
            MODULE.plot_reproducibility_heatmap(
                matrix, ["GENE1", "GENE2"], self.replicates, 2, metadata
            ),
            MODULE.plot_global_reproducibility_heatmap(
                matrix, self.consensus_rows, self.replicates, 5.5, metadata
            ),
        ]
        try:
            for figure in figures:
                labels = [tick.get_text() for tick in figure.axes[0].get_xticklabels()]
                self.assertIn("23.0M P-sites", labels[0])
                self.assertIn("F0 93.3% · A0 3.28", labels[0])
        finally:
            for figure in figures:
                MODULE.plt.close(figure)

    def test_heatmap_disambiguates_multiple_candidates_from_the_same_gene(self) -> None:
        rows = [dict(row) for row in self.consensus_rows]
        rows[1]["gene_name"] = "GENE1"
        labels = MODULE.candidate_display_labels(rows)
        self.assertEqual(labels, ["GENE1 (F1; h...abcdef)", "GENE1 (F2; h...abcdef)"])

    def test_heatmap_label_falls_back_to_gene_id_then_candidate_key(self) -> None:
        rows = [
            {"candidate_key": "candidate-1", "gene_id": "ENSG000003.4"},
            {"candidate_key": "candidate-2"},
        ]
        self.assertEqual(
            MODULE.candidate_display_labels(rows),
            ["ENSG000003.4", "candidate-2"],
        )


if __name__ == "__main__":
    unittest.main()
