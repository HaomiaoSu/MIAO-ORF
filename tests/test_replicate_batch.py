#!/usr/bin/env python3
"""Tests for replicate-batch parsing and failure-isolated summaries."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "miao_orf.py"
SPEC = importlib.util.spec_from_file_location("miao_orf_batch_under_test", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class ReplicateBatchTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def write_batch(self, rows: list[str]) -> Path:
        path = self.root / "replicates.tsv"
        path.write_text(
            "replicate_id\tbam\tinput_mode\tlength_offsets\t"
            "ribotish_para\tribotish_quality\n"
            + "\n".join(rows)
            + "\n",
            encoding="utf-8",
        )
        return path

    def test_reads_explicit_and_ribotish_replicates(self) -> None:
        path = self.write_batch(
            [
                "Rep1\tdata/rep1.bam\texplicit\t28:12,29:12,30:12\t\t",
                "Rep2\tdata/rep2.bam\tribotish\t\tqc/rep2.para.py\tqc/rep2.qual.txt",
            ]
        )
        rows = MODULE.read_batch_replicates(path)
        self.assertEqual([row.replicate_id for row in rows], ["Rep1", "Rep2"])
        self.assertEqual(rows[0].length_offsets, {28: 12, 29: 12, 30: 12})
        self.assertEqual(rows[0].bam, self.root / "data" / "rep1.bam")
        self.assertEqual(rows[1].ribotish_para, self.root / "qc" / "rep2.para.py")
        self.assertIsNone(rows[1].length_offsets)

    def test_rejects_duplicate_ids_and_mixed_row_inputs(self) -> None:
        duplicate = self.write_batch(
            [
                "Rep1\ta.bam\texplicit\t28:12\t\t",
                "rep1\tb.bam\texplicit\t29:12\t\t",
            ]
        )
        with self.assertRaisesRegex(ValueError, "duplicate replicate_id"):
            MODULE.read_batch_replicates(duplicate)

        mixed = self.write_batch(
            ["Rep2\tb.bam\tribotish\t28:12\tb.para.py\tb.qual.txt"]
        )
        with self.assertRaisesRegex(ValueError, "cannot include length_offsets"):
            MODULE.read_batch_replicates(mixed)

    def test_failed_setup_still_has_a_batch_summary_row(self) -> None:
        outcome = MODULE.PipelineOutcome(
            sample="Rep1",
            input_mode="explicit",
            status="failed",
            error="bad input",
            failed_stage="setup",
            manifest=None,
            mode_paths={},
        )
        rows = MODULE.build_batch_summary_rows([outcome])
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["replicate_id"], "Rep1")
        self.assertEqual(rows[0]["status"], "failed")

        batch = self.write_batch(
            ["Rep1\tmissing.bam\texplicit\t28:12\t\t"]
        )
        args = SimpleNamespace(out_root=str(self.root / "results"), batch=str(batch))
        summary, manifest = MODULE.write_batch_records(
            args, [outcome], [], "2026-01-01T00:00:00+00:00"
        )
        self.assertTrue(summary.is_file())
        self.assertTrue(manifest.is_file())
        self.assertIn("Rep1\texplicit\t\tfailed", summary.read_text(encoding="utf-8"))

    def test_top_level_calls_independent_integration_module(self) -> None:
        out_root = self.root / "results"
        dm_results = self.root / "Rep1.intorf_dm_results.tsv"
        dm_results.write_text("header\n", encoding="utf-8")
        paths = MODULE.Paths(
            project=ROOT,
            source=ROOT / "src" / "miao_orf",
            out_root=out_root,
            reference_name="ref",
            sample="Rep1",
            orf_prefix=out_root / "01_orfscan" / "ref",
            per_chrom_dir=out_root / "01_orfscan" / "per_chrom",
            psite_prefix=out_root / "02_psite" / "Rep1" / "Rep1",
            qc_prefix=out_root / "03_metagene_qc" / "Rep1" / "Rep1",
            dm_prefix=out_root / "04_intorf_dm" / "Rep1" / "accurate" / "Rep1",
            abundance_prefix=out_root / "04_intorf_dm" / "Rep1" / "accurate" / "Rep1",
            context_prefix=out_root / "04_intorf_dm" / "Rep1" / "accurate" / "Rep1",
            codon_prefix=out_root / "04_intorf_dm" / "Rep1" / "accurate" / "Rep1",
            visualization_prefix=out_root / "05_visualization" / "Rep1" / "accurate" / "Rep1",
            torf=out_root / "01_orfscan" / "ref.torf.tsv",
            psite_bam=out_root / "02_psite" / "Rep1" / "Rep1.psite.bam",
            dm_background=out_root / "03_metagene_qc" / "Rep1" / "Rep1.dm_background.tsv",
            dm_results=dm_results,
            dm_run_manifest=out_root / "04_intorf_dm" / "Rep1" / "accurate" / "Rep1.run_manifest.json",
            gene_context_results=out_root / "04_intorf_dm" / "Rep1" / "accurate" / "Rep1.gene_cds_context.tsv",
        )
        dm_run_manifest = paths.dm_run_manifest
        dm_run_manifest.parent.mkdir(parents=True, exist_ok=True)
        dm_run_manifest.write_text("{}\n", encoding="utf-8")
        paths.torf.parent.mkdir(parents=True, exist_ok=True)
        paths.torf.write_text("torf_id\tgorf_id\n", encoding="utf-8")
        abundance_path = Path(f"{paths.abundance_prefix}.intorf_abundance.tsv")
        abundance_path.write_text(
            "gorf_id\toverlap_type\tintorf_pFPKM\tabundance_status\n"
            "g1\tintORF_altframe\t1.0\tquantified\n",
            encoding="utf-8",
        )
        outcome = MODULE.PipelineOutcome(
            "Rep1", "explicit", "completed", "", "", None, {"accurate": paths}
        )
        args = SimpleNamespace(
            out_root=str(out_root),
            batch=str(self.root / "replicates.tsv"),
            python=sys.executable,
            dm_mode="accurate",
            integration_min_replicates=2,
            integration_min_fraction=0.5,
            fdr_threshold=0.05,
            formats="png,pdf",
            dpi=220,
        )

        def fake_stream(command: list[str], _log: Path) -> int:
            prefix = Path(command[command.index("--out-prefix") + 1])
            script_name = Path(command[1]).name
            if script_name == "integrate_replicates.py":
                suffixes = (
                    ".replicate_long.tsv", ".consensus.tsv", ".summary.tsv", ".manifest.json"
                )
            elif script_name == "annotate_gene_cds_context.py":
                suffixes = (
                    ".gene_cds_context.tsv", ".gene_cds_context_only.tsv",
                    ".gene_cds_context_summary.tsv", ".gene_cds_context_manifest.json",
                )
            else:
                figure_names = [
                    "call_counts", "support_combinations", "lambda_concordance",
                    "lambda_heatmap", "primary_reproducibility",
                ]
                if "--abundance" in command:
                    figure_names.extend([
                        "pfpkm_correlation", "pfpkm_correlation_matrix", "pfpkm_heatmap",
                    ])
                suffixes = tuple(
                    f".{figure}.{extension}"
                    for figure in figure_names
                    for extension in ("png", "pdf")
                ) + (
                    ".lambda_matrix.tsv", ".pfpkm_matrix.tsv",
                    ".pfpkm_correlations.tsv", ".plot_summary.tsv", ".plot_manifest.json",
                )
            for suffix in suffixes:
                output = Path(f"{prefix}{suffix}")
                output.parent.mkdir(parents=True, exist_ok=True)
                if suffix.endswith(".json"):
                    output.write_text("{}\n", encoding="utf-8")
                elif suffix.endswith(".png"):
                    output.write_bytes(b"\x89PNG\r\n\x1a\n")
                elif suffix.endswith(".pdf"):
                    output.write_bytes(b"%PDF-1.4\n")
                elif suffix == ".summary.tsv":
                    output.write_text("metric\tvalue\nstatus\tok\n", encoding="utf-8")
                elif suffix in (".gene_cds_context.tsv", ".gene_cds_context_only.tsv"):
                    output.write_text(
                        "gorf_id\tgene_id\tgene_level_orf_class\t"
                        "gene_level_pure_intorf_eligible\tgene_cds_nterm_match\t"
                        "gene_cds_nterm_coordinate_prefix_complete_codons\n",
                        encoding="utf-8",
                    )
                elif suffix == ".gene_cds_context_summary.tsv":
                    output.write_text("metric\tvalue\nstatus\tok\n", encoding="utf-8")
                elif suffix == ".lambda_matrix.tsv":
                    output.write_text(
                        "candidate_key\tgorf_id\toverlap_type\tgene_name\t"
                        "n_primary_credible\tmax_core_reads\n",
                        encoding="utf-8",
                    )
                elif suffix == ".pfpkm_matrix.tsv":
                    output.write_text(
                        "candidate_key\tgorf_id\toverlap_type\tgene_name\t"
                        "n_primary_credible\n",
                        encoding="utf-8",
                    )
                elif suffix == ".pfpkm_correlations.tsv":
                    output.write_text(
                        "replicate_left\treplicate_right\tshared_candidates\t"
                        "pearson_log2p1_pFPKM\tspearman_log2p1_pFPKM\n",
                        encoding="utf-8",
                    )
                else:
                    output.write_text("header\n", encoding="utf-8")
            return 0

        with mock.patch.object(MODULE, "stream_command", side_effect=fake_stream) as called:
            records = MODULE.run_replicate_integrations(args, [outcome])
        self.assertEqual(records[0].status, "completed")
        self.assertEqual(records[0].context_status, "completed")
        self.assertEqual(records[0].visualization_status, "completed")
        self.assertEqual(called.call_count, 3)
        integration_command = called.call_args_list[0].args[0]
        context_command = called.call_args_list[1].args[0]
        visualization_command = called.call_args_list[2].args[0]
        self.assertEqual("integrate_replicates.py", Path(integration_command[1]).name)
        self.assertIn(f"Rep1={dm_results}", integration_command)
        self.assertIn(f"Rep1={dm_run_manifest}", integration_command)
        self.assertEqual("annotate_gene_cds_context.py", Path(context_command[1]).name)
        self.assertIn(str(paths.torf), context_command)
        self.assertIn(str(Path(f"{out_root / '06_replicate_integration' / 'accurate' / 'replicates'}.consensus.tsv")), context_command)
        self.assertEqual("visualize_replicate_integration.py", Path(visualization_command[1]).name)
        self.assertIn("--gene-level-pure-intorf-only", visualization_command)
        self.assertIn(
            str(Path(f"{out_root / '06_replicate_integration' / 'accurate' / 'replicates'}.consensus.gene_cds_context.tsv")),
            visualization_command,
        )
        self.assertIn("--sample-metadata", visualization_command)
        self.assertIn("--abundance", visualization_command)
        self.assertIn(f"Rep1={abundance_path}", visualization_command)
        metadata_path = Path(
            visualization_command[visualization_command.index("--sample-metadata") + 1]
        )
        self.assertTrue(metadata_path.is_file())
        self.assertIn(
            "replicate_id\tselected_length_offsets\tpsite_alignments\tframe0_prop\tA0",
            metadata_path.read_text(),
        )


if __name__ == "__main__":
    unittest.main()
