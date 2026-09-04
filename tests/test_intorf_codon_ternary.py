from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src" / "miao_orf" / "export_intorf_codon_ternary.py"
SPEC = importlib.util.spec_from_file_location("test_intorf_codon_ternary_module", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class ExactTernaryRatioTests(unittest.TestCase):
    def test_reduces_proportional_counts_without_rounding(self) -> None:
        self.assertEqual(MODULE.reduced_phase_ratio((1, 1, 0)), (1, 1, 0))
        self.assertEqual(MODULE.reduced_phase_ratio((2, 2, 0)), (1, 1, 0))
        self.assertEqual(MODULE.reduced_phase_ratio((1, 2, 3)), (1, 2, 3))


class OrfSelectionTests(unittest.TestCase):
    def write_results(self, root: Path) -> Path:
        path = root / "results.tsv"
        columns = [
            "torf_id", "gorf_id", "overlap_type", "primary_credible_call",
            "gene_level_pure_intorf_eligible",
        ]
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            writer.writerows([
                {
                    "torf_id": "t1", "gorf_id": "g1", "overlap_type": "F1",
                    "primary_credible_call": "1", "gene_level_pure_intorf_eligible": "1",
                },
                {
                    "torf_id": "t2", "gorf_id": "g2", "overlap_type": "F2",
                    "primary_credible_call": "0", "gene_level_pure_intorf_eligible": "0",
                },
            ])
        return path

    def test_default_selects_credible_rows(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            rows, _ = MODULE.load_selected_rows(
                str(self.write_results(Path(temporary))), False
            )
            self.assertEqual([row["torf_id"] for row in rows], ["t1"])

    def test_explicit_single_or_batch_ids_override_credibility(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = self.write_results(Path(temporary))
            rows, _ = MODULE.load_selected_rows(str(path), True, requested_ids={"t2"})
            self.assertEqual([row["gorf_id"] for row in rows], ["g2"])
            rows, _ = MODULE.load_selected_rows(
                str(path), False, requested_ids={"g1|F1", "g2"}
            )
            self.assertEqual({row["torf_id"] for row in rows}, {"t1", "t2"})

    def test_text_and_tsv_orf_lists(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            plain = root / "plain.txt"
            plain.write_text("t1\ng2\n", encoding="utf-8")
            self.assertEqual(MODULE.load_orf_list(str(plain)), {"t1", "g2"})
            table = root / "table.tsv"
            table.write_text("gorf_id\ttorf_id\ng1\t\n\tt2\n", encoding="utf-8")
            self.assertEqual(MODULE.load_orf_list(str(table)), {"g1", "t2"})

    def test_missing_requested_id_is_an_error(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = self.write_results(Path(temporary))
            with self.assertRaisesRegex(SystemExit, "requested ORF ID"):
                MODULE.load_selected_rows(str(path), False, requested_ids={"missing"})


@unittest.skipUnless(importlib.util.find_spec("pysam"), "pysam is not installed")
class IntorfCodonTernaryTests(unittest.TestCase):
    def test_exports_and_reconciles_credible_candidate(self) -> None:
        import pysam

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bam_path = root / "sample.psite.bam"
            header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
            positions = [103, 103, 103, 106, 107, 111, 111]
            with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
                for index, position in enumerate(positions):
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

            torf_path = root / "reference.torf.tsv"
            torf_columns = [
                "torf_id", "gorf_id", "gene_id", "gene_name", "transcript_id",
                "chrom", "strand", "blockSizes", "genomic_block_starts1",
                "t_start", "t_end", "aa_len", "peptide_len", "cds_first_t",
                "orf_biotype", "has_mane",
            ]
            with torf_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=torf_columns, delimiter="\t")
                writer.writeheader()
                writer.writerow({
                    "torf_id": "t1", "gorf_id": "g1", "gene_id": "GENE1",
                    "gene_name": "Gene1", "transcript_id": "tx1", "chrom": "chr1",
                    "strand": "+", "blockSizes": "15", "genomic_block_starts1": "101",
                    "t_start": "1", "t_end": "16", "aa_len": "4", "peptide_len": "4",
                    "cds_first_t": "0", "orf_biotype": "intORF_altframe", "has_mane": "1",
                })

            results_path = root / "sample.gene_cds_context.tsv"
            result_columns = [
                "torf_id", "gorf_id", "gene_id", "gene_name", "transcript_id",
                "overlap_type", "classification", "q_BH", "lambda_hat",
                "primary_credible_call", "gene_level_pure_intorf_eligible",
                "n_core_codons", "n_active_core_codons", "core_reads",
                "analyzed_core_reads", "observed_phase0_core_reads",
                "observed_phase1_core_reads", "observed_phase2_core_reads",
            ]
            with results_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=result_columns, delimiter="\t")
                writer.writeheader()
                writer.writerow({
                    "torf_id": "t1", "gorf_id": "g1", "gene_id": "GENE1",
                    "gene_name": "Gene1", "transcript_id": "tx1", "overlap_type": "F1",
                    "classification": "credible_translated_intORF", "q_BH": "0.01",
                    "lambda_hat": "0.4", "primary_credible_call": "1",
                    "gene_level_pure_intorf_eligible": "1", "n_core_codons": "3",
                    "n_active_core_codons": "3", "core_reads": "7",
                    "analyzed_core_reads": "7", "observed_phase0_core_reads": "4",
                    "observed_phase1_core_reads": "1", "observed_phase2_core_reads": "2",
                })

            manifest_path = root / "sample.run_manifest.json"
            manifest_path.write_text(json.dumps({
                "parameters": {
                    "exclude_start_codons": 1,
                    "exclude_stop_codons": 1,
                    "rl": None,
                    "rl_source": "auto",
                }
            }), encoding="utf-8")
            prefix = root / "out" / "sample"
            self.assertEqual(MODULE.main([
                "--input-results", str(results_path),
                "--torf", str(torf_path),
                "--psite-bam", str(bam_path),
                "--dm-run-manifest", str(manifest_path),
                "--out-prefix", str(prefix),
                "--sample", "sample",
                "--gene-level-pure-intorf-only",
            ]), 0)

            with gzip.open(
                f"{prefix}.orf_psite_codons.tsv.gz", "rt", encoding="utf-8"
            ) as handle:
                codons = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(codons), 3)
            self.assertEqual(
                [(row["phase0_reads"], row["phase1_reads"], row["phase2_reads"]) for row in codons],
                [("3", "0", "0"), ("1", "1", "0"), ("0", "0", "2")],
            )
            with Path(f"{prefix}.orf_psite_ternary.tsv").open(
                "r", encoding="utf-8", newline=""
            ) as handle:
                ternary = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(ternary), 3)
            self.assertTrue(all(row["candidate_percentage"].startswith("33.333") for row in ternary))
            exact_points = {
                (
                    row["phase_ratio_P0"], row["phase_ratio_P1"], row["phase_ratio_P2"],
                    round(float(row["exact_P0_percent"]), 6),
                    round(float(row["exact_P1_percent"]), 6),
                    round(float(row["exact_P2_percent"]), 6),
                )
                for row in ternary
            }
            self.assertEqual(exact_points, {
                ("1", "0", "0", 100.0, 0.0, 0.0),
                ("1", "1", "0", 50.0, 50.0, 0.0),
                ("0", "0", "1", 0.0, 0.0, 100.0),
            })
            with Path(f"{prefix}.orf_psite_summary.tsv").open(
                "r", encoding="utf-8", newline=""
            ) as handle:
                summary = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(summary["count_reconciliation_status"], "exact_match")
            self.assertAlmostEqual(float(summary["codon_equal_mean_P0_percent"]), 50.0)
            self.assertAlmostEqual(float(summary["codon_equal_mean_P1_percent"]), 100.0 / 6.0)
            self.assertAlmostEqual(float(summary["codon_equal_mean_P2_percent"]), 100.0 / 3.0)
            manifest = json.loads(Path(
                f"{prefix}.orf_psite_manifest.json"
            ).read_text(encoding="utf-8"))
            self.assertNotIn("grid_percent", manifest["parameters"])
            self.assertIn("no binning", manifest["semantics"]["exact_point_aggregation"])


if __name__ == "__main__":
    unittest.main()
