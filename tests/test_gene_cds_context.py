#!/usr/bin/env python3
"""Tests for gene-level annotated-CDS N-terminal reuse postprocessing."""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src" / "miao_orf" / "annotate_gene_cds_context.py"
SPEC = importlib.util.spec_from_file_location("gene_cds_context_under_test", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


TORF_FIELDS = [
    "torf_id", "gorf_id", "gene_id", "gene_name", "transcript_id",
    "transcript_name", "tx_tags", "chrom", "strand", "is_annotated_cds",
    "genomic_block_starts1", "genomic_block_ends1", "peptide_len", "peptide",
]


class GeneCdsContextTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def write_tsv(self, path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def atp1b3_rows(self) -> list[dict[str, object]]:
        return [
            {
                "torf_id": "t|ENST00000462082.1|118|235|f1",
                "gorf_id": "g|chr3|+|blocks|hb9ceccfd347d11ef",
                "gene_id": "ENSG00000069849.11",
                "gene_name": "ATP1B3",
                "transcript_id": "ENST00000462082.1",
                "transcript_name": "ATP1B3-202",
                "tx_tags": "basic",
                "chrom": "chr3",
                "strand": "+",
                "is_annotated_cds": "0",
                "genomic_block_starts1": "141903682,141925531",
                "genomic_block_ends1": "141903748,141925580",
                "peptide_len": "38",
                "peptide": "MWVMLQTLNDEVPKYRDQIPSPGWVSTAIGCCSGQLCS",
            },
            {
                "torf_id": "t|ENST00000286371.8|CDS",
                "gorf_id": "g|ATP1B3|canonical_CDS",
                "gene_id": "ENSG00000069849.11",
                "gene_name": "ATP1B3",
                "transcript_id": "ENST00000286371.8",
                "transcript_name": "ATP1B3-201",
                "tx_tags": "MANE_Select,CCDS,Ensembl_canonical",
                "chrom": "chr3",
                "strand": "+",
                "is_annotated_cds": "1",
                "genomic_block_starts1": "141876802,141903620,141907167,141913652,141915970,141921977,141925531",
                "genomic_block_ends1": "141876910,141903748,141907274,141913836,141916020,141922063,141925701",
                "peptide_len": "279",
                "peptide": "A" * 57 + "MWVMLQTLNDEVPKYRDQIPSPG" + "L" * 199,
            },
        ]

    def test_atp1b3_like_splice_divergence_is_gene_level_nterm_reuse(self) -> None:
        torf = self.root / "ref.torf.tsv"
        result = self.root / "sample.intorf_dm_results.tsv"
        prefix = self.root / "sample"
        self.write_tsv(torf, self.atp1b3_rows(), TORF_FIELDS)
        self.write_tsv(
            result,
            [{
                "gorf_id": "g|chr3|+|blocks|hb9ceccfd347d11ef",
                "overlap_type": "intORF_altframe",
                "gene_id": "ENSG00000069849.11",
                "gene_name": "ATP1B3",
                "primary_credible_call": "1",
                "classification": "credible_extra_ORF_like_signal",
            }],
            ["gorf_id", "overlap_type", "gene_id", "gene_name", "primary_credible_call", "classification"],
        )
        outputs = MODULE.run(SimpleNamespace(
            torf=str(torf), input_tsv=str(result), out_prefix=str(prefix), min_prefix_codons=5
        ))
        with outputs["annotated"].open("r", encoding="utf-8", newline="") as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(row["gene_level_orf_class"], "annotated_cds_nterm_reuse_with_splice_derived_cterm")
        self.assertEqual(row["gene_level_pure_intorf_eligible"], "0")
        self.assertEqual(row["gene_cds_nterm_best_transcript_id"], "ENST00000286371.8")
        self.assertEqual(row["gene_cds_nterm_coordinate_prefix_nt"], "67")
        self.assertEqual(row["gene_cds_nterm_coordinate_prefix_complete_codons"], "22")
        self.assertEqual(row["gene_cds_nterm_peptide_prefix_aa"], "23")
        self.assertEqual(row["gene_cds_nterm_break_reason"], "splice_or_coordinate_path_divergence")
        self.assertEqual(row["classification"], "credible_extra_ORF_like_signal")
        manifest = json.loads(outputs["manifest"].read_text(encoding="utf-8"))
        self.assertFalse(manifest["algorithm"]["transcript_level_classification_modified"])
        self.assertFalse(manifest["algorithm"]["interval_overlap_alone_is_match"])

    def test_codon_phase_mismatch_is_not_a_match(self) -> None:
        candidate = MODULE.to_torf_record(self.atp1b3_rows()[0], 2)
        cds_row = self.atp1b3_rows()[1]
        cds_row["genomic_block_starts1"] = "141876801,141903620,141907167,141913652,141915970,141921977,141925531"
        cds_row["genomic_block_ends1"] = "141876910,141903748,141907274,141913836,141916020,141922063,141925701"
        cds_row["peptide_len"] = "279"
        cds_row["peptide"] = "A" * 279
        cds = MODULE.to_torf_record(cds_row, 3)
        self.assertIsNone(MODULE.compare(candidate, cds))

    def test_minus_strand_translation_path_and_full_match(self) -> None:
        candidate_row = {
            "torf_id": "candidate", "gorf_id": "g1", "gene_id": "gene1", "gene_name": "G",
            "transcript_id": "tx-candidate", "transcript_name": "", "tx_tags": "",
            "chrom": "chr1", "strand": "-", "is_annotated_cds": "0",
            "genomic_block_starts1": "100,200", "genomic_block_ends1": "105,205",
            "peptide_len": "4", "peptide": "MPEP",
        }
        cds_row = {
            **candidate_row,
            "torf_id": "cds", "gorf_id": "cds-g", "transcript_id": "tx-cds",
            "is_annotated_cds": "1",
        }
        candidate = MODULE.to_torf_record(candidate_row, 2)
        cds = MODULE.to_torf_record(cds_row, 3)
        self.assertEqual(candidate.coding_path, (205, 204, 203, 202, 201, 200, 105, 104, 103, 102, 101, 100))
        match = MODULE.compare(candidate, cds)
        self.assertIsNotNone(match)
        assert match is not None
        self.assertEqual(match.coordinate_prefix_nt, 12)
        context = MODULE.context_for_candidate([candidate], [cds], 2)
        self.assertEqual(context["gene_level_orf_class"], "annotated_cds_derived_full_coordinate_path")


if __name__ == "__main__":
    unittest.main()
