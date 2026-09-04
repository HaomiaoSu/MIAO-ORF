#!/usr/bin/env python3
"""Tests for evidence-preserving replicate integration."""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src" / "miao_orf" / "integrate_replicates.py"
SPEC = importlib.util.spec_from_file_location("replicate_integration_under_test", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class ReplicateIntegrationTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def write_result(self, name: str, rows: list[tuple[str, str, str, int, float, float, int]]) -> Path:
        path = self.root / f"{name}.tsv"
        header = [
            "gorf_id", "overlap_type", "gene_id", "gene_name", "chrom", "strand",
            "orf_biotype", "aa_len", "peptide_len", "classification",
            "primary_credible_call", "q_BH", "q_BY", "p_final", "lambda_hat",
            "core_reads", "n_active_core_codons", "active_core_codon_frac",
            "mixture_geometry_consistent", "distance_to_mixture_segment",
        ]
        lines = ["\t".join(header)]
        for gorf, overlap, gene, primary, q_value, lam, reads in rows:
            classification = (
                MODULE.PRIMARY_CREDIBLE_CLASS if primary else "not_significant"
            )
            values = [
                gorf, overlap, gene, gene, "chr1", "+", "intORF_altframe", "20", "20",
                classification, str(primary), str(q_value), str(q_value), str(q_value),
                str(lam), str(reads), "8", "0.5", "1", "0.01",
            ]
            lines.append("\t".join(values))
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return path

    def write_manifest(
        self, replicate: str, result: Path, parameter_override: tuple[str, object] | None = None
    ) -> Path:
        parameters = {key: f"same::{key}" for key in MODULE.SEMANTIC_PARAMETER_KEYS}
        if parameter_override is not None:
            parameters[parameter_override[0]] = parameter_override[1]
        manifest = {
            "program": "miao-orf-intorf-dm",
            "version": "1.0.0",
            "output_schema_version": "1.0",
            "statistical_engine_id": "dm-formal-adaptive-importance-v1",
            "benchmark_certification_id": "benchmark-520-dual-20260813",
            "parameters": parameters,
            "source": {"caller": {"sha256": "caller-sha256"}},
            "inputs": {
                "torf_files": [
                    {
                        "path": "/reference/shared.torf.tsv", "exists": True,
                        "size_bytes": 1234, "mtime_utc": "2026-01-01T00:00:00+00:00",
                        "sha256": "torf-sha256",
                    }
                ]
            },
            "output_schema": result.read_text(encoding="utf-8").splitlines()[0].split("\t"),
            "outputs": [
                {
                    "path": str(result.absolute()),
                    "size_bytes": result.stat().st_size,
                    "sha256": MODULE.sha256_file(result),
                }
            ],
            "result_rows": 1,
        }
        path = self.root / f"{replicate}.run_manifest.json"
        path.write_text(json.dumps(manifest), encoding="utf-8")
        return path

    def test_reproducible_consensus_uses_support_not_combined_pvalues(self) -> None:
        rep1 = self.write_result(
            "rep1",
            [("g1", "intORF_altframe", "GENE1", 1, 0.01, 0.8, 20),
             ("g2", "intORF_altframe", "GENE2", 0, 0.8, 0.1, 5)],
        )
        rep2 = self.write_result(
            "rep2",
            [("g1", "intORF_altframe", "GENE1", 1, 0.02, 0.7, 18),
             ("g2", "intORF_altframe", "GENE2", 1, 0.03, 0.6, 12)],
        )
        rep3 = self.write_result(
            "rep3",
            [("g1", "intORF_altframe", "GENE1", 0, 0.4, 0.2, 8)],
        )
        long_rows, consensus, summary = MODULE.integrate(
            ["Rep1", "Rep2", "Rep3", "RepUnavailable"],
            {"Rep1": rep1, "Rep2": rep2, "Rep3": rep3},
            min_replicates=2,
            min_fraction=0.5,
            fdr_threshold=0.05,
        )
        by_key = {row["candidate_key"]: row for row in consensus}
        self.assertEqual(
            by_key["g1|intORF_altframe"]["consensus_classification"],
            "reproducible_primary_credible",
        )
        self.assertEqual(by_key["g1|intORF_altframe"]["n_primary_credible"], 2)
        self.assertEqual(by_key["g1|intORF_altframe"]["n_bh_significant"], 2)
        self.assertEqual(
            by_key["g2|intORF_altframe"]["consensus_classification"],
            "primary_credible_not_reproducible",
        )
        self.assertEqual(len(long_rows), 5)
        self.assertEqual(summary["n_replicates_available"], 3)
        self.assertEqual(summary["unavailable_replicates"], "RepUnavailable")

    def test_single_available_replicate_is_not_called_reproducible(self) -> None:
        rep1 = self.write_result(
            "rep1", [("g1", "intORF_altframe", "GENE1", 1, 0.01, 0.8, 20)]
        )
        _long, consensus, _summary = MODULE.integrate(
            ["Rep1", "Rep2"], {"Rep1": rep1}, 2, 0.5, 0.05
        )
        self.assertEqual(
            consensus[0]["consensus_classification"],
            "single_replicate_primary_credible",
        )

    def test_run_manifests_must_share_the_formal_contract(self) -> None:
        rep1 = self.write_result(
            "rep1_compatible", [("g1", "intORF_altframe", "GENE1", 1, 0.01, 0.8, 20)]
        )
        rep2 = self.write_result(
            "rep2_compatible", [("g1", "intORF_altframe", "GENE1", 1, 0.02, 0.7, 18)]
        )
        manifest1 = self.write_manifest("Rep1", rep1)
        manifest2 = self.write_manifest("Rep2", rep2)
        _manifests, report = MODULE.validate_manifest_compatibility(
            ["Rep1", "Rep2"],
            {"Rep1": rep1, "Rep2": rep2},
            {"Rep1": manifest1, "Rep2": manifest2},
        )
        self.assertEqual(report["status"], "passed")

        manifest2 = self.write_manifest("Rep2_changed", rep2, ("lambda_min", 0.99))
        with self.assertRaisesRegex(ValueError, "semantic_parameters"):
            MODULE.validate_manifest_compatibility(
                ["Rep1", "Rep2"],
                {"Rep1": rep1, "Rep2": rep2},
                {"Rep1": manifest1, "Rep2": manifest2},
            )

    def test_result_must_still_match_its_run_manifest(self) -> None:
        result = self.write_result(
            "changed_after_manifest", [("g1", "intORF_altframe", "GENE1", 1, 0.01, 0.8, 20)]
        )
        manifest = self.write_manifest("Rep1_result_guard", result)
        result.write_text(result.read_text(encoding="utf-8") + "corrupt\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "row count no longer matches"):
            MODULE.validate_manifest_compatibility(
                ["Rep1"], {"Rep1": result}, {"Rep1": manifest}
            )

    def test_cli_writes_compatibility_audit(self) -> None:
        rep1 = self.write_result(
            "rep1_cli", [("g1", "intORF_altframe", "GENE1", 1, 0.01, 0.8, 20)]
        )
        rep2 = self.write_result(
            "rep2_cli", [("g1", "intORF_altframe", "GENE1", 1, 0.02, 0.7, 18)]
        )
        manifest1 = self.write_manifest("Rep1_cli", rep1)
        manifest2 = self.write_manifest("Rep2_cli", rep2)
        prefix = self.root / "integrated"
        argv = [
            "integrate_replicates.py",
            "--expected-replicate", "Rep1",
            "--expected-replicate", "Rep2",
            "--result", f"Rep1={rep1}",
            "--result", f"Rep2={rep2}",
            "--run-manifest", f"Rep1={manifest1}",
            "--run-manifest", f"Rep2={manifest2}",
            "--out-prefix", str(prefix),
        ]
        with mock.patch.object(sys, "argv", argv):
            self.assertEqual(MODULE.main(), 0)
        audit = json.loads(Path(f"{prefix}.manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(audit["compatibility"]["status"], "passed")
        self.assertEqual(set(audit["result_inputs"]), {"Rep1", "Rep2"})


if __name__ == "__main__":
    unittest.main()
