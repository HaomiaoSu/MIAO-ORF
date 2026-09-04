#!/usr/bin/env python3
"""Static release-branding checks that do not require scientific dependencies."""

from __future__ import annotations

import json
import hashlib
import re
import unittest
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10
    import tomli as tomllib


ROOT = Path(__file__).resolve().parents[1]


class PublicBrandingTest(unittest.TestCase):
    def test_public_release_identity(self) -> None:
        expected = "1.0.0"
        files_and_patterns = {
            "miao_orf.py": r'__version__\s*=\s*"([^"]+)"',
            "src/miao_orf/__init__.py": r'__version__\s*=\s*"([^"]+)"',
            "src/miao_orf/orf_scan_transcriptome.py": r'__version__\s*=\s*"([^"]+)"',
            "src/miao_orf/ribo_metagene_qc.py": r'VERSION\s*=\s*"([^"]+)"',
            "src/miao_orf/quantify_intorf_abundance.py": r'VERSION\s*=\s*"([^"]+)"',
            "src/miao_orf/annotate_gene_cds_context.py": r'VERSION\s*=\s*"([^"]+)"',
            "src/miao_orf/psite-caller.py": r'VERSION\s*=\s*"([^"]+)"',
            "src/miao_orf/visualize_intorf_dm_results.py": r'VERSION\s*=\s*"([^"]+)"',
            "src/miao_orf/ribotish_offsets.py": r'VERSION\s*=\s*"([^"]+)"',
            "src/miao_orf/integrate_replicates.py": r'VERSION\s*=\s*"([^"]+)"',
            "src/miao_orf/visualize_replicate_integration.py": r'VERSION\s*=\s*"([^"]+)"',
        }
        for relative, pattern in files_and_patterns.items():
            with self.subTest(file=relative):
                text = (ROOT / relative).read_text(encoding="utf-8")
                match = re.search(pattern, text)
                self.assertIsNotNone(match)
                self.assertEqual(match.group(1), expected)

        # The benchmark-certified DM engine and public release are 1.0.0.
        dm_text = (ROOT / "src/miao_orf/ribo_intorf_dm_caller.py").read_text(encoding="utf-8")
        self.assertRegex(dm_text, r'__version__\s*=\s*"1\.0\.0"')

        pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        self.assertEqual(pyproject["project"]["name"], "miao-orf")
        self.assertEqual(pyproject["project"]["version"], expected)
        self.assertEqual(pyproject["project"]["license"], "MIT")
        self.assertEqual(pyproject["project"]["license-files"], ["LICENSE"])
        self.assertEqual(
            pyproject["project"]["authors"],
            [{"name": "Haomiao Su", "email": "suhaomiao@csu.edu.cn"}],
        )
        self.assertIn(
            "Development Status :: 5 - Production/Stable",
            pyproject["project"]["classifiers"],
        )
        self.assertEqual(
            pyproject["project"]["scripts"]["miao-orf"],
            "miao_orf.launcher:main",
        )
        self.assertEqual(
            pyproject["project"]["scripts"]["miao-orf-offsets"],
            "miao_orf.ribotish_offsets:main",
        )
        citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
        self.assertRegex(citation, r"(?m)^version:\s*1\.0\.0$")
        self.assertIn('email: "suhaomiao@csu.edu.cn"', citation)
        self.assertTrue((ROOT / "LICENSE").read_text(encoding="utf-8").startswith("MIT License"))
        self.assertTrue((ROOT / "environment.yml").is_file())
        self.assertTrue((ROOT / "CHANGELOG.md").is_file())
        manifest = (ROOT / "MANIFEST.in").read_text(encoding="utf-8")
        self.assertIn("include CITATION.cff", manifest)
        self.assertIn("recursive-include docs *.md", manifest)

    def test_release_lock_uses_public_names(self) -> None:
        lock = json.loads(
            (ROOT / "config/releases/benchmark_520_dual_v1.0.0.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(lock["software_version"], "1.0.0")
        self.assertEqual(lock["output_schema_version"], "1.0")
        self.assertTrue(lock["release_id"].startswith("miao-orf-1.0.0-"))

        caller = ROOT / lock["release_source"]["caller"]
        caller_sha256 = hashlib.sha256(caller.read_bytes()).hexdigest()
        self.assertEqual(
            lock["release_source"]["caller_sha256"], caller_sha256
        )

    def test_new_user_examples_use_portable_paths(self) -> None:
        checked = [
            ROOT / "docs/MANUAL.md",
            ROOT / "docs/MIAO_ORF_PIPELINE.md",
            *sorted((ROOT / "examples").rglob("*.sh")),
        ]
        for path in checked:
            with self.subTest(file=str(path.relative_to(ROOT))):
                text = path.read_text(encoding="utf-8").casefold()
                self.assertNotRegex(text, r"[a-z]:\\")
                self.assertNotIn("/home/", text)

    def test_component_machine_names_are_branded(self) -> None:
        top_level = (ROOT / "miao_orf.py").read_text(encoding="utf-8")
        self.assertIn('PROGRAM = "MIAO"', top_level)
        self.assertIn('prog="miao-orf"', top_level)

        expected = {
            "src/miao_orf/orf_scan_transcriptome.py": "miao-orf-orfscan",
            "src/miao_orf/psite-caller.py": "miao-orf-psite",
            "src/miao_orf/ribo_metagene_qc.py": "miao-orf-metagene-qc",
            "src/miao_orf/ribo_intorf_dm_caller.py": "miao-orf-intorf-dm",
            "src/miao_orf/quantify_intorf_abundance.py": "miao-orf-quantify-intorf-abundance",
            "src/miao_orf/annotate_gene_cds_context.py": "miao-orf-gene-cds-context",
            "src/miao_orf/visualize_intorf_dm_results.py": "miao-orf-visualize",
            "src/miao_orf/ribotish_offsets.py": "miao-orf-offsets",
            "src/miao_orf/integrate_replicates.py": "miao-orf-integrate-replicates",
            "src/miao_orf/visualize_replicate_integration.py": "miao-orf-visualize-replicates",
        }
        for relative, name in expected.items():
            with self.subTest(file=relative):
                self.assertIn(name, (ROOT / relative).read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
