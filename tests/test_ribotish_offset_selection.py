#!/usr/bin/env python3
"""Tests for safe Ribo-TISH quality/offset import."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src" / "miao_orf" / "ribotish_offsets.py"
SPEC = importlib.util.spec_from_file_location("ribotish_offsets_under_test", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class RiboTishOffsetSelectionTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        root = Path(self.tempdir.name)
        self.para = root / "sample.mapped.para.py"
        self.quality = root / "sample.mapped_qual.txt"
        self.para.write_text(
            "offdict = {25: 9, 26: 10, 27: 11, 28: 12, "
            "'m0': {25: 9, 26: 10}}\n",
            encoding="utf-8",
        )
        dictionaries = [
            {25: 100, 26: 90, 27: 200, 28: 300},
            {25: [0], 26: [0], 27: [0], 28: [0]},
            {25: [0], 26: [0], 27: [0], 28: [0]},
            {
                25: [80, 10, 10],
                26: [15, 15, 60],
                27: [10, 75, 15],
                28: [75, 15, 10],
            },
            {25: [[0]], 26: [[0]], 27: [[0]], 28: [[0]]},
            {25: 10, 26: 10, 27: 10, 28: 10},
            {25: [0], 26: [0], 27: [0], 28: [0]},
            {25: [0], 26: [0], 27: [0], 28: [0]},
            {25: [1, 8, 1], 26: [1, 8, 1], 27: [1, 8, 1], 28: [1, 8, 1]},
            {25: [[0]], 26: [[0]], 27: [[0]], 28: [[0]]},
        ]
        self.quality.write_text(
            "\n".join(repr(value) for value in dictionaries) + "\n",
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_rotates_qc_frames_using_length_specific_offset(self) -> None:
        self.assertEqual(MODULE.corrected_frame_counts((10, 70, 20), 11), (70, 20, 10))
        self.assertEqual(MODULE.corrected_frame_counts((20, 20, 60), 10), (60, 20, 20))

    def test_default_contract_is_strictly_above_two_thirds_and_dominant(self) -> None:
        selected, rows = MODULE.select_ribotish_offsets(self.para, self.quality)
        self.assertEqual(selected, {27: 11, 28: 12})
        by_length = {row.length: row for row in rows}
        self.assertAlmostEqual(by_length[26].corrected_frame0_proportion, 2 / 3)
        self.assertFalse(by_length[26].selected)
        self.assertEqual(by_length[26].reason, "frame0_not_above_threshold")
        self.assertFalse(by_length[25].selected)
        self.assertEqual(by_length[25].reason, "outside_dominant_contiguous_block")

    def test_uses_match_group_not_m0_group(self) -> None:
        selected, _rows = MODULE.select_ribotish_offsets(
            self.para, self.quality, length_selection_policy="all_passing"
        )
        self.assertIn(25, selected)
        self.assertNotIn(26, selected)

    def test_all_passing_keeps_isolated_quality_lengths(self) -> None:
        selected, _rows = MODULE.select_ribotish_offsets(
            self.para, self.quality, length_selection_policy="all_passing"
        )
        self.assertEqual(selected, {25: 9, 27: 11, 28: 12})

    def test_dominant_block_uses_total_matched_read_support(self) -> None:
        selected = MODULE.dominant_contiguous_block(
            [25, 28, 29, 30], {25: 2_000_000, 28: 18_000_000, 29: 5_000_000, 30: 300_000}
        )
        self.assertEqual(selected, {28, 29, 30})

    def test_traditional_length_offset_mode_accepts_comma_or_space_separation(self) -> None:
        parsed = MODULE.parse_length_offset_specs(["28:12,29:12", "30:12"])
        selected, rows = MODULE.select_explicit_length_offsets(parsed)
        self.assertEqual(selected, {28: 12, 29: 12, 30: 12})
        self.assertTrue(all(row.selected for row in rows))
        self.assertTrue(
            all(row.selection_mode == "explicit_length_offsets" for row in rows)
        )

    def test_traditional_length_offset_mode_rejects_invalid_pairs(self) -> None:
        with self.assertRaisesRegex(ValueError, "offset for read length 28"):
            MODULE.parse_length_offset_specs(["28:28"])
        with self.assertRaisesRegex(ValueError, "duplicate read length"):
            MODULE.parse_length_offset_specs(["28:12", "28:13"])

    def test_keep_lengths_is_an_additional_whitelist(self) -> None:
        selected, rows = MODULE.select_ribotish_offsets(
            self.para, self.quality, keep_lengths=[28]
        )
        self.assertEqual(selected, {28: 12})
        by_length = {row.length: row for row in rows}
        self.assertEqual(by_length[25].reason, "not_in_keep_lengths")

    def test_explicit_override_remains_available_for_frozen_runs(self) -> None:
        selected, rows = MODULE.select_ribotish_offsets(
            self.para, quality_path=None, keep_lengths=[26, 28]
        )
        self.assertEqual(selected, {26: 10, 28: 12})
        self.assertTrue(all(row.corrected_frame0_proportion is None for row in rows))


if __name__ == "__main__":
    unittest.main()
