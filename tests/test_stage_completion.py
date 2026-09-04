#!/usr/bin/env python3
"""Tests for content-validated pipeline completion records."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "miao_orf.py"
SPEC = importlib.util.spec_from_file_location("miao_orf_completion_under_test", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class StageCompletionTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.program = self.root / "stage.py"
        self.program.write_text("print('ok')\n", encoding="utf-8")
        self.source = self.root / "input.tsv"
        self.source.write_text("id\n1\n", encoding="utf-8")
        self.output = self.root / "output.tsv"
        self.output.write_text("id\n1\n", encoding="utf-8")
        self.spec = MODULE.StageSpec(
            name="test",
            command=[sys.executable, str(self.program), "--input", str(self.source)],
            required_inputs=[self.source],
            expected_outputs=[self.output],
            log_path=self.root / "run.log",
        )

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_completion_binds_code_inputs_command_and_outputs(self) -> None:
        completion = MODULE.write_stage_completion(self.spec)
        self.assertTrue(completion.is_file())
        self.assertEqual(MODULE.stage_completion_state(self.spec), (True, "validated completion record"))

        self.output.write_text("id\n2\n", encoding="utf-8")
        complete, reason = MODULE.stage_completion_state(self.spec)
        self.assertFalse(complete)
        self.assertIn("expected_outputs", reason)

    def test_changed_input_or_program_invalidates_reuse(self) -> None:
        MODULE.write_stage_completion(self.spec)
        self.source.write_text("id\n1\n2\n", encoding="utf-8")
        self.assertIn("required_inputs", MODULE.stage_completion_state(self.spec)[1])

        MODULE.write_stage_completion(self.spec)
        self.program.write_text("print('changed')\n", encoding="utf-8")
        self.assertIn("stage_program", MODULE.stage_completion_state(self.spec)[1])

    def test_nonempty_malformed_output_is_rejected(self) -> None:
        malformed = self.root / "bad.json"
        malformed.write_text("{not-json}\n", encoding="utf-8")
        errors = MODULE.validate_output_file(malformed)
        self.assertEqual(len(errors), 1)
        self.assertIn("invalid", errors[0])

    def test_outputs_without_completion_are_stale(self) -> None:
        complete, reason = MODULE.stage_completion_state(self.spec)
        self.assertFalse(complete)
        self.assertIn("completion record is missing", reason)


if __name__ == "__main__":
    unittest.main()
