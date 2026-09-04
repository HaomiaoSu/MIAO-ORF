"""Installed console launcher for the repository-level MIAO pipeline."""

from __future__ import annotations

import importlib.util
import os
import sys
import sysconfig
from pathlib import Path
from types import ModuleType


def pipeline_script() -> Path:
    """Locate the orchestrator in an editable checkout or installed data directory."""
    override = os.environ.get("MIAO_ORF_PIPELINE_SCRIPT", "").strip()
    package_path = Path(__file__).resolve()
    installed_prefix = (
        package_path.parents[4]
        if len(package_path.parents) > 4
        else Path(sysconfig.get_path("data"))
    )
    candidates = []
    if override:
        candidates.append(Path(override).expanduser())
    candidates.extend(
        [
            Path(__file__).resolve().parents[2] / "miao_orf.py",
            installed_prefix / "share" / "miao-orf" / "miao_orf.py",
            Path(sysconfig.get_path("data")) / "share" / "miao-orf" / "miao_orf.py",
        ]
    )
    for candidate in candidates:
        if candidate.is_file() and candidate.stat().st_size > 0:
            return candidate
    rendered = "\n".join(f"  - {path}" for path in candidates)
    raise RuntimeError(f"cannot locate installed miao_orf.py; checked:\n{rendered}")


def load_pipeline(path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location("_miao_orf_pipeline", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load pipeline module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    module = load_pipeline(pipeline_script())
    return int(module.main())
