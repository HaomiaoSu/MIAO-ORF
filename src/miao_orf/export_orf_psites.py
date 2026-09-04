#!/usr/bin/env python3
"""Optional internal exporter for single-ORF or batch ORF P-site inspection.

The implementation remains importable through the historical codon-ternary
module name so existing internal commands continue to work.
"""

from __future__ import annotations

try:
    from .export_intorf_codon_ternary import main
except ImportError:  # Direct script execution.
    from export_intorf_codon_ternary import main


if __name__ == "__main__":
    raise SystemExit(main())
