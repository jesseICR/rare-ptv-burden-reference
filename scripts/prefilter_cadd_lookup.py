#!/usr/bin/env python3
from __future__ import annotations

import argparse

from ptvlib import parse_float, read_tsv, write_tsv


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variants", required=True)
    parser.add_argument("--max-gnomad-af", type=float, required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rows = []
    fields = None
    for row in read_tsv(args.variants):
        fields = list(row.keys())
        af = parse_float(row.get("gnomad_af_selected"))
        if af is not None and af < args.max_gnomad_af:
            rows.append(row)
    write_tsv(args.out, rows, fields or ["variant_id"])
    print(f"wrote {len(rows)} rows for CADD lookup")


if __name__ == "__main__":
    main()
