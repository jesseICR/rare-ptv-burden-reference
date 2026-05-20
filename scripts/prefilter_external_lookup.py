#!/usr/bin/env python3
from __future__ import annotations

import argparse

from ptvlib import parse_float, read_tsv, truthy, write_tsv


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotations", required=True)
    parser.add_argument("--max-ac", type=int, required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rows = []
    fields = None
    for row in read_tsv(args.annotations):
        fields = list(row.keys())
        ac = parse_float(row.get("AC_1000G_EUR_unrelated"))
        if truthy(row.get("is_candidate_ptv")) and ac is not None and ac <= args.max_ac:
            rows.append(row)
    write_tsv(args.out, rows, fields or ["variant_id"])
    print(f"wrote {len(rows)} post-1000G-singleton rows for external lookup")


if __name__ == "__main__":
    main()
