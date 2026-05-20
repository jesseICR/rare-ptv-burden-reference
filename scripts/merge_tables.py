#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv

from ptvlib import open_text


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    fieldnames = None
    n = 0
    with open_text(args.out, "wt") as out:
        writer = None
        for path in args.inputs:
            with open_text(path) as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if fieldnames is None:
                    fieldnames = reader.fieldnames or []
                    writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
                    writer.writeheader()
                for row in reader:
                    writer.writerow(row)
                    n += 1
    print(f"merged {n} rows into {args.out}")


if __name__ == "__main__":
    main()
