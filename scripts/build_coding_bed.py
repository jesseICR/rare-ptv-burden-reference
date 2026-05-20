#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

from ptvlib import open_text


def attrs_to_dict(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for match in re.finditer(r'(\S+) "([^"]*)"', text):
        attrs[match.group(1)] = match.group(2)
    return attrs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--splice-window", type=int, default=2)
    args = parser.parse_args()

    intervals: set[tuple[str, int, int]] = set()
    with open_text(args.gtf) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            chrom, _source, feature, start_s, end_s, _score, _strand, _phase, attrs_s = line.rstrip("\n").split("\t")
            attrs = attrs_to_dict(attrs_s)
            gene_type = attrs.get("gene_type") or attrs.get("gene_biotype")
            transcript_type = attrs.get("transcript_type") or attrs.get("transcript_biotype")
            if gene_type != "protein_coding" and transcript_type != "protein_coding":
                continue
            start = int(start_s)
            end = int(end_s)
            if feature == "CDS":
                intervals.add((chrom, max(0, start - 1), end))
            elif feature == "exon":
                w = args.splice_window
                intervals.add((chrom, max(0, start - 1 - w), start - 1 + w + 1))
                intervals.add((chrom, max(0, end - w), end + w))

    merged = []
    for chrom in sorted({x[0] for x in intervals}, key=lambda c: (c.replace("chr", "").isdigit() is False, c)):
        cur_start = cur_end = None
        for _chrom, start, end in sorted((x for x in intervals if x[0] == chrom), key=lambda x: (x[1], x[2])):
            if cur_start is None:
                cur_start, cur_end = start, end
            elif start <= cur_end:
                cur_end = max(cur_end, end)
            else:
                merged.append((chrom, cur_start, cur_end))
                cur_start, cur_end = start, end
        if cur_start is not None:
            merged.append((chrom, cur_start, cur_end))

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.out, "wt") as handle:
        for chrom, start, end in merged:
            if end > start:
                handle.write(f"{chrom}\t{start}\t{end}\n")
    print(f"wrote {len(merged)} merged intervals to {args.out}")


if __name__ == "__main__":
    main()
