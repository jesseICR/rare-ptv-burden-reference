#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import os
import re
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


PTV_CONSEQUENCES = {
    "stop_gained",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
}


def open_text(path: str | os.PathLike, mode: str = "rt"):
    path = str(path)
    if path.endswith((".gz", ".bgz")):
        return gzip.open(path, mode)
    return open(path, mode, newline="" if "t" in mode else None)


def read_tsv(path: str | os.PathLike) -> Iterator[Dict[str, str]]:
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            yield {k: ("" if v is None else v) for k, v in row.items()}


def write_tsv(path: str | os.PathLike, rows: Iterable[Dict[str, object]], fieldnames: Sequence[str]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(fieldnames), lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: "" if row.get(k) is None else row.get(k) for k in fieldnames})


def parse_info(info: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not info or info == ".":
        return out
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
        else:
            out[item] = "1"
    return out


def variant_key(chrom: str, pos: str | int, ref: str, alt: str) -> str:
    return f"{chrom}:{pos}:{ref}:{alt}"


def strip_chr(chrom: str) -> str:
    return chrom[3:] if chrom.startswith("chr") else chrom


def with_chr(chrom: str) -> str:
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def parse_float(value: object) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip()
    if text in {"", ".", "NA", "nan", "None"}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def clean_gene_id(gene_id: str) -> str:
    return re.sub(r"\.\d+$", "", gene_id or "")


def split_consequences(value: str) -> List[str]:
    parts: List[str] = []
    for chunk in (value or "").replace("&", ",").split(","):
        chunk = chunk.strip()
        if chunk:
            parts.append(chunk)
    return parts


def is_ptv_consequence(value: str) -> bool:
    return any(c in PTV_CONSEQUENCES for c in split_consequences(value))


def parse_gt(gt: str) -> Tuple[int, int, bool]:
    if not gt or gt in {".", "./.", ".|."}:
        return 0, 0, False
    alleles = re.split(r"[\/|]", gt)
    called = [a for a in alleles if a != "."]
    if not called:
        return 0, 0, False
    alt_count = sum(1 for a in called if a not in {"0", "."})
    return alt_count, len(called), True


def format_float(value: Optional[float], digits: int = 8) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}g}"
