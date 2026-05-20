#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from ptvlib import read_tsv, strip_chr, with_chr, write_tsv


def candidate_files(cadd_dir: str, ref: str, alt: str) -> list[Path]:
    base = Path(cadd_dir)
    if len(ref) == 1 and len(alt) == 1:
        return [base / "whole_genome_SNVs.tsv.gz"]
    return [base / "gnomad.genomes.r4.0.indel.tsv.gz", base / "InDels.tsv.gz"]


def query_cadd(path: Path, chrom: str, pos: str, ref: str, alt: str) -> tuple[str, str, str]:
    if not path.exists():
        return "", "", "missing_resource"
    regions = [f"{chrom}:{pos}-{pos}", f"{strip_chr(chrom)}:{pos}-{pos}", f"{with_chr(chrom)}:{pos}-{pos}"]
    seen = set()
    for region in regions:
        if region in seen:
            continue
        seen.add(region)
        proc = subprocess.run(["tabix", str(path), region], text=True, capture_output=True)
        if proc.returncode not in {0, 1}:
            continue
        for line in proc.stdout.splitlines():
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            c_chrom, c_pos, c_ref, c_alt = parts[:4]
            if strip_chr(c_chrom) == strip_chr(chrom) and c_pos == str(pos) and c_ref == ref and c_alt == alt:
                return parts[4], parts[5], path.name
    return "", "", "not_found"


def fetch_cadd_api_json(url: str) -> object | None:
    request = urllib.request.Request(url, headers={"User-Agent": "ptv-workbox/0.1"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return []
            if attempt == 4:
                return None
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
            if attempt == 4:
                return None
        time.sleep(min(30, 2 ** attempt))
    return None


def cadd_api_match(payload: object, chrom: str, pos: str, ref: str, alt: str) -> tuple[str, str, str] | None:
    if not isinstance(payload, list):
        return None
    for item in payload:
        if not isinstance(item, dict):
            continue
        if (
            strip_chr(str(item.get("Chrom", ""))) == strip_chr(chrom)
            and str(item.get("Pos", "")) == str(pos)
            and str(item.get("Ref", "")) == ref
            and str(item.get("Alt", "")) == alt
        ):
            return str(item.get("RawScore", "")), str(item.get("PHRED", "")), "CADD_API_GRCh38-v1.7"
    return None


def query_cadd_api(chrom: str, pos: str, ref: str, alt: str) -> tuple[str, str, str]:
    if not (len(ref) == 1 and len(alt) == 1):
        return "", "", "api_snv_only"
    api_base = "https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.7"
    coord = urllib.parse.quote(f"{strip_chr(chrom)}:{pos}_{ref}_{alt}", safe=":_")
    payload = fetch_cadd_api_json(f"{api_base}/{coord}")
    match = cadd_api_match(payload, chrom, pos, ref, alt)
    if match:
        return match

    # The single-SNV endpoint can intermittently return an empty list under
    # load. The position endpoint returns all three SNV substitutions and lets
    # us verify whether the requested allele really is absent.
    pos_coord = urllib.parse.quote(f"{strip_chr(chrom)}:{pos}", safe=":")
    payload = fetch_cadd_api_json(f"{api_base}/{pos_coord}")
    match = cadd_api_match(payload, chrom, pos, ref, alt)
    if match:
        return match
    if payload is None:
        return "", "", "api_error"
    return "", "", "api_not_found"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variants", required=True)
    parser.add_argument("--cadd-dir", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--api-workers", type=int, default=2)
    args = parser.parse_args()

    rows = list(read_tsv(args.variants))
    out_rows = []
    unique_rows: dict[str, dict[str, str]] = {}
    for row in rows:
        unique_rows.setdefault(row["variant_id"], row)

    def annotate(row: dict[str, str]) -> tuple[str, tuple[str, str, str]]:
        raw = phred = source = ""
        for path in candidate_files(args.cadd_dir, row["ref"], row["alt"]):
            raw, phred, source = query_cadd(path, row["chrom"], row["pos"], row["ref"], row["alt"])
            if phred:
                break
        if not phred:
            raw, phred, source = query_cadd_api(row["chrom"], row["pos"], row["ref"], row["alt"])
        return row["variant_id"], (raw, phred, source)

    cache: dict[str, tuple[str, str, str]] = {}
    with ThreadPoolExecutor(max_workers=max(1, args.api_workers)) as executor:
        futures = [executor.submit(annotate, row) for row in unique_rows.values()]
        for future in as_completed(futures):
            key, value = future.result()
            cache[key] = value

    for row in rows:
        key = row["variant_id"]
        raw, phred, source = cache[key]
        merged = dict(row)
        merged.update({"cadd_raw": raw, "cadd_phred": phred, "cadd_source": source, "cadd_version": "1.7"})
        out_rows.append(merged)

    fields = (list(rows[0].keys()) if rows else ["variant_id"]) + ["cadd_raw", "cadd_phred", "cadd_source", "cadd_version"]
    write_tsv(args.out, out_rows, fields)
    print(f"wrote CADD annotations for {len(out_rows)} rows")
    for source, count in Counter(value[2] for value in cache.values()).most_common():
        print(f"CADD source {source}: {count} unique variants")


if __name__ == "__main__":
    main()
