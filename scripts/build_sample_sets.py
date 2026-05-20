#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

EUR_POPS = {"CEU", "FIN", "GBR", "IBS", "TSI"}


def norm(name: str) -> str:
    return "".join(ch.lower() for ch in name if ch.isalnum())


def pick(row: dict[str, str], *candidates: str) -> str:
    by_norm = {norm(k): k for k in row}
    for cand in candidates:
        key = by_norm.get(norm(cand))
        if key is not None:
            return row.get(key, "")
    return ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    with open(args.metadata, newline="") as handle:
        first = handle.readline().strip()
        handle.seek(0)
        if "\t" in first:
            reader = csv.DictReader(handle, delimiter="\t")
        else:
            header = handle.readline().strip().split()

            def whitespace_rows():
                for line in handle:
                    if line.strip():
                        values = line.strip().split()
                        yield dict(zip(header, values))

            reader = whitespace_rows()
        for row in reader:
            sample_id = pick(row, "Sample name", "sample", "sample_id", "IID")
            pop = pick(row, "Population code", "population", "pop")
            superpop = pick(row, "Superpopulation code", "superpopulation", "super_pop")
            sex = pick(row, "Sex")
            father = pick(row, "FatherID", "father", "paternal_id")
            mother = pick(row, "MotherID", "mother", "maternal_id")
            if not superpop and pop in EUR_POPS:
                superpop = "EUR"
            has_listed_parent = father not in {"", "0", "."} or mother not in {"", "0", "."}
            if sample_id and (superpop == "EUR" or pop in EUR_POPS):
                if not has_listed_parent:
                    rows.append(
                        {
                            "sample_id": sample_id,
                            "population": pop,
                            "superpopulation": "EUR",
                            "sex": sex,
                            "relatedness_filter": "parents_unlisted",
                        }
                    )

    if not rows:
        raise SystemExit(f"No EUR samples found in {args.metadata}")

    samples = out_dir / "eur_unrelated.samples.txt"
    with samples.open("w") as handle:
        for row in rows:
            handle.write(row["sample_id"] + "\n")

    metadata = out_dir / "eur_unrelated.metadata.tsv"
    with metadata.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=["sample_id", "population", "superpopulation", "sex", "relatedness_filter"],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    by_pop: dict[str, int] = {}
    for row in rows:
        by_pop[row["population"]] = by_pop.get(row["population"], 0) + 1

    qc_dir = out_dir.parent / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    with (qc_dir / "sample_selection_report.tsv").open("w") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"n_eur_unrelated\t{len(rows)}\n")
        for pop, count in sorted(by_pop.items()):
            handle.write(f"n_{pop}\t{count}\n")
        handle.write("relatedness_policy\tremove_samples_with_listed_parents\n")

    print(f"wrote {len(rows)} EUR samples to {samples}")


if __name__ == "__main__":
    main()
