#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional


STRUCTURE_EXTS = {".pdb", ".ent", ".cif", ".mmcif"}
TSV_EXT = ".tsv"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Discover AFFlecto input structures and optional region TSV files."
    )
    ap.add_argument("--input-dir", required=True, help="Directory containing PDB/CIF/mmCIF/ENT files and optional TSV files.")
    ap.add_argument("--region-mode", required=True, choices=["auto", "tsv"], help="auto = classify regions; tsv = use provided TSV regions.")
    ap.add_argument("--manifest", default="input_manifest.tsv", help="Output manifest TSV.")
    ap.add_argument("--summary", default="input_discovery_summary.json", help="Output summary JSON.")
    return ap.parse_args()


def collect_files(input_dir: Path, allowed_exts: set[str]) -> Dict[str, List[Path]]:
    files_by_stem: Dict[str, List[Path]] = defaultdict(list)

    for p in sorted(input_dir.iterdir()):
        if not p.is_file():
            continue
        if p.suffix.lower() in allowed_exts:
            files_by_stem[p.stem].append(p.resolve())

    return files_by_stem


def fail_duplicates(kind: str, files_by_stem: Dict[str, List[Path]]) -> None:
    duplicates = {stem: paths for stem, paths in files_by_stem.items() if len(paths) > 1}
    if not duplicates:
        return

    print(f"\nERROR: duplicate {kind} accession/stem detected.", file=sys.stderr)
    for stem, paths in duplicates.items():
        print(f"  {stem}:", file=sys.stderr)
        for p in paths:
            print(f"    - {p}", file=sys.stderr)

    raise SystemExit(2)


def one_or_none(files_by_stem: Dict[str, List[Path]], stem: str) -> Optional[Path]:
    paths = files_by_stem.get(stem, [])
    if not paths:
        return None
    return paths[0]


def write_manifest(rows: List[dict], manifest_path: Path) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    fields = ["accession", "structure_path", "tsv_path", "status"]

    with open(manifest_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_summary(summary: dict, summary_path: Path) -> None:
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=2)


def print_report(summary: dict) -> None:
    print("\nAFFlecto input discovery")
    print("------------------------")
    print(f"Input directory: {summary['input_dir']}")
    print(f"Region mode: {summary['region_mode']}")
    print(f"Structure files found: {summary['n_structure_files']}")
    print(f"TSV files found: {summary['n_tsv_files']}")

    if summary["region_mode"] == "auto":
        print(f"Usable proteins: {summary['n_usable']}")
        print(f"Ignored TSV files: {summary['n_ignored_tsv']}")
        print(f"Skipped structures: {summary['n_skipped_structures']}")
    else:
        print(f"Usable protein/TSV pairs: {summary['n_usable']}")
        print(f"Structures skipped because TSV missing: {summary['n_structures_without_tsv']}")
        print(f"TSV files ignored because structure missing: {summary['n_tsv_without_structure']}")

    print(f"Manifest: {summary['manifest']}")
    print(f"Summary JSON: {summary['summary_json']}\n")


def main() -> None:
    args = parse_args()

    input_dir = Path(args.input_dir).resolve()
    manifest_path = Path(args.manifest).resolve()
    summary_path = Path(args.summary).resolve()

    if not input_dir.exists():
        raise SystemExit(f"ERROR: input directory does not exist: {input_dir}")

    if not input_dir.is_dir():
        raise SystemExit(f"ERROR: input path is not a directory: {input_dir}")

    structures_by_stem = collect_files(input_dir, STRUCTURE_EXTS)
    tsv_by_stem = collect_files(input_dir, {TSV_EXT})

    fail_duplicates("structure", structures_by_stem)
    fail_duplicates("TSV", tsv_by_stem)

    structure_stems = set(structures_by_stem)
    tsv_stems = set(tsv_by_stem)

    rows: List[dict] = []

    if args.region_mode == "auto":
        for stem in sorted(structure_stems):
            structure = one_or_none(structures_by_stem, stem)
            rows.append(
                {
                    "accession": stem,
                    "structure_path": str(structure),
                    "tsv_path": "NA",
                    "status": "usable_auto",
                }
            )

    else:
        for stem in sorted(structure_stems):
            structure = one_or_none(structures_by_stem, stem)
            tsv = one_or_none(tsv_by_stem, stem)

            if tsv is None:
                rows.append(
                    {
                        "accession": stem,
                        "structure_path": str(structure),
                        "tsv_path": "NA",
                        "status": "skipped_missing_tsv",
                    }
                )
            else:
                rows.append(
                    {
                        "accession": stem,
                        "structure_path": str(structure),
                        "tsv_path": str(tsv),
                        "status": "usable_tsv",
                    }
                )

        for stem in sorted(tsv_stems - structure_stems):
            tsv = one_or_none(tsv_by_stem, stem)
            rows.append(
                {
                    "accession": stem,
                    "structure_path": "NA",
                    "tsv_path": str(tsv),
                    "status": "ignored_missing_structure",
                }
            )

    usable_rows = [r for r in rows if r["status"] in {"usable_auto", "usable_tsv"}]

    if not usable_rows:
        raise SystemExit(
            f"ERROR: no usable inputs found in {input_dir} with region_mode={args.region_mode}"
        )

    summary = {
        "input_dir": str(input_dir),
        "region_mode": args.region_mode,
        "n_structure_files": sum(len(v) for v in structures_by_stem.values()),
        "n_tsv_files": sum(len(v) for v in tsv_by_stem.values()),
        "n_usable": len(usable_rows),
        "n_ignored_tsv": len(tsv_by_stem) if args.region_mode == "auto" else 0,
        "n_skipped_structures": 0 if args.region_mode == "auto" else len(structure_stems - tsv_stems),
        "n_structures_without_tsv": len(structure_stems - tsv_stems),
        "n_tsv_without_structure": len(tsv_stems - structure_stems),
        "manifest": str(manifest_path),
        "summary_json": str(summary_path),
    }

    write_manifest(rows, manifest_path)
    write_summary(summary, summary_path)
    print_report(summary)


if __name__ == "__main__":
    main()
