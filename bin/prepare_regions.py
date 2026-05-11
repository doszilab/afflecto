#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

from afflecto_core import (
    LinkersParams,
    accession_from_filename,
    classify_protein,
    get_plddt,
    read_classes_from_tsv,
    validate_classes,
    write_region_tsv,
)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Prepare AFFlecto region TSV in auto or TSV mode.")

    ap.add_argument("--accession", required=True)
    ap.add_argument("--structure", required=True)
    ap.add_argument("--region-mode", required=True, choices=["auto", "tsv"])
    ap.add_argument("--tsv", default=None)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--meta-json", required=True)

    ap.add_argument("--dssp-bin", default="/usr/bin/dssp")
    ap.add_argument("--dssp-version", default="4.0.4")

    ap.add_argument("--plddt-cutoff", type=float, default=0.7)
    ap.add_argument("--contact-dist", type=float, default=8.0)
    ap.add_argument("--shortest-rigid", type=int, default=3)
    ap.add_argument("--shortest-flexible", type=int, default=3)
    ap.add_argument("--short-rigid-threshold", type=int, default=25)
    ap.add_argument("--contacts-cutoff", type=float, default=10)
    ap.add_argument("--contacts-cutoff-normed", type=float, default=0.15)
    ap.add_argument("--tertcont-filt-cutoff", type=float, default=0.474)
    ap.add_argument("--max-prot-len", type=int, default=2700)
    ap.add_argument("--ca-continuity-max", type=float, default=6.0)

    return ap.parse_args()


def safe_copy_or_keep(src: Path, dst: Path) -> None:
    """
    Nextflow may stage the input TSV with the same filename as the desired output TSV.
    If src and dst resolve to the same file, do nothing.
    """
    dst.parent.mkdir(parents=True, exist_ok=True)

    try:
        if src.resolve() == dst.resolve():
            return
    except FileNotFoundError:
        pass

    shutil.copy2(src, dst)


def main() -> None:
    args = parse_args()

    accession = args.accession
    structure = Path(args.structure)
    out_tsv = Path(args.out_tsv)
    meta_json = Path(args.meta_json)

    if accession_from_filename(structure) != accession:
        raise SystemExit(
            f"ERROR: accession '{accession}' does not match structure stem '{accession_from_filename(structure)}'"
        )

    params = LinkersParams(
        plddt_cutoff=args.plddt_cutoff,
        contact_dist=args.contact_dist,
        shortest_rigid=args.shortest_rigid,
        shortest_flexible=args.shortest_flexible,
        short_rigid_threshold=args.short_rigid_threshold,
        contacts_cutoff=args.contacts_cutoff,
        contacts_cutoff_normed=args.contacts_cutoff_normed,
        tertcont_filt_cutoff=args.tertcont_filt_cutoff,
        max_prot_len=args.max_prot_len,
        ca_continuity_max=args.ca_continuity_max,
    )

    if args.region_mode == "auto":
        acc, classes, trs, meta = classify_protein(
            structure,
            params=params,
            dssp_bin=args.dssp_bin,
            dssp_version=args.dssp_version,
        )

        write_region_tsv(acc, classes, trs, out_tsv)

        meta.update(
            {
                "accession": acc,
                "structure": str(structure),
                "region_mode": "auto",
                "source_tsv": None,
                "out_tsv": str(out_tsv),
                "n_regions": len(classes),
                "n_trs": len(trs),
            }
        )

    else:
        if not args.tsv or args.tsv == "NA":
            raise SystemExit(f"ERROR: region_mode=tsv but no TSV provided for {accession}")

        tsv = Path(args.tsv)

        protein_len = len(get_plddt(structure, params))
        classes, trs = read_classes_from_tsv(tsv, accession)
        validate_classes(classes, trs, accession, protein_len)

        safe_copy_or_keep(tsv, out_tsv)

        meta = {
            "accession": accession,
            "structure": str(structure),
            "region_mode": "tsv",
            "source_tsv": str(tsv),
            "out_tsv": str(out_tsv),
            "len": protein_len,
            "n_regions": len(classes),
            "n_trs": len(trs),
        }

    meta_json.parent.mkdir(parents=True, exist_ok=True)

    with open(meta_json, "w") as fh:
        json.dump(meta, fh, indent=2)

    print(f"Prepared regions for {accession}")
    print(f"  mode: {args.region_mode}")
    print(f"  output TSV: {out_tsv}")
    print(f"  metadata: {meta_json}")


if __name__ == "__main__":
    main()
