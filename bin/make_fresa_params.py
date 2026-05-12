#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path
from typing import Dict, List, Tuple


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Create FReSa/AFFlecto parameter file from structure + region TSV.")

    ap.add_argument("--accession", required=True)
    ap.add_argument("--structure", required=True)
    ap.add_argument("--regions-tsv", required=True)
    ap.add_argument("--out-params", required=True)
    ap.add_argument("--out-structure", required=True)
    ap.add_argument("--meta-json", required=True)

    ap.add_argument("--db-trs", default="/database/trs")
    ap.add_argument("--db-srs", default="/database/srs")

    ap.add_argument("--n-conformers", type=int, default=5)
    ap.add_argument("--fresa-threads", type=int, default=10)

    ap.add_argument("--writepdbs", default="true")
    ap.add_argument("--write-advanced-results", default="false")
    ap.add_argument("--write-linker-map", default="false")
    ap.add_argument("--write-backtrack-report", default="true")

    ap.add_argument("--max-attempt-by-conformation", type=int, default=1000)
    ap.add_argument("--max-backtracks-allowed", type=int, default=300)
    ap.add_argument("--nb-failures-before-recoil", type=int, default=5)
    ap.add_argument("--up-helix", type=int, default=10)
    ap.add_argument("--down-helix", type=int, default=8)

    ap.add_argument("--check-chains", default="false")
    ap.add_argument("--place-sidechains", default="false")
    ap.add_argument("--cb-for-globular", default="true")

    ap.add_argument("--use-TRS", default="true", choices=["true", "false"])

    return ap.parse_args()


def parse_regions_tsv(tsv_path: Path, accession: str) -> Tuple[List[List], List[List]]:
    classes: List[List] = []
    trs: List[List] = []

    with open(tsv_path, "r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) != 4:
                continue

            acc, s, e, t = parts

            if acc != accession:
                raise ValueError(f"TSV accession mismatch: expected {accession}, found {acc}")

            s_i = int(s)
            e_i = int(e)

            if t.upper() == "TRS":
                trs.append([s_i, e_i])
            else:
                classes.append([acc, s_i, e_i, t])

    if not classes:
        raise ValueError(f"No non-TRS regions found in {tsv_path}")

    return classes, trs


def yes_no(value: str) -> str:
    value = str(value).lower().strip()
    if value not in {"true", "false"}:
        raise ValueError(f"Expected true/false value, got: {value}")
    return value


def build_region_fields(accession: str, classes: List[List], trs: List[List]) -> Dict[str, List[str]]:
    regions: Dict[str, List[str]] = {}

    for _acc, s, e, t in classes:
        t_low = str(t).lower()

        if t_low == "flexible_prot":
            regions.setdefault("IDP", []).append(f"{accession}_0,A,A-{s},A-{e}")

        elif t_low in {"linker", "tail", "loop"}:
            key = t_low.capitalize()
            regions.setdefault(key, []).append(f"{accession}_0,A,A-{s},A-{e}")

        elif t_low == "rigid":
            continue

        else:
            raise ValueError(f"Unsupported region type for FReSa params: {t}")

    for s, e in trs:
        regions.setdefault("TRS", []).append(f"{accession}_0,A,A-{s},A-{e}")

    return regions


def write_params(
    args: argparse.Namespace,
    structure_out: Path,
    params_path: Path,
    classes: List[List],
    trs_for_params: List[List],
) -> None:
    regions = build_region_fields(args.accession, classes, trs_for_params)

    path_output_abs = str(params_path.parent.resolve()) + "/"
    structure_abs = str(structure_out.resolve())

    with open(params_path, "w") as param:
        param.write(f'DB_filename_TRS = "{args.db_trs}"\n')
        param.write(f'DB_filename_SRS = "{args.db_srs}"\n')
        param.write(f'pdb_filename = "{structure_abs}"\n')
        param.write(f"writepdbs = {yes_no(args.writepdbs)}\n")
        param.write(f"write_advanced_results = {yes_no(args.write_advanced_results)}\n")
        param.write(f"write_linker_map = {yes_no(args.write_linker_map)}\n")
        param.write(f"write_backtrack_report = {yes_no(args.write_backtrack_report)}\n")
        param.write('outpdb_dirname = "PDBs"\n')
        param.write(f'path_output = "{path_output_abs}"\n')

        for key, lst in regions.items():
            label = "TRS" if key == "TRS" else ("IDPs" if key == "IDP" else f"{key}s")
            param.write(f'{label} = "' + ";".join(lst) + '"\n')

            if key == "IDP":
                param.write("recenter = true\n")

        param.write("all_in_TRS = false\n")
        param.write(f"number_of_structures = int {args.n_conformers}\n")
        param.write(f"nb_threads = int {args.fresa_threads}\n")
        param.write(f"max_attempt_by_conformation = int {args.max_attempt_by_conformation}\n")
        param.write(f"max_backtracks_allowed = int {args.max_backtracks_allowed}\n")
        param.write(f"nb_failures_before_recoil = int {args.nb_failures_before_recoil}\n")
        param.write(f"up_helix = int {args.up_helix}\n")
        param.write(f"down_helix = int {args.down_helix}\n")
        param.write(f"check_chains = {yes_no(args.check_chains)}\n")
        param.write(f"place_sidechains = {yes_no(args.place_sidechains)}\n")
        param.write(f"CB_for_globular = {yes_no(args.cb_for_globular)}\n")


def main() -> None:
    args = parse_args()

    structure_in = Path(args.structure)
    regions_tsv = Path(args.regions_tsv)
    params_path = Path(args.out_params)
    structure_out = Path(args.out_structure)
    meta_json = Path(args.meta_json)

    if structure_in.stem != args.accession:
        raise SystemExit(f"ERROR: structure stem '{structure_in.stem}' does not match accession '{args.accession}'")

    params_path.parent.mkdir(parents=True, exist_ok=True)
    structure_out.parent.mkdir(parents=True, exist_ok=True)

    if structure_in.resolve() != structure_out.resolve():
        shutil.copy2(structure_in, structure_out)

    classes, trs = parse_regions_tsv(regions_tsv, args.accession)

    trs_for_params = trs if args.use_TRS == "true" else []

    write_params(
        args=args,
        structure_out=structure_out,
        params_path=params_path,
        classes=classes,
        trs_for_params=trs_for_params,
    )

    meta = {
        "accession": args.accession,
        "structure_in": str(structure_in),
        "structure_out": str(structure_out),
        "regions_tsv": str(regions_tsv),
        "params_file": str(params_path),
        "db_trs": args.db_trs,
        "db_srs": args.db_srs,
        "n_conformers": args.n_conformers,
        "fresa_threads": args.fresa_threads,
        "n_regions": len(classes),
        "n_trs": len(trs),
        "use_TRS": args.use_TRS,
        "n_trs_written_to_params": len(trs_for_params),
    }

    with open(meta_json, "w") as fh:
        json.dump(meta, fh, indent=2)

    print(f"Created FReSa params for {args.accession}")
    print(f"  params: {params_path}")
    print(f"  structure: {structure_out}")
    print(f"  metadata: {meta_json}")


if __name__ == "__main__":
    main()
