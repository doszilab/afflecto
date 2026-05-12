#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import time
from pathlib import Path


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Run FReSa/AFFlecto container from a params file.")

    ap.add_argument("--accession", required=True)
    ap.add_argument("--params-file", required=True)
    ap.add_argument("--run-dir", required=True)

    ap.add_argument("--use-apptainer", action="store_true")
    ap.add_argument("--apptainer-bin", default="apptainer")
    ap.add_argument("--container", required=True)
    ap.add_argument("--fresa-binary", default="./fresa")

    ap.add_argument("--binds", default="/ceph:/ceph,/tmp:/tmp,/scratch:/scratch")
    ap.add_argument("--requested-conformers", type=int, required=True)
    ap.add_argument("--timeout-min", type=int, default=100)

    ap.add_argument("--check-file", default="fresa_check.tsv")
    ap.add_argument("--meta-json", required=True)

    return ap.parse_args()


def existing_binds(bind_string: str) -> list[str]:
    binds = []

    for raw in bind_string.split(","):
        raw = raw.strip()
        if not raw:
            continue

        src = raw.split(":", 1)[0]

        if Path(src).exists():
            binds.append(raw)
        else:
            print(f"[WARN] Skipping bind because source does not exist: {raw}")

    return binds


def count_conformers(run_dir: Path) -> int:
    pdbs_dir = run_dir / "PDBs"

    if not pdbs_dir.exists():
        return 0

    n = 0

    for p in pdbs_dir.iterdir():
        if not p.is_file():
            continue
        if p.suffix.lower() != ".pdb":
            continue
        if p.name.lower() == "reference_conformation.pdb":
            continue
        n += 1

    return n


def write_check(path: Path, got: int, requested: int, status: str, elapsed_min: float = 0.0) -> None:
    with open(path, "w") as fh:
        fh.write(f"{got}\t{requested}\t{status}\t{elapsed_min}\n")


def build_command(
    args: argparse.Namespace,
    params_file_abs: Path,
) -> list[str]:
    if args.use_apptainer:
        cmd = [args.apptainer_bin, "run"]

        for bind in existing_binds(args.binds):
            cmd.extend(["-B", bind])

        cmd.extend([args.container, str(params_file_abs)])
        return cmd

    return [args.fresa_binary, str(params_file_abs)]


def main() -> None:
    args = parse_args()

    run_dir = Path(args.run_dir).resolve()
    params_file = Path(args.params_file).resolve()
    meta_json = Path(args.meta_json).resolve()

    run_dir.mkdir(parents=True, exist_ok=True)

    stdout_log = run_dir / "fresa_stdout.log"
    stderr_log = run_dir / "fresa_stderr.log"
    check_file = run_dir / args.check_file

    requested = int(args.requested_conformers)
    timeout_sec = int(args.timeout_min * 60)

    if not params_file.exists():
        raise SystemExit(f"ERROR: params file does not exist: {params_file}")

    cmd = build_command(args, params_file)

    start_ts = time.time()
    status = "pending"
    write_check(check_file, 0, requested, status, 0.0)

    print("Running FReSa command:")
    print(" ".join(cmd))
    print(f"Run directory: {run_dir}")
    print(f"Params file: {params_file}")

    try:
        with open(stdout_log, "w") as so, open(stderr_log, "w") as se:
            proc = subprocess.Popen(
                cmd,
                cwd=str(run_dir),
                stdout=so,
                stderr=se,
                preexec_fn=os.setsid,
            )

            while True:
                rc = proc.poll()
                got = count_conformers(run_dir)
                elapsed = time.time() - start_ts
                elapsed_min = round(elapsed / 60.0, 2)

                write_check(check_file, got, requested, "running", elapsed_min)

                if got >= requested and rc is None:
                    status = "ok_reached_target"
                    try:
                        os.killpg(proc.pid, signal.SIGTERM)
                        try:
                            proc.wait(timeout=10)
                        except subprocess.TimeoutExpired:
                            os.killpg(proc.pid, signal.SIGKILL)
                    except Exception:
                        pass
                    break

                if rc is not None:
                    status = "ok" if rc == 0 else f"fresa_returncode_{rc}"
                    break

                if elapsed >= timeout_sec:
                    status = "aborted_timeout"
                    try:
                        os.killpg(proc.pid, signal.SIGTERM)
                        try:
                            proc.wait(timeout=10)
                        except subprocess.TimeoutExpired:
                            os.killpg(proc.pid, signal.SIGKILL)
                    except Exception:
                        pass
                    break

                time.sleep(15)

        got = count_conformers(run_dir)

        try:
            err_txt = stderr_log.read_text(errors="ignore")
            if "Segmentation fault" in err_txt:
                status = "segfault"
            elif "std::bad_alloc" in err_txt:
                status = "memory_error"
            elif "core dumped" in err_txt:
                status = "core_dumped"
        except Exception:
            pass

        if status == "aborted_timeout" and got >= requested:
            status = "ok_timeout_reached_target"
        elif status.startswith("fresa_returncode_") and got >= requested:
            status = "ok"

    except Exception as e:
        got = count_conformers(run_dir)
        status = f"error:{e}"

    elapsed_min = round((time.time() - start_ts) / 60.0, 2)
    write_check(check_file, got, requested, status, elapsed_min)

    cleanup_files = [
        # f"{args.accession}_params.txt",
        "backtrack_report.txt",
        "fresa_check.tsv",
        "Reference_conformation.pdb",
        "report.txt",
    ]

    for fname in cleanup_files:
        fpath = run_dir / fname
        try:
            if fpath.exists():
                fpath.unlink()
        except Exception as e:
            print(f"[WARN] Could not delete {fpath}: {e}")

    meta = {
        "accession": args.accession,
        "run_dir": str(run_dir),
        "params_file": str(params_file),
        "use_apptainer": args.use_apptainer,
        "container": args.container,
        "command": cmd,
        "requested_conformers": requested,
        "generated_conformers": got,
        "status": status,
        "elapsed_min": elapsed_min,
        "stdout_log": str(stdout_log),
        "stderr_log": str(stderr_log),
        "check_file": str(check_file),
    }

    with open(meta_json, "w") as fh:
        json.dump(meta, fh, indent=2)

    if not (status == "ok" or status.startswith("ok_")):
        raise SystemExit(f"FReSa failed for {args.accession}: {status}")

    print(f"FReSa finished for {args.accession}: {status}, {got}/{requested} conformers")


if __name__ == "__main__":
    main()
