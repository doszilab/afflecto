from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.DSSP import dssp_dict_from_pdb_file


SSE_CHARS_SHEET = set("BE")
SSE_CHARS_ALPHA = set("HGI")


@dataclass
class LinkersParams:
    plddt_cutoff: float = 0.7
    contact_dist: float = 8.0
    shortest_rigid: int = 3
    shortest_flexible: int = 3
    short_rigid_threshold: int = 25
    contacts_cutoff: float = 10
    contacts_cutoff_normed: float = 0.15
    tertcont_filt_cutoff: float = 0.474
    max_prot_len: int = 2700
    ca_continuity_max: float = 6.0


class PDBException(Exception):
    pass


def is_mmcif(path: str | Path) -> bool:
    return str(path).lower().endswith((".cif", ".mmcif"))


def accession_from_filename(path: str | Path) -> str:
    return Path(path).stem


def parse_structure(path: str | Path):
    path = str(path)

    if is_mmcif(path):
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("structure", path)
    else:
        parser = PDBParser(QUIET=True)
        with open(path, "rt") as handle:
            structure = parser.get_structure("structure", handle)

    model = structure[0]
    chains = [c for c in model]

    if len(chains) != 1:
        raise PDBException(f"Multiple chains detected ({len(chains)}). Single-chain input required: {path}")

    return model, chains[0]


def get_plddt(path: str | Path, params: LinkersParams) -> List[float]:
    model, chain = parse_structure(path)
    plddt_values: List[float] = []
    prev_idx = 0

    for residue in chain.get_residues():
        resseq = residue.get_id()[1]

        if resseq != prev_idx + 1 and prev_idx != 0:
            raise PDBException(f"Indexing discrepancy after position {prev_idx}; found {resseq}")

        try:
            bfactor = float(residue["CA"].get_bfactor())
        except KeyError:
            raise PDBException(f"Missing CA atom near residue {resseq}")

        if not (0.0 <= bfactor <= 100.0):
            raise PDBException(f"pLDDT/B-factor not in [0,100] at residue {resseq}")

        plddt_values.append(bfactor / 100.0)
        prev_idx = resseq

    if len(plddt_values) > params.max_prot_len:
        raise PDBException(f"Protein too long: len={len(plddt_values)}, max={params.max_prot_len}")

    return plddt_values


def parse_coordinates(path: str | Path, params: LinkersParams) -> List[Optional[np.ndarray]]:
    model, chain = parse_structure(path)
    coords: List[Optional[np.ndarray]] = []

    prev_ca: Optional[np.ndarray] = None
    prev_resseq: Optional[int] = None

    for residue in chain.get_residues():
        resseq = residue.get_id()[1]

        try:
            ca = residue["CA"].coord
        except KeyError:
            raise PDBException(f"Missing CA atom near residue {resseq}")

        if prev_ca is not None:
            dist = float(np.linalg.norm(ca - prev_ca))
            if dist > params.ca_continuity_max:
                raise PDBException(
                    f"Continuity discrepancy: CA distance {dist:.2f} Å between residues {prev_resseq} and {resseq}"
                )

        prev_ca = ca
        prev_resseq = resseq

        try:
            coords.append(residue["CB"].coord)
        except KeyError:
            coords.append(None)

    return coords


def ensure_cryst1_for_dssp(path: str | Path) -> Path:
    path = Path(path)

    if is_mmcif(path):
        return path

    lines = path.read_text().splitlines(keepends=True)
    coord_lines = [
        line for line in lines
        if line.startswith(("ATOM  ", "HETATM", "TER", "END"))
    ]

    if not any(line.startswith(("ATOM  ", "HETATM")) for line in coord_lines):
        raise ValueError(f"No ATOM/HETATM records found in {path}")

    tmp_path = path.with_suffix(path.suffix + ".dssp_tmp.pdb")

    cryst1 = (
        "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 "
        "P 1           1\n"
    )

    with open(tmp_path, "w") as out:
        out.write(cryst1)
        out.writelines(coord_lines)
        if not coord_lines[-1].startswith("END"):
            out.write("END\n")

    return tmp_path


def get_dssp_positions(
    path: str | Path,
    dssp_bin: Optional[str],
    dssp_version: Optional[str],
) -> List[int]:
    if not dssp_bin:
        return []

    tmp_path: Optional[Path] = None

    try:
        dssp_input = ensure_cryst1_for_dssp(path)

        if Path(dssp_input) != Path(path):
            tmp_path = Path(dssp_input)

        kwargs = {"DSSP": dssp_bin}
        if dssp_version:
            kwargs["dssp_version"] = dssp_version

        dssp = dssp_dict_from_pdb_file(str(dssp_input), **kwargs)

    finally:
        if tmp_path is not None and tmp_path.exists():
            try:
                tmp_path.unlink()
            except Exception:
                pass

    seq_chars: List[str] = []
    items = list(dssp[0].items())
    items.sort(key=lambda kv: kv[0][1][1])

    for (_, _res_id), dat in items:
        c = dat[1]
        if c in SSE_CHARS_SHEET:
            seq_chars.append("E")
        elif c in SSE_CHARS_ALPHA:
            seq_chars.append("H")
        else:
            seq_chars.append("-")

    sse_positions: List[int] = []
    i = 0

    while i < len(seq_chars):
        if seq_chars[i] in ("E", "H"):
            j = i
            while j < len(seq_chars) and seq_chars[j] == seq_chars[i]:
                j += 1
            if (j - i) >= 3:
                sse_positions.extend(range(i, j))
            i = j
        else:
            i += 1

    return sse_positions


def positions_to_ranges(pos: List[int]) -> List[List[int]]:
    if not pos:
        return []

    pos = sorted(pos)
    ranges: List[List[int]] = []
    s = pos[0]
    e = pos[0]

    for x in pos[1:]:
        if x == e + 1:
            e = x
        else:
            ranges.append([s, e])
            s = e = x

    ranges.append([s, e])
    return ranges


def find_segments_by_threshold(values: List[float], cutoff: float) -> List[List[int]]:
    segs: List[List[int]] = []
    start: Optional[int] = None

    for i, value in enumerate(values):
        if value <= cutoff:
            if start is None:
                start = i
        else:
            if start is not None and i - 1 > start:
                segs.append([start, i - 1])
            start = None

    if start is not None and len(values) - 1 > start:
        segs.append([start, len(values) - 1])

    return segs


def merge_ranges(ranges: List[List[int]], shortest_rigid_gap: int) -> List[List[int]]:
    if not ranges:
        return []

    ranges = sorted(ranges, key=lambda r: r[0])
    cur = ranges[0][:]
    out: List[List[int]] = []

    for nxt in ranges[1:]:
        if nxt[0] - cur[1] <= shortest_rigid_gap:
            cur[1] = max(cur[1], nxt[1])
        else:
            out.append(cur)
            cur = nxt[:]

    out.append(cur)
    return out


def complement_ranges(ranges: List[List[int]], end_idx: int) -> List[List[int]]:
    if not ranges:
        return [[0, end_idx]] if end_idx >= 0 else []

    ranges = sorted(ranges, key=lambda r: r[0])
    out: List[List[int]] = []
    cur_end = -1

    for s, e in ranges:
        if s > cur_end + 1:
            out.append([cur_end + 1, s - 1])
        cur_end = max(cur_end, e)

    if cur_end < end_idx:
        out.append([cur_end + 1, end_idx])

    return out


def overlaps(a: List[int], b: List[int]) -> Optional[List[int]]:
    s = max(a[0], b[0])
    e = min(a[1], b[1])
    return [s, e] if s <= e else None


def region_contacts(
    reg1: List[int],
    reg2: List[int],
    coords: List[Optional[np.ndarray]],
    contact_dist: float,
) -> int:
    n = 0

    for i in range(reg1[0], reg1[1] + 1):
        for j in range(reg2[0], reg2[1] + 1):
            if i == j:
                continue
            if coords[i] is None or coords[j] is None:
                continue
            if float(np.linalg.norm(coords[i] - coords[j])) <= contact_dist:
                n += 1

    return n


def build_rigid_contact_graph(
    rigid: List[List[int]],
    coords: List[Optional[np.ndarray]],
    params: LinkersParams,
) -> Dict[int, set[int]]:
    adj: Dict[int, set[int]] = {i: set() for i in range(len(rigid))}

    for i in range(len(rigid)):
        for j in range(i + 1, len(rigid)):
            c = region_contacts(rigid[i], rigid[j], coords, params.contact_dist)
            cutoff = params.contacts_cutoff

            if (rigid[i][1] - rigid[i][0] + 1) < params.short_rigid_threshold:
                c = c / max(1, (rigid[i][1] - rigid[i][0] + 1))
                cutoff = params.contacts_cutoff_normed
            elif (rigid[j][1] - rigid[j][0] + 1) < params.short_rigid_threshold:
                c = c / max(1, (rigid[j][1] - rigid[j][0] + 1))
                cutoff = params.contacts_cutoff_normed

            if c >= cutoff:
                adj[i].add(j)
                adj[j].add(i)

    return adj


def connected(adj: Dict[int, set[int]], a: int, b: int) -> bool:
    if a == b:
        return True

    seen = {a}
    stack = [a]

    while stack:
        u = stack.pop()
        for v in adj.get(u, ()):
            if v == b:
                return True
            if v not in seen:
                seen.add(v)
                stack.append(v)

    return False


def classify_protein(
    path: str | Path,
    params: LinkersParams,
    dssp_bin: Optional[str],
    dssp_version: Optional[str],
) -> Tuple[str, List[List], List[List], Dict]:
    accession = accession_from_filename(path)

    plddt = get_plddt(path, params)
    coords = parse_coordinates(path, params)
    length = len(plddt)

    flexible_unmerged = find_segments_by_threshold(plddt, params.plddt_cutoff)
    flexible_unfiltered = merge_ranges(flexible_unmerged, params.shortest_rigid)
    flexible: List[List[int]] = [
        fu for fu in flexible_unfiltered
        if (fu[1] - fu[0] + 1) >= params.shortest_flexible
    ]

    rigid = complement_ranges(flexible, length - 1)

    if len(rigid) >= 2 and len(flexible) != 0 and (rigid[-1][1] - rigid[-1][0] + 1) < params.shortest_rigid:
        del_reg = rigid.pop()
        flexible[-1][1] = del_reg[1]

    trs_regions: List[List[int]] = []
    sse_pos = get_dssp_positions(path, dssp_bin=dssp_bin, dssp_version=dssp_version)
    sse_ranges = positions_to_ranges(sse_pos)

    if sse_ranges:
        for sse in sse_ranges:
            if (sse[1] - sse[0] + 1) < 3:
                continue

            cont = 0
            for other in sse_ranges:
                if other is sse:
                    continue
                cont += region_contacts(sse, other, coords, params.contact_dist)

            perpos = round(cont / max(1, (sse[1] - sse[0] + 1)), 3)

            if perpos <= params.tertcont_filt_cutoff:
                flexible.append(sse[:])
                trs_regions.append(sse[:])

        flexible = merge_ranges(sorted(flexible, key=lambda x: x[0]), params.shortest_rigid)
        rigid = complement_ranges(flexible, length - 1)

        newly_flex = []
        sse_ranges_sorted = sorted(sse_ranges, key=lambda r: r[0])

        for rig in rigid:
            has_sse = any(overlaps(rig, s) for s in sse_ranges_sorted)
            if not has_sse:
                newly_flex.append(rig)

        if newly_flex:
            flexible.extend(newly_flex)
            flexible = merge_ranges(sorted(flexible, key=lambda x: x[0]), params.shortest_rigid)
            rigid = complement_ranges(flexible, length - 1)

    classes_out: List[List] = []
    trs_sorted = sorted(trs_regions, key=lambda r: r[0])

    if len(flexible) == 1 and len(rigid) == 0:
        classes_out.append([accession, flexible[0][0] + 1, flexible[0][1] + 1, "flexible_prot"])
        return accession, classes_out, trs_sorted, {"len": length}

    if len(rigid) == 1 and len(flexible) == 0:
        classes_out.append([accession, rigid[0][0] + 1, rigid[0][1] + 1, "rigid"])
        return accession, classes_out, trs_sorted, {"len": length}

    rigid_graph = build_rigid_contact_graph(rigid, coords, params)
    rigid_out_pos: Dict[int, int] = {}

    def add_block(is_rigid: bool, idx: int) -> None:
        if is_rigid:
            s, e = rigid[idx]
            classes_out.append([accession, s + 1, e + 1, "rigid"])
            rigid_out_pos[idx] = len(classes_out) - 1
            return

        s, e = flexible[idx]

        if s == 0 or e == length - 1:
            classes_out.append([accession, s + 1, e + 1, "tail"])
            return

        idx_n = None
        idx_c = None

        for r_i_tmp, rig in enumerate(rigid):
            if rig[1] == s - 1:
                idx_n = r_i_tmp
            if rig[0] == e + 1:
                idx_c = r_i_tmp

        if idx_n is None or idx_c is None:
            classes_out.append([accession, s + 1, e + 1, "linker"])
            return

        c = region_contacts(rigid[idx_n], rigid[idx_c], coords, params.contact_dist)
        cutoff = params.contacts_cutoff

        if (rigid[idx_n][1] - rigid[idx_n][0] + 1) < params.short_rigid_threshold:
            c = c / max(1, (rigid[idx_n][1] - rigid[idx_n][0] + 1))
            cutoff = params.contacts_cutoff_normed
        elif (rigid[idx_c][1] - rigid[idx_c][0] + 1) < params.short_rigid_threshold:
            c = c / max(1, (rigid[idx_c][1] - rigid[idx_c][0] + 1))
            cutoff = params.contacts_cutoff_normed

        if c >= cutoff:
            t = "loop"
        else:
            t = "loop" if connected(rigid_graph, idx_n, idx_c) else "linker"

        if t == "loop":
            loop_len = e - s + 1
            if loop_len < 5:
                if idx_n is not None and s > 0 and rigid[idx_n][1] == s - 1:
                    if (rigid[idx_n][1] - rigid[idx_n][0] + 1) > 1:
                        s -= 1
                        flexible[idx][0] -= 1
                        rigid[idx_n][1] -= 1

                        pos = rigid_out_pos.get(idx_n)
                        if pos is not None:
                            classes_out[pos][2] -= 1

                if idx_c is not None and e < length - 1 and rigid[idx_c][0] == e + 1:
                    if (rigid[idx_c][1] - rigid[idx_c][0] + 1) > 1:
                        e += 1
                        flexible[idx][1] += 1
                        rigid[idx_c][0] += 1

        classes_out.append([accession, s + 1, e + 1, t])

    r_i = 0
    f_i = 0

    if rigid and rigid[0][0] == 0:
        add_block(True, 0)
        r_i = 1

    while True:
        candidates = []

        if r_i < len(rigid):
            candidates.append((rigid[r_i][0], "R"))
        if f_i < len(flexible):
            candidates.append((flexible[f_i][0], "F"))

        if not candidates:
            break

        typ = min(candidates, key=lambda x: x[0])[1]

        if typ == "R":
            add_block(True, r_i)
            r_i += 1
        else:
            add_block(False, f_i)
            f_i += 1

    return accession, classes_out, trs_sorted, {"len": length}


def read_classes_from_tsv(tsv_path: str | Path, expected_acc: str) -> Tuple[List[List], List[List]]:
    classes = []
    trs = []

    with open(tsv_path, "r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) != 4:
                continue

            acc, s, e, t = parts

            if expected_acc and acc != expected_acc:
                continue

            s_i = int(s)
            e_i = int(e)

            if t.upper() == "TRS":
                trs.append([s_i - 1, e_i - 1])
            else:
                classes.append([acc, s_i, e_i, t])

    return classes, trs

def trs_fully_contained_in_non_rigid(trs_region: List[int], classes: List[List]) -> bool:
    trs_s0, trs_e0 = trs_region

    for _acc, s, e, t in classes:
        if str(t).lower() == "rigid":
            continue

        class_s0 = int(s) - 1
        class_e0 = int(e) - 1

        if class_s0 <= trs_s0 and trs_e0 <= class_e0:
            return True

    return False

def validate_classes(
    classes: List[List],
    trs: List[List],
    accession: str,
    protein_len: int,
) -> None:
    allowed = {"rigid", "tail", "loop", "linker", "flexible_prot"}

    if not classes:
        raise ValueError(f"[{accession}] TSV contains no non-TRS regions")

    for acc, s, e, t in classes:
        if acc != accession:
            raise ValueError(f"[{accession}] TSV accession mismatch: {acc}")

        if t not in allowed:
            raise ValueError(f"[{accession}] Unknown region type: {t}")

        if not (1 <= int(s) <= int(e) <= protein_len):
            raise ValueError(f"[{accession}] Invalid region coordinates: {s}-{e}; length={protein_len}")

    for s0, e0 in trs:
        if not (0 <= s0 <= e0 < protein_len):
            raise ValueError(f"[{accession}] Invalid TRS coordinates: {s0 + 1}-{e0 + 1}; length={protein_len}")

        if not trs_fully_contained_in_non_rigid([s0, e0], classes):
            raise ValueError(
                f"[{accession}] Invalid TRS region {s0 + 1}-{e0 + 1}: "
                "TRS regions must be fully contained within a non-rigid region "
                "(tail, loop, linker, or flexible_prot)."
            )


    last_end = max(int(c[2]) for c in classes)

    if last_end != protein_len:
        raise ValueError(f"[{accession}] TSV ends at {last_end}, but protein length is {protein_len}")


def write_region_tsv(accession: str, classes: List[List], trs: List[List], out_tsv: str | Path) -> None:
    out_tsv = Path(out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    with open(out_tsv, "w") as fho:
        for acc, s, e, t in classes:
            fho.write(f"{acc}\t{s}\t{e}\t{t}\n")

        for s, e in trs:
            fho.write(f"{accession}\t{s + 1}\t{e + 1}\tTRS\n")
