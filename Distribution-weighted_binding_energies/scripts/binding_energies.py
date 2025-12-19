#!/usr/bin/env python3
"""
Compute binding energies and Boltzmann-weighted averages
from CREST-like conformers using xtb (GFN2-xTB) or gxtb.

- The multi-XYZ input is assumed to be ordered by energy (lowest first),
  all energies at the same level of theory. Solved using rerank_xtb.py output as input here.
- The total energy Etot of each conformer is read from the *comment line*
  (the line below the atom count). The user must ensure that this line
  contains a numeric total energy in Hartree.
- Etot is NOT recomputed by the script.
- Boltzmann populations are computed from these Etot values.
- Conformers with very low population (below a threshold) are skipped:
  no Eprot/Esub are computed for them.
- Only conformers with population >= threshold are used to compute the
  Boltzmann-weighted average binding energy.

Assumptions:
- Input is a multi-structure XYZ file (e.g. crest_conformers.xyz).
- User MUST indicate non-protein atoms using --n_sub or --sub_sel.
- Remaining atoms belong to the protein cluster.
- Charges of substrate and protein are provided by the user.
- Cluster charge is q_tot = q_sub + q_prot.
- For GFN2 (xtb): uses implicit solvent ALPB with water.
- For gxtb: single-point only, charge provided via .CHRG file.

Outputs:
- A log file with a table of:
  index, Etot (from comment), Etot_rel, Eprot, Esub, Ebind, population
- Final Boltzmann-weighted binding energy (using only conformers
  with population >= threshold).

"""
import argparse
import math
import os
import subprocess
import sys
import tempfile
import re
from pathlib import Path
from typing import List, Tuple, Optional

HARTREE_TO_KJMOL = 2625.499638  # kJ/mol per Eh
R_KJMOLK = 0.008314462618       # kJ/(mol K)

# Default population threshold: 0.0001% = 1e-6 as a fraction
POP_THRESHOLD_DEFAULT = 1.0e-6


# ---------- XYZ parsing helpers ----------

def read_multi_xyz_with_energy(
    path: Path
) -> List[Tuple[List[str], List[Tuple[float, float, float]], str, float]]:
    """
    Read a multi-structure XYZ file and return a list of:
      (symbols, coords, raw_comment_line, Etot_from_comment)

    The comment line (2nd line of each block) must contain a numeric
    total energy. We take the *last* float-like token on that line
    as Etot.

    symbols: list of element symbols (str)
    coords:  list of (x, y, z) floats
    raw_comment_line: the full comment string as in the file
    Etot_from_comment: float (unit determined by --Etot_unit)
    """
    structures = []
    with path.open() as f:
        lines = f.readlines()

    i = 0
    nlines = len(lines)
    while i < nlines:
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        try:
            nat = int(line)
        except ValueError:
            raise ValueError(f"Expected atom count at line {i+1}, got: {line!r}")

        if i + 1 >= nlines:
            raise ValueError("Unexpected end of file while reading comment line.")
        comment = lines[i+1].rstrip("\n")

        # Extract total energy from comment (last float-like token)
        tokens = comment.replace("=", " ").replace(":", " ").split()
        Etot_val: Optional[float] = None
        for tok in reversed(tokens):
            try:
                Etot_val = float(tok)
                break
            except ValueError:
                continue
        if Etot_val is None:
            raise ValueError(
                f"Could not parse a numeric energy from comment line:\n{comment!r}\n"
                f"File: {path}, block starting at line {i+1}."
            )

        atoms = []
        coords = []
        start = i + 2
        end = start + nat
        if end > nlines:
            raise ValueError("Unexpected end of file while reading atom block.")
        for j in range(start, end):
            parts = lines[j].split()
            if len(parts) < 4:
                raise ValueError(f"Bad XYZ atom line at {j+1}: {lines[j]!r}")
            sym = parts[0]
            x, y, z = map(float, parts[1:4])
            atoms.append(sym)
            coords.append((x, y, z))

        structures.append((atoms, coords, comment, Etot_val))
        i = end

    return structures

def parse_substrate_selection(sel: str, natoms: int) -> Tuple[int, int]:
    """
    Parse a substrate selection string of the form 'start:end' (1-based, inclusive end)
    into Python slice bounds (0-based start, end-exclusive).

    Examples (1-based, inclusive):
      ':60'     -> atoms 1..60     -> (0, 60)
      '40:60'   -> atoms 40..60    -> (39, 60)
      '580:'    -> atoms 580..end  -> (579, natoms)

    Returns:
      (start0, end_excl)
    """
    if ":" not in sel:
        raise ValueError(
            f"Invalid --sub_sel {sel!r}. Expected 'start:end' like ':60', '40:60', '580:'."
        )

    left, right = sel.split(":", 1)
    left = left.strip()
    right = right.strip()

    # start (1-based)
    if left == "":
        start1 = 1
    else:
        try:
            start1 = int(left)
        except ValueError:
            raise ValueError(f"Invalid --sub_sel start {left!r} in {sel!r}. Must be an integer.")

    # end (1-based, inclusive)
    if right == "":
        end1 = natoms
    else:
        try:
            end1 = int(right)
        except ValueError:
            raise ValueError(f"Invalid --sub_sel end {right!r} in {sel!r}. Must be an integer.")

    if start1 < 1:
        raise ValueError(f"Invalid --sub_sel {sel!r}: start must be >= 1 (got {start1}).")
    if end1 < start1:
        raise ValueError(f"Invalid --sub_sel {sel!r}: end must be >= start (got {end1} < {start1}).")
    if end1 > natoms:
        raise ValueError(
            f"Invalid --sub_sel {sel!r}: end={end1} exceeds number of atoms ({natoms})."
        )

    start0 = start1 - 1
    end_excl = end1  # because end is inclusive in input
    return start0, end_excl



def write_xyz(path: Path, symbols: List[str],
              coords: List[Tuple[float, float, float]],
              comment: str = "") -> None:
    """Write a single-structure XYZ file."""
    if len(symbols) != len(coords):
        raise ValueError("symbols and coords lengths differ")
    with path.open("w") as f:
        f.write(f"{len(symbols)}\n")
        f.write(comment + "\n")
        for sym, (x, y, z) in zip(symbols, coords):
            f.write(f"{sym:2s}  {x:16.8f}  {y:16.8f}  {z:16.8f}\n")


# ---------- xtb / gxtb runners ----------

def run_xtb_sp(
    xyz_path: Path,
    charge: int,
    xtb_cmd: str = "xtb",
    solvent: str = "water",
    workdir: Path = Path("."),
) -> float:
    """
    Run a single-point GFN2-xTB calculation with ALPB solvent
    using xtb, and return the total energy in Hartree.

    Command (conceptually):
        xtb geom.xyz --sp --gfn 2 --alpb water --chrg <charge>

    Energy is read primarily from the 'energy' file, falling back
    to parsing stdout if needed.
    """
    cmd = [
        xtb_cmd,
        str(xyz_path),
        "--sp",
        "--gfn", "2",
        "--alpb", solvent,
        "--chrg", str(charge),
    ]

    result = subprocess.run(
        cmd,
        cwd=str(workdir),
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        sys.stderr.write("xtb failed with return code "
                         f"{result.returncode}.\nCommand: {' '.join(cmd)}\n")
        sys.stderr.write("xtb stdout:\n" + result.stdout + "\n")
        sys.stderr.write("xtb stderr:\n" + result.stderr + "\n")
        raise RuntimeError("xtb single-point failed")

    # Try energy file first
    energy_file = workdir / "energy"
    if energy_file.exists():
        with energy_file.open() as f:
            lines = [ln.strip() for ln in f if ln.strip()]
        if not lines:
            raise RuntimeError("xtb energy file is empty")
        last = lines[-1].split()
        for token in reversed(last):
            try:
                return float(token)
            except ValueError:
                continue

    # Fallback: parse from stdout summary (:: total energy ...)
    for line in reversed(result.stdout.splitlines()):
        low = line.lower()
        if " total energy" in low:
            parts = line.replace(":", " ").replace("|", " ").split()
            # Take last float-ish token
            for token in reversed(parts):
                try:
                    return float(token)
                except ValueError:
                    continue

    raise RuntimeError("Could not parse total energy from xtb output")


def write_chrg_file(workdir: Path, charge: int) -> None:
    """
    Write .CHRG file with given integer charge in workdir.
    Used for gxtb (and could also work for xtb).
    """
    chrg_path = workdir / ".CHRG"
    with chrg_path.open("w") as f:
        f.write(f"{charge}\n")


def run_gxtb_sp(
    xyz_path: Path,
    charge: int,
    gxtb_cmd: str = "gxtb",
    workdir: Path = Path("."),
) -> float:
    """
    Run a single-point g-xTB (gxtb) calculation and return
    the total energy in Hartree.

    - gxtb is called as: gxtb -c <coords.xyz>
    - Charge is provided via a .CHRG file in the working directory.
    - We remove stale restart/coord files in the workdir because gxtb will
      try to read them and can crash if they don't match the current system.

    Energy parsing:
      - Prefer a footer line starting with 'total' (as in your example):
            total   -40.49351309
      - Also accept a line containing 'total energy' (some builds).
      - Parse from stdout+stderr (some builds print to stderr).
    """
    # ensure .CHRG is correct for this fragment
    write_chrg_file(workdir, charge)

    # gxtb can read restart/coord files implicitly if present; remove them
    # to avoid cross-talk between substrate/protein runs in the same directory.
    for fn in ("gxtbrestart", "gxtb_restart", "restart", "coord"):
        p = workdir / fn
        try:
            p.unlink()
        except FileNotFoundError:
            pass

    cmd = [gxtb_cmd, "-c", str(xyz_path)]

    result = subprocess.run(
        cmd,
        cwd=str(workdir),
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        sys.stderr.write("gxtb failed with return code "
                         f"{result.returncode}.\nCommand: {' '.join(cmd)}\n")
        sys.stderr.write("gxtb stdout:\n" + result.stdout + "\n")
        sys.stderr.write("gxtb stderr:\n" + result.stderr + "\n")
        raise RuntimeError("gxtb single-point failed")

    # Combine stdout + stderr for parsing
    lines = (result.stdout + "\n" + result.stderr).splitlines()

    total_re = re.compile(r"^\s*total\s+([-+]?\d+(?:\.\d+)?(?:[eEdD][-+]?\d+)?)\s*$", re.IGNORECASE)

    for line in reversed(lines):
        m = total_re.match(line)
        if m:
            return float(m.group(1).replace("D", "E").replace("d", "e"))

    raise RuntimeError(
        "Could not parse total energy from gxtb output. "
        "Expected either a 'total energy' line or a footer line starting with 'total'."
    )


# ---------- Boltzmann helper ----------

def compute_boltzmann_populations(energies_kj: List[float], temp: float) -> List[float]:
    """
    Given a list of total energies (kJ/mol) and temperature (K),
    compute Boltzmann populations.
    """
    if len(energies_kj) == 0:
        return []

    Emin = min(energies_kj)
    beta = 1.0 / (R_KJMOLK * temp)

    weights = []
    for E in energies_kj:
        dE = E - Emin  # kJ/mol
        w = math.exp(-beta * dE)
        weights.append(w)

    Z = sum(weights)
    if Z == 0.0:
        # edge case: all energies identical or extremely large?
        n = len(weights)
        return [1.0 / n] * n

    return [w / Z for w in weights]


# ---------- Main workflow ----------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compute binding energies and Boltzmann-weighted averages from "
            "multi-conformer XYZ (sorted by Etot), using xtb (GFN2-xTB) or gxtb.\n"
            "Etot is read from the comment line of each structure; low-population "
            "conformers are skipped."
        )
    )
    parser.add_argument(
        "xyz",
        help="Multi-structure XYZ file (e.g. crest_conformers.xyz), "
             "sorted from lowest to highest total energy.",
    )
    parser.add_argument(
        "--engine",
        choices=["gfn2", "gxtb"],
        required=True,
        help="Electronic structure backend: 'gfn2' (xtb, GFN2-xTB + ALPB water) "
             "or 'gxtb' (g-xTB binary with .CHRG).",
    )
    parser.add_argument(
        "--q_sub",
        type=int,
        required=True,
        help="Charge of the substrate fragment (integer).",
    )
    parser.add_argument(
        "--q_prot",
        type=int,
        required=True,
        help="Charge of the protein cluster fragment (integer).",
    )
    parser.add_argument(
        "--n_sub",
        type=int,
        default=60,
        help="Number of atoms belonging to the substrate (default: 60).",
    )
    parser.add_argument(
        "--sub_sel",
        type=str,
        default=None,
        help=(
            "Substrate atom selection as 'start:end' (1-based, inclusive end). "
            "Examples: ':60' (1..60), '40:60' (40..60), '580:' (580..end). "
            "If provided, this overrides --n_sub."
        ),
    )

    parser.add_argument(
        "--temp",
        type=float,
        default=298.15,
        help="Temperature in K for Boltzmann populations (default: 298.15 K).",
    )
    parser.add_argument(
        "--xtb_cmd",
        default="xtb",
        help="Name or path of xtb executable (used when --engine gfn2).",
    )
    parser.add_argument(
        "--gxtb_cmd",
        default="gxtb",
        help="Name or path of gxtb executable (used when --engine gxtb).",
    )
    parser.add_argument(
        "--log",
        default=None,
        help="Output log filename (default: <input>_<engine>_binding.log).",
    )
    parser.add_argument(
        "--Etot_unit",
        choices=["Eh", "kjmol", "kcalmol"],
        default="Eh",
        help=(
            "Unit of the Etot values extracted from the XYZ comment lines "
            "(default: Eh = Hartree)."
        ),
    )
    parser.add_argument(
        "--pop_threshold",
        type=float,
        default=POP_THRESHOLD_DEFAULT,
        help=(
            "Population threshold (fraction) below which conformers are ignored "
            "for Eprot/Esub and binding-energy averaging. "
            "Default: 1e-6 (i.e. 0.0001%%)."
        ),
    )

    args = parser.parse_args()

    xyz_path = Path(args.xyz).resolve()
    if not xyz_path.exists():
        sys.stderr.write(f"Input file not found: {xyz_path}\n")
        sys.exit(1)

    engine = args.engine
    q_sub = args.q_sub
    q_prot = args.q_prot
    q_tot = q_sub + q_prot
    n_sub = args.n_sub
    T = args.temp
    Etot_unit = args.Etot_unit
    pop_thr = args.pop_threshold

    if args.log is None:
        log_name = f"{xyz_path.stem}_{engine}_binding.log"
    else:
        log_name = args.log
    log_path = Path(log_name).resolve()

    # Read conformers + Etot from comment
    conformers = read_multi_xyz_with_energy(xyz_path)
    n_conf = len(conformers)
    
    # Determine substrate selection (either --sub_sel or default first --n_sub atoms)
    nat0 = len(conformers[0][0])
    if args.sub_sel is not None:
        sub_start0, sub_end_excl = parse_substrate_selection(args.sub_sel, nat0)
        n_sub = sub_end_excl - sub_start0
    else:
        sub_start0, sub_end_excl = 0, n_sub


    if n_conf == 0:
        sys.stderr.write("No conformers found in XYZ file.\n")
        sys.exit(1)

    print(f"Found {n_conf} conformers in {xyz_path.name}")
    if args.sub_sel is not None:
        print(f"Substrate atoms: selection '{args.sub_sel}' (1-based, inclusive) -> "
              f"{sub_start0+1}:{sub_end_excl} (1-based, inclusive), n_sub={n_sub}")
    else:
        print(f"Substrate atoms: first {n_sub}")
    print(f"Charges: q_sub = {q_sub}, q_prot = {q_prot}, q_tot = {q_tot}")
    print(f"Engine: {engine}")
    print(f"Temperature: {T:.2f} K")
    print(f"Etot_unit from comments: {Etot_unit}")
    print(f"Population threshold: {pop_thr:.2e} (fraction)")
    print(f"Log file: {log_path.name}")

    # Convert Etot to kJ/mol depending on unit
    Etot_raw: List[float] = [conf[3] for conf in conformers]  # unit as given
    if Etot_unit == "Eh":
        Etot_kj = [E * HARTREE_TO_KJMOL for E in Etot_raw]
    elif Etot_unit == "kjmol":
        Etot_kj = Etot_raw[:]
    elif Etot_unit == "kcalmol":
        Etot_kj = [E * 4.184 for E in Etot_raw]  # 1 kcal/mol ≈ 4.184 kJ/mol
    else:
        raise ValueError(f"Unsupported Etot_unit: {Etot_unit}")

    # Relative total energies (kJ/mol)
    Emin_kj = min(Etot_kj)
    Etot_rel_kj = [E - Emin_kj for E in Etot_kj]

    # Boltzmann populations from Etot_kj
    pops = compute_boltzmann_populations(Etot_kj, T)

    # Identify which conformers to evaluate (population above threshold)
    active_indices = [i for i, p in enumerate(pops) if p >= pop_thr]

    if len(active_indices) == 0:
        sys.stderr.write(
            "All conformers have population below the threshold; "
            "nothing to evaluate. Lower --pop_threshold.\n"
        )
        sys.exit(1)

    print(f"\nConformers to evaluate (pop >= {pop_thr:.2e}): "
          f"{len(active_indices)} out of {n_conf}")

    # Prepare arrays (use None for skipped conformers)
    Esub_Eh: List[Optional[float]] = [None] * n_conf
    Eprot_Eh: List[Optional[float]] = [None] * n_conf
    Ebind_kj: List[Optional[float]] = [None] * n_conf

    # Temporary working directory for all calculations
    with tempfile.TemporaryDirectory(prefix="binding_xtb_") as tmpdir_str:
        tmpdir = Path(tmpdir_str)

        for global_idx in active_indices:
            idx = global_idx  # 0-based
            symbols, coords, comment, Etot_val = conformers[idx]

            nat = len(symbols)

            # Validate selection for this conformer (in case nat differs across blocks)
            if args.sub_sel is not None:
                if sub_end_excl > nat:
                    raise ValueError(
                        f"Conformer {idx+1} has {nat} atoms, but --sub_sel {args.sub_sel!r} "
                        f"requires at least {sub_end_excl} atoms."
                    )
            else:
                if nat < n_sub:
                    raise ValueError(
                        f"Conformer {idx+1} has only {nat} atoms, but n_sub={n_sub}."
                    )

            # Split substrate and protein using the selected contiguous range
            sub_symbols = symbols[sub_start0:sub_end_excl]
            sub_coords = coords[sub_start0:sub_end_excl]

            prot_symbols = symbols[:sub_start0] + symbols[sub_end_excl:]
            prot_coords = coords[:sub_start0] + coords[sub_end_excl:]

            # Conformer-specific workdir
            conf_dir = tmpdir / f"conf_{idx+1:04d}"
            conf_dir.mkdir(parents=True, exist_ok=True)
            geom_path = conf_dir / "geom.xyz"

            # Substrate energy (frozen geometry)
            write_xyz(geom_path, sub_symbols, sub_coords,
                      comment=f"conf {idx+1} substrate")
            if engine == "gfn2":
                E_sub = run_xtb_sp(
                    geom_path,
                    charge=q_sub,
                    xtb_cmd=args.xtb_cmd,
                    solvent="water",
                    workdir=conf_dir,
                )
            else:
                E_sub = run_gxtb_sp(
                    geom_path,
                    charge=q_sub,
                    gxtb_cmd=args.gxtb_cmd,
                    workdir=conf_dir,
                )

            # Protein energy (frozen geometry)
            write_xyz(geom_path, prot_symbols, prot_coords,
                      comment=f"conf {idx+1} protein")
            if engine == "gfn2":
                E_prot = run_xtb_sp(
                    geom_path,
                    charge=q_prot,
                    xtb_cmd=args.xtb_cmd,
                    solvent="water",
                    workdir=conf_dir,
                )
            else:
                E_prot = run_gxtb_sp(
                    geom_path,
                    charge=q_prot,
                    gxtb_cmd=args.gxtb_cmd,
                    workdir=conf_dir,
                )

            Esub_Eh[idx] = E_sub
            Eprot_Eh[idx] = E_prot

            # Binding energy in Hartree and kJ/mol
            # Etot is NOT recomputed; taken from comment, so we must convert
            if Etot_unit == "Eh":
                E_tot_Eh = Etot_val
            elif Etot_unit == "kjmol":
                E_tot_Eh = Etot_val / HARTREE_TO_KJMOL
            elif Etot_unit == "kcalmol":
                E_tot_Eh = (Etot_val * 4.184) / HARTREE_TO_KJMOL
            else:
                raise ValueError(f"Unsupported Etot_unit: {Etot_unit}")

            E_bind_Eh = E_tot_Eh - (E_prot + E_sub)
            E_bind_kj = E_bind_Eh * HARTREE_TO_KJMOL
            Ebind_kj[idx] = E_bind_kj

            print(
                f"Conf {idx+1:3d}: pop = {pops[idx]:.3e}, "
                f"Etot(comment) = {Etot_val: .8f} {Etot_unit}, "
                f"E_bind = {E_bind_kj: .2f} kJ/mol"
            )

    # Compute weighted average binding energy over evaluated conformers
    sum_pop_used = 0.0
    sum_popE_used = 0.0
    for i in range(n_conf):
        if Ebind_kj[i] is not None:
            p = pops[i]
            Eb = Ebind_kj[i]
            sum_pop_used += p
            sum_popE_used += p * Eb

    if sum_pop_used == 0.0:
        sys.stderr.write("No conformers with binding energies; this should not happen.\n")
        sys.exit(1)

    weighted_Ebind = sum_popE_used / sum_pop_used

    # Write log
    with log_path.open("w") as out:
        out.write("# Binding energy analysis (xtb/gxtb) with population filtering\n")
        out.write(f"# Input XYZ : {xyz_path.name}\n")
        out.write(f"# Engine    : {engine}\n")
        out.write(f"# n_conf    : {n_conf}\n")
        out.write(f"# n_sub     : {n_sub}\n")
        if args.sub_sel is not None:
            out.write(f"# sub_sel   : {args.sub_sel}\n")
        out.write(f"# q_sub     : {q_sub}\n")
        out.write(f"# q_prot    : {q_prot}\n")
        out.write(f"# q_tot     : {q_tot}\n")
        out.write(f"# Temperature (K): {T:.2f}\n")
        out.write(f"# Etot_unit from comments: {Etot_unit}\n")
        out.write(f"# Population threshold (fraction): {pop_thr:.2e}\n")
        out.write("# Energies in Eh and kJ/mol; populations from Etot(comment)\n")
        out.write("# Conformers with pop < threshold are skipped for Eprot/Esub/Ebind.\n")
        out.write("#\n")
        header = (
            " idx  "
            "Etot_raw        "
            "Etot_unit  "
            "Etot_rel(kJ/mol)  "
            "Eprot(Eh)          "
            "Esub(Eh)           "
            "Ebind(kJ/mol)   "
            "Pop"
        )
        out.write(header + "\n")
        out.write("-" * len(header) + "\n")

        for idx in range(n_conf):
            Etot_val = Etot_raw[idx]
            rel = Etot_rel_kj[idx]
            pop = pops[idx]
            Eprot_val = Eprot_Eh[idx]
            Esub_val = Esub_Eh[idx]
            Ebind_val = Ebind_kj[idx]

            Eprot_str = f"{Eprot_val: .10f}" if Eprot_val is not None else "    NA       "
            Esub_str = f"{Esub_val: .10f}" if Esub_val is not None else "    NA       "
            Ebind_str = f"{Ebind_val: 12.3f}" if Ebind_val is not None else "     NA     "

            out.write(
                f"{idx+1:4d}  "
                f"{Etot_val: .10f}  "
                f"{Etot_unit:>7s}  "
                f"{rel: 12.3f}  "
                f"{Eprot_str}  "
                f"{Esub_str}  "
                f"{Ebind_str}  "
                f"{pop: 8.4e}\n"
            )

        out.write("\n")
        out.write(
            f"# Boltzmann-weighted average binding energy (over conformers with "
            f"pop >= {pop_thr:.2e}): {weighted_Ebind:.3f} kJ/mol\n"
        )
        out.write(f"# Effective sum of populations used: {sum_pop_used:.6f}\n")
        out.write("# E_bind = Etot(comment) - (Eprot + Esub), all at frozen complex geometries.\n")

    print()
    print(f"Boltzmann-weighted average binding energy: {weighted_Ebind:.3f} kJ/mol")
    print(f"(using conformers with pop >= {pop_thr:.2e}, sum(pop_used) = {sum_pop_used:.6f})")
    print(f"Detailed table written to: {log_path}")


if __name__ == "__main__":
    main()
