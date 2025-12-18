#!/usr/bin/env python3
import os
import re
import sys
import tempfile
import subprocess
import resource
from typing import List, Tuple

from tqdm import tqdm

# ---------------- Stack size (ulimit -s unlimited equivalent) --------------

try:
    resource.setrlimit(
        resource.RLIMIT_STACK,
        (resource.RLIM_INFINITY, resource.RLIM_INFINITY)
    )
except Exception:
    # If we can't change it (e.g. not supported), just ignore
    pass

# ---------------- Energy parsing from xtb/gxtb output ----------------------

TOTAL_ENERGY_RE = re.compile(
    r"total\s+energy\s+([-+]?\d+\.\d+(?:[eEdD][-+]?\d+)?)",
    re.IGNORECASE,
)

TOTAL_FOOTER_RE = re.compile(
    r"^\s*total\s+([-+]?\d+\.\d+(?:[eEdD][-+]?\d+)?)\s*$",
    re.IGNORECASE | re.MULTILINE,
)

def parse_energy_from_output(output: str) -> float:
    """
    Extract energy (Eh) from xTB/g-xTB output.

    For gxtb, the final energy is printed as a footer line:
        total                          -40.49351309
    """
    # 1) gxtb footer style: "total   <Eh>"
    m_footer = TOTAL_FOOTER_RE.search(output)
    if m_footer:
        return float(m_footer.group(1))

    # 2) xtb outputs: "total energy <Eh>"
    m_total = TOTAL_ENERGY_RE.search(output)
    if m_total:
        return float(m_total.group(1))

    raise RuntimeError(
        "Could not find energy in output "
        "(no footer 'total', no 'total energy')"
    )

# ---------------- XYZ handling ---------------------------------------------


def read_xyz_blocks(lines: List[str]) -> List[List[str]]:
    """
    Read a multi-structure XYZ file and return a list of blocks.
    Each block is [natoms_line, comment_line, coord1, ..., coordN].

    Any incomplete final block is ignored.
    """
    blocks: List[List[str]] = []
    i = 0
    n_lines = len(lines)

    while i < n_lines:
        # Skip empty lines between blocks
        while i < n_lines and not lines[i].strip():
            i += 1
        if i >= n_lines:
            break

        first_line = lines[i].rstrip("\n")
        parts = first_line.split()
        try:
            natoms = int(parts[0])
        except ValueError:
            # Not a valid block start (footer etc.) -> stop
            break

        block_start = i
        i += 1

        # Comment line
        if i >= n_lines:
            break
        i += 1  # skip comment line

        coords_end = i + natoms
        if coords_end > n_lines:
            break

        block = lines[block_start:coords_end]
        blocks.append(block)

        i = coords_end

    return blocks

# ---------------- Single-point runners -------------------------------------


def run_single_point(
    method: str,
    block: List[str],
    charge: int,
    exe: str,
) -> float:
    """
    Run a single-point calculation on one XYZ block using either:
      - method='gfn2'  -> xtb --gfn 2 --sp
      - method='gxtb'  -> gxtb -c coord.xyz (+ .CHRG if needed)

    Returns energy in Eh.
    """
    method = method.lower()
    if method not in ("gfn2", "gxtb"):
        raise ValueError("method must be 'gfn2' or 'gxtb'")

    # Environment for this run: threads and memory
    env = os.environ.copy()
    env.setdefault("OMP_MAX_ACTIVE_LEVELS", "1")

    with tempfile.TemporaryDirectory() as tmpdir:
        xyz_name = "coord.xyz"
        xyz_path = os.path.join(tmpdir, xyz_name)

        # Write block to coord.xyz
        with open(xyz_path, "w") as f:
            for line in block:
                if not line.endswith("\n"):
                    line = line + "\n"
                f.write(line)

        # gfn2 via xtb
        if method == "gfn2":
            cmd = [exe, xyz_name, "--gfn", "2", "--sp", "--alpb", "water"]
            if charge != 0:
                cmd += ["-c", str(charge)]


        else:  # method == "gxtb"
            # charge via .CHRG
            if charge != 0:
                chrg_path = os.path.join(tmpdir, ".CHRG")
                with open(chrg_path, "w") as fch:
                    fch.write(f"{charge}\n")

            cmd = [
                exe,
                "-c", xyz_name
            ]

        result = subprocess.run(
            cmd,
            cwd=tmpdir,
            env=env,
            capture_output=True,
            text=True,
        )
        output = result.stdout + "\n" + result.stderr


        # Try to parse energy regardless of return code
        try:
            energy = parse_energy_from_output(output)
        except Exception as e:
            raise RuntimeError(
                f"{method} run failed for one structure.\n"
                f"Command: {' '.join(cmd)}\n"
                f"Return code: {result.returncode}\n"
                f"Error while parsing energy: {e}\n\n"
                f"=== stdout ===\n{result.stdout[:2000]}\n\n"
                f"=== stderr ===\n{result.stderr[:2000]}"
            )

        return energy

# ---------------- Reranking -------------------------------------------------


def rerank_xyz(
    method: str,
    input_path: str,
    output_path: str,
    charge: int,
    exe: str,
) -> None:
    """
    Rerank conformations in an XYZ ensemble using the specified method
    ('gfn2' or 'gxtb'), recompute energies, and write a new XYZ file with:

      line 1 : number of atoms
      line 2 : updated energy (in Eh)
      lines 3.. : coordinates
    """
    with open(input_path, "r") as f:
        lines = f.readlines()

    blocks = read_xyz_blocks(lines)
    if not blocks:
        raise RuntimeError("No valid XYZ blocks found in input file.")

    blocks_with_energy: List[Tuple[float, List[str]]] = []

    desc = f"Calculating energies ({method})"
    for block in tqdm(blocks, desc=desc):
        e = run_single_point(method, block, charge, exe)

        # Make a copy of the block and overwrite comment line with NEW energy
        new_block = block[:]  # shallow copy
        # Format similar to your example: just a number on the line
        new_block[1] = f"{e:.10f}\n"

        blocks_with_energy.append((e, new_block))

    # Sort by energy (lowest = most stable)
    blocks_with_energy.sort(key=lambda x: x[0])

    # Write out XYZ without extra blank lines
    with open(output_path, "w") as f:
        for _, block in blocks_with_energy:
            for line in block:
                if not line.endswith("\n"):
                    line = line + "\n"
                f.write(line)


def main():
    # Usage:
    #   python rerank_xtb.py gfn2 input.xyz output.xyz [charge] [xtb_executable]
    #   python rerank_xtb.py gxtb input.xyz output.xyz [charge] [gxtb_executable]
    #
    if not (4 <= len(sys.argv) <= 6):
        print(
            "Usage:\n"
            f"  {sys.argv[0]} gfn2 input.xyz output.xyz [charge] [xtb_executable]\n"
            f"  {sys.argv[0]} gxtb input.xyz output.xyz [charge] [gxtb_executable]\n\n"
            "Arguments:\n"
            "  method            'gfn2' for xTB GFN2, 'gxtb' for g-xTB\n"
            "  input.xyz         Multi-structure XYZ file to rerank\n"
            "  output.xyz        Output XYZ sorted by recomputed energies\n"
            "  charge            (optional) molecular charge (integer), default: 0\n"
            "  executable        (optional) path/name of xtb or gxtb binary\n"
        )
        sys.exit(1)

    method = sys.argv[1].lower()
    if method not in ("gfn2", "gxtb"):
        print("Error: method must be 'gfn2' or 'gxtb'")
        sys.exit(1)

    input_path = sys.argv[2]
    output_path = sys.argv[3]

    charge = 0
    if len(sys.argv) >= 5:
        charge = int(sys.argv[4])

    if len(sys.argv) == 6:
        exe = sys.argv[5]
    else:
        exe = "xtb" if method == "gfn2" else "gxtb"

    rerank_xyz(method, input_path, output_path, charge, exe)


if __name__ == "__main__":
    main()