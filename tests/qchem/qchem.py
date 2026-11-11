#!/usr/bin/env python3
"""Extract unperturbed properties for PMM/PyMM from a Q-Chem TDDFT output.

Usage:
    python qchem.py output.log [-o OUTDIR] [-q]

Writes to OUTDIR (default: current directory):
    energies.pmm       total energy of each state (Hartree), ground state first
    dipole_matrix.pmm  'i j mu_x mu_y mu_z' for every pair of states (a.u.)
    esp_matrix.pmm     ESP charges, one row per state, one column per atom
                       (only if present in the output)
    geometry.pmm       number of atoms, then 'mass x y z' (Angstrom) per atom

Required $rem options:
    CIS_N_ROOTS    n        number of excited states
    CIS_SINGLETS   TRUE
    CIS_TRIPLETS   FALSE
    STS_MOM        TRUE     dipoles between all pairs of states
    SYM_IGNORE     TRUE     keep the input orientation
    ESP_CHARGES    2        optional: ESP charges
    STATE_ANALYSIS TRUE     optional: needed for excited-state ESP charges

Dipole convention
-----------------
The STS_MOM section prints, in atomic units, matrix elements of +sum_k r_k over
the electrons (no electron charge sign, no nuclear term), from unrelaxed
densities. They are converted to physical dipoles about the centre of mass:

    mu_ij = (mu_nuc - q R_cm) delta_ij - STS_ij

with mu_nuc = sum_A Z_A R_A and q the total charge. The diagonal then equals
the total dipole printed in Q-Chem's state analysis.

If the file contains several jobs, the last one is used.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

BOHR = 0.529177210903        # Angstrom per bohr
AU_TO_DEBYE = 2.541746473
HARTREE_TO_EV = 27.211386245988

ELEMENTS = ("X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar "
            "K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr").split()
MASSES = [0.0, 1.008, 4.0026, 6.94, 9.0122, 10.81, 12.011, 14.007, 15.999, 18.998,
          20.180, 22.990, 24.305, 26.982, 28.085, 30.974, 32.06, 35.45, 39.948,
          39.098, 40.078, 44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933,
          58.693, 63.546, 65.38, 69.723, 72.630, 74.922, 78.971, 79.904, 83.798]

# Section markers in the Q-Chem output
BANNER = "Welcome to Q-Chem"
GS_ENERGY = "Total energy in the final basis set"
ES_ENERGY = "Total energy for state"
GEOMETRY = "Standard Nuclear Orientation"
GS_DIPOLE = "Electron Dipole Moments of Ground State"
ES_DIPOLE = "Electron Dipole Moments of Singlet Excited State"
GS_ES_TRANS = "Transition Moments Between Ground and Singlet Excited States"
ES_ES_TRANS = "Transition Moments Between Singlet Excited States"
GS_ESP = "Merz-Kollman ESP Net Atomic Charges"
ES_ESP = "ESP charges for excited states"


class ParserError(Exception):
    """The output cannot be parsed or lacks required data."""


# ============================================================ main function


def parse(path):
    """Parse a Q-Chem TDDFT output and return a dictionary of properties.

    energies     (n_states,)              Hartree
    dipoles      (n_states, n_states, 3)  physical dipoles, a.u., about the
                                          centre of mass
    esp          (n_states, n_atoms)      ESP charges, or None
    symbols, coordinates (Angstrom), charge, method, basis, tda
    """
    with open(path, "r", errors="replace") as f:
        lines = f.read().splitlines()

    if not any(BANNER in line for line in lines[:1000]):
        raise ParserError(f"'{path}' does not look like a Q-Chem output")

    tda = check_tddft(lines)
    rem = parse_rem(lines)
    charge = parse_charge(lines)
    symbols, coords = parse_geometry(lines)
    energies = parse_energies(lines)
    sts = parse_sts_moments(lines, len(energies))

    # Physical dipoles about the centre of mass (see the docstring at the top)
    Z = np.array([ELEMENTS.index(s.capitalize()) for s in symbols])
    masses = np.array([MASSES[z] for z in Z])
    mu_nuc = Z @ coords / BOHR
    r_cm = masses @ coords / masses.sum() / BOHR
    dipoles = -sts
    for k in range(len(energies)):
        dipoles[k, k] += mu_nuc - charge * r_cm

    return {
        "energies": energies,
        "dipoles": dipoles,
        "esp": parse_esp(lines, len(symbols), len(energies)),
        "symbols": symbols,
        "coordinates": coords,
        "charge": charge,
        "method": rem.get("METHOD", rem.get("EXCHANGE", "?")),
        "basis": rem.get("BASIS", "?"),
        "tda": tda,
    }


# ============================================================ sections


def check_tddft(lines):
    """Ensure this is a TDDFT run; return True for TDA, False for full TDDFT."""
    if not any("TDDFT" in line for line in lines):
        raise ParserError("no TDDFT calculation found; only DFT/TDDFT is supported")
    for line in reversed(lines):
        if "Excitation Energies" in line and "TDDFT" in line:
            return "TDA" in line
    return None


def parse_rem(lines):
    rem = {}
    start = find_last(lines, "$rem")
    if start is None:
        return rem
    for line in lines[start + 1:]:
        if line.strip().lower().startswith("$end"):
            break
        tokens = line.split("!")[0].replace("=", " ").split()
        if len(tokens) >= 2:
            rem[tokens[0].upper()] = tokens[1]
    return rem


def parse_charge(lines):
    start = find_last(lines, "$molecule")
    if start is not None:
        for line in lines[start + 1:]:
            tokens = line.split()
            if tokens:
                if is_int(tokens[0]):
                    return int(tokens[0])
                break
    raise ParserError("could not read the molecular charge from the $molecule block")


def parse_geometry(lines):
    idx = find_last(lines, GEOMETRY)
    if idx is None:
        raise ParserError("molecular geometry not found")
    rows = table_rows(lines, idx)
    symbols = [r[1] for r in rows]
    coords = np.array([[float(x) for x in r[2:5]] for r in rows])
    return symbols, coords


def parse_energies(lines):
    gs = find_last(lines, GS_ENERGY)
    if gs is None:
        raise ParserError("ground-state energy not found")
    energies = [float(lines[gs].split("=")[1].split()[0])]
    for line in lines[gs + 1:]:
        if ES_ENERGY in line:
            energies.append(float(line.split(":")[1].split()[0]))
    if len(energies) == 1:
        raise ParserError("no excited-state energies found")
    return np.array(energies)


def parse_sts_moments(lines, n_states):
    """Dipole matrix exactly as printed in the STS_MOM section (a.u.)."""
    n_es = n_states - 1
    sts = np.zeros((n_states, n_states, 3))

    def section(marker, expected, what):
        idx = find_last(lines, marker)
        if idx is None:
            raise ParserError(f"section '{marker}' not found (is STS_MOM = TRUE set?)")
        rows = table_rows(lines, idx)
        if len(rows) != expected:
            raise ParserError(
                f"found {len(rows)} {what}, expected {expected} from the number of "
                "excited-state energies (only singlet excited states are supported)")
        return rows

    r = section(GS_DIPOLE, 1, "ground-state dipoles")[0]
    sts[0, 0] = [float(x) for x in r[1:4]]

    for r in section(ES_DIPOLE, n_es, "excited-state dipoles"):
        k = int(r[0])
        sts[k, k] = [float(x) for x in r[1:4]]

    for r in section(GS_ES_TRANS, n_es, "ground-to-excited transition dipoles"):
        k = int(r[1])
        sts[0, k] = sts[k, 0] = [float(x) for x in r[2:5]]

    for r in section(ES_ES_TRANS, n_es * (n_es - 1) // 2,
                     "excited-to-excited transition dipoles"):
        i, j = int(r[0]), int(r[1])
        sts[i, j] = sts[j, i] = [float(x) for x in r[2:5]]

    return sts


def parse_esp(lines, n_atoms, n_states):
    """ESP charges of all states, or None if not in the output."""
    gs_idx = find_last(lines, GS_ESP)
    if gs_idx is None:
        return None

    esp = np.zeros((n_states, n_atoms))
    rows = table_rows(lines, gs_idx)
    if len(rows) != n_atoms:
        raise ParserError(f"found {len(rows)} ground-state ESP charges, expected {n_atoms}")
    esp[0] = [float(r[2]) for r in rows]

    es_idx = find_last(lines, ES_ESP)
    if es_idx is None:
        print("warning: no excited-state ESP charges (set STATE_ANALYSIS = TRUE); "
              "esp_matrix.pmm will not be written", file=sys.stderr)
        return None

    # Blocks of: a header line of state indices, then "atom q_1 q_2 ..." rows
    filled = np.zeros((n_states, n_atoms), dtype=bool)
    filled[0] = True
    columns = []
    for line in lines[es_idx + 1:]:
        tokens = line.split()
        if not tokens:
            continue
        if is_dashes(line):
            break
        if all(is_int(t) for t in tokens):
            columns = [int(t) for t in tokens]
        elif columns and is_int(tokens[0]) and len(tokens) == len(columns) + 1:
            atom = int(tokens[0]) - 1
            for state, value in zip(columns, tokens[1:]):
                esp[state, atom] = float(value)
                filled[state, atom] = True
        else:
            break

    if not filled.all():
        raise ParserError("incomplete excited-state ESP charge table")
    return esp


# ============================================================ helpers


def find_last(lines, marker):
    """Index of the last line containing marker, or None."""
    for i in range(len(lines) - 1, -1, -1):
        if marker in lines[i]:
            return i
    return None


def is_int(token):
    try:
        int(token)
        return True
    except ValueError:
        return False


def is_dashes(line):
    s = line.strip()
    return len(s) >= 3 and set(s) == {"-"}


def table_rows(lines, title_idx):
    """Rows starting with an integer in the dash-delimited table after a title."""
    rows = []
    seen_dashes = False
    for line in lines[title_idx + 1:]:
        if is_dashes(line):
            if rows:
                break
            seen_dashes = True
            continue
        tokens = line.split()
        if not seen_dashes or not tokens:
            if rows:
                break
            continue
        if is_int(tokens[0]):
            rows.append(tokens)
        elif rows:
            break
    return rows


# ============================================================ output


def write_pmm(props, outdir="."):
    """Write the .pmm files and return their paths."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    n = len(props["energies"])
    written = []

    path = outdir / "energies.pmm"
    np.savetxt(path, props["energies"], fmt="%.6f")
    written.append(path)

    path = outdir / "dipole_matrix.pmm"
    with open(path, "w") as f:
        for i in range(n):
            for j in range(n):
                x, y, z = props["dipoles"][i, j]
                f.write(f"{i} {j} {x:.6f} {y:.6f} {z:.6f} \n")
    written.append(path)

    if props["esp"] is not None:
        path = outdir / "esp_matrix.pmm"
        np.savetxt(path, props["esp"], fmt="%.6f")
        written.append(path)

    path = outdir / "geometry.pmm"
    with open(path, "w") as f:
        f.write(f"{len(props['symbols'])}\n")
        for s, (x, y, z) in zip(props["symbols"], props["coordinates"]):
            f.write(f"{MASSES[ELEMENTS.index(s.capitalize())]:.4f} {x:.6f} {y:.6f} {z:.6f}\n")
    written.append(path)

    return written


def print_summary(props, written):
    kind = {True: "TDA", False: "full TDDFT", None: "unknown"}[props["tda"]]
    e = props["energies"]
    print(f"Method / basis : {props['method']} / {props['basis']}")
    print(f"Charge         : {props['charge']}")
    print(f"Excited states : {len(e) - 1} ({kind})")
    print(f"Atoms          : {len(props['symbols'])}")
    print(f"ESP charges    : {'yes' if props['esp'] is not None else 'no'}\n")
    print("  State   Energy (Eh)        dE (eV)   |mu| (D)")
    for k in range(len(e)):
        mu = np.linalg.norm(props["dipoles"][k, k]) * AU_TO_DEBYE
        print(f"  {k:>5}   {e[k]:<16.8f}   {(e[k] - e[0]) * HARTREE_TO_EV:8.4f}   {mu:8.4f}")
    print("\nWritten:")
    for path in written:
        print(f"  {path}")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Extract unperturbed properties for PMM/PyMM from a Q-Chem TDDFT output.")
    ap.add_argument("output", help="Q-Chem output file")
    ap.add_argument("-o", "--outdir", default=".",
                    help="directory for the .pmm files (default: current directory)")
    ap.add_argument("-q", "--quiet", action="store_true", help="only report errors")
    args = ap.parse_args(argv)

    try:
        props = parse(args.output)
        written = write_pmm(props, args.outdir)
    except (ParserError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    if not args.quiet:
        print_summary(props, written)
    return 0


if __name__ == "__main__":
    sys.exit(main())
