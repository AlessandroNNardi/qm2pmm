#!/usr/bin/env python3
"""Extract unperturbed properties for PMM/PyMM from Gaussian 16 TDDFT outputs.

Usage:
    python gaussian.py GROUND.log ROOT1.log ROOT2.log ... [-o OUTDIR] [-q]

e.g. python gaussian.py gua_n_CAM-B3LYP.log gua_n_CAM-B3LYP_esp*.log

One output per state is needed: a ground-state calculation and one TD
calculation per excited state (Root=1, Root=2, ...). The file order on the
command line does not matter: each file is identified from its content.

Writes to OUTDIR (default: current directory):
    energies.pmm       total energy of each state (Hartree), ground state first
    dipole_matrix.pmm  'i j mu_x mu_y mu_z' for every pair of states (a.u.)
    esp_matrix.pmm     ESP charges, one row per state, one column per atom
                       (only if present in all outputs)
    geometry.pmm       number of atoms, then 'mass x y z' (Angstrom) per atom

Required route keywords (same functional, basis and geometry in all files):
    ground state:   #P <functional>/<basis> Pop=(MK,dipole)
    excited state:  #P <functional>/<basis> TD=(NStates=m,Root=n)
                       Pop=(MK,dipole) Density=Current
                    with m >= the highest root used
    Without NoSymm, Gaussian rotates the molecule to its standard orientation.
    All properties are read in that frame, and geometry.pmm is written in the
    same frame: use it as the reference structure for PMM.

Dipole convention
-----------------
Gaussian prints physical total dipoles (nuclear + electronic). With
Density=Current, the excited-state dipole and ESP charges come from the
relaxed density of the root (Z-vector), while Q-Chem's STS_MOM values are
unrelaxed. Dipoles are referred to the centre of mass: mu = mu_printed - q R_cm.

ESP charges are taken from the last dipole-constrained fit of each file. In an
Opt+Freq output, the Freq job repeats the fit without the constraint.

Excited-to-excited transition dipoles are not printed in a standard Gaussian
TD output, so they are written as zero. With only ground-to-excited couplings,
the sign convention of the transition dipoles does not affect PMM results.
"""

import argparse
import re
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

# Section markers in the Gaussian output
BANNER = "Entering Gaussian System"
SCF_ENERGY = "SCF Done:"
TD_ENERGY = "Total Energy, E(TD-HF/TD-DFT)"
ROOT_MARKER = "This state for optimization"
EXCITED_STATE = "Excited State"
TRANS_DIPOLES = "Ground to excited state transition electric dipole moments (Au):"
DIPOLE = "Dipole moment (field-independent basis, Debye):"
POPULATION = "Population analysis using"
ESP = "ESP charges:"
ESP_CONSTRAINT = "dipole moment will be constrained"
GEOMETRIES = ("Standard orientation:", "Input orientation:")

GEOMETRY_TOLERANCE = 1e-4    # Angstrom
ENERGY_TOLERANCE = 1e-6      # Hartree
ESP_DIPOLE_TOLERANCE = 0.01  # Debye


class ParserError(Exception):
    """An output cannot be parsed, lacks required data, or the set is inconsistent."""


# ============================================================ main function


def parse(paths):
    """Parse a ground-state output and one TD output per excited state.

    Returns a dictionary with:
    energies     (n_states,)              Hartree
    dipoles      (n_states, n_states, 3)  physical dipoles, a.u., about the
                                          centre of mass
    esp          (n_states, n_atoms)      ESP charges, or None
    symbols, coordinates (Angstrom), charge, method, basis, tda, files
    """
    outputs = [read_output(p) for p in paths]

    ground = [o for o in outputs if o["root"] is None]
    excited = sorted((o for o in outputs if o["root"] is not None), key=lambda o: o["root"])
    if len(ground) != 1:
        raise ParserError(f"expected one ground-state output (no TD), found {len(ground)}: "
                          + ", ".join(o["path"] for o in ground))
    ground = ground[0]
    roots = [o["root"] for o in excited]
    if not roots:
        raise ParserError("no TD outputs given")
    if roots != list(range(1, len(roots) + 1)):
        raise ParserError(f"excited-state roots must be 1..n without gaps or repeats, got {roots}")
    states = [ground] + excited

    # Geometry: must be the same in all files
    symbols, coords = parse_geometry(ground)
    for o in excited:
        s, c = parse_geometry(o)
        if s != symbols or np.abs(c - coords).max() > GEOMETRY_TOLERANCE:
            raise ParserError(f"geometry in {o['path']} differs from {ground['path']}")

    # Energies
    e0 = last_scf_energy(ground)
    for o in excited:
        if abs(last_scf_energy(o) - e0) > ENERGY_TOLERANCE:
            raise ParserError(f"ground-state SCF energy in {o['path']} differs from "
                              f"{ground['path']}: different level of theory or geometry?")
    energies = np.array([e0] + [td_energy(o) for o in excited])

    # Dipoles: physical, about the centre of mass
    charge = parse_charge(ground)
    Z = np.array([ELEMENTS.index(s) for s in symbols])
    masses = np.array([MASSES[z] for z in Z])
    r_cm = masses @ coords / masses.sum() / BOHR
    n = len(states)
    dipoles = np.zeros((n, n, 3))
    for k, o in enumerate(states):
        dipoles[k, k] = state_dipole(o) - charge * r_cm
    trans = transition_dipoles(excited[0], n - 1)
    for k in range(1, n):
        dipoles[0, k] = dipoles[k, 0] = trans[k - 1]

    # ESP charges: written only if present for all states
    esp = [parse_esp(o, len(symbols)) for o in states]
    if all(q is None for q in esp):
        esp = None
    elif any(q is None for q in esp):
        missing = [o["path"] for o, q in zip(states, esp) if q is None]
        raise ParserError("ESP charges missing in: " + ", ".join(missing)
                          + " (use Pop=(MK,dipole) in all calculations)")
    else:
        esp = np.array(esp)
        # The charges of each state must reproduce its dipole
        for k, o in enumerate(states):
            mu_q = esp[k] @ coords / BOHR
            mu = dipoles[k, k] + charge * r_cm
            diff = np.abs(mu_q - mu).max() * AU_TO_DEBYE
            if diff > ESP_DIPOLE_TOLERANCE:
                print(f"warning: ESP charges of state {k} ({o['path']}) differ from its "
                      f"dipole by {diff:.3f} D (fit without the dipole constraint?)",
                      file=sys.stderr)

    route = excited[0]["route"]
    method_basis = next((t for t in route.split() if "/" in t), "?/?").split("/", 1)
    return {
        "energies": energies,
        "dipoles": dipoles,
        "esp": esp,
        "symbols": symbols,
        "coordinates": coords,
        "charge": charge,
        "method": method_basis[0],
        "basis": method_basis[1],
        "tda": bool(re.search(r"\bTDA\b", route.upper())),
        "files": [o["path"] for o in states],
    }


# ============================================================ one output file


def read_output(path):
    """Read a Gaussian output and identify which state it describes."""
    path = str(path)
    with open(path, "r", errors="replace") as f:
        lines = f.read().splitlines()
    if not any(BANNER in line for line in lines[:200]):
        raise ParserError(f"'{path}' does not look like a Gaussian output")
    if not any("Normal termination" in line for line in lines[-10:]):
        raise ParserError(f"'{path}' did not terminate normally")

    out = {"path": path, "lines": lines, "route": parse_route(lines), "root": None}
    if any(EXCITED_STATE in line for line in lines):
        out["root"] = find_root(out)
        # The last population analysis must use the excited-state density
        idx = find_last(lines, POPULATION)
        if idx is None or "SCF" in lines[idx]:
            raise ParserError(f"no excited-state population analysis in '{path}' "
                              "(is Density=Current set?)")
    return out


def parse_route(lines):
    """Route section (#...), joined if split over several lines."""
    for i, line in enumerate(lines):
        if line.strip().startswith("#"):
            route = []
            for l in lines[i:]:
                if is_dashes(l):
                    break
                route.append(l.strip())
            return "".join(route)
    return ""


def find_root(output):
    """Number of the excited state whose density was used (Root=n)."""
    lines = output["lines"]
    idx = find_last(lines, ROOT_MARKER)
    if idx is not None:
        for line in reversed(lines[:idx]):
            if EXCITED_STATE in line:
                return int(line.split()[2].rstrip(":"))
    match = re.search(r"ROOT\s*=\s*(\d+)", output["route"].upper())
    if match:
        return int(match.group(1))
    raise ParserError(f"cannot determine the root of '{output['path']}'")


def parse_geometry(output):
    lines = output["lines"]
    for marker in GEOMETRIES:
        idx = find_last(lines, marker)
        if idx is not None:
            rows = table_rows(lines, idx)
            symbols = [ELEMENTS[int(r[1])] for r in rows]
            coords = np.array([[float(x) for x in r[3:6]] for r in rows])
            return symbols, coords
    raise ParserError(f"no geometry found in '{output['path']}'")


def parse_charge(output):
    for line in reversed(output["lines"]):
        match = re.search(r"Charge\s*=\s*(-?\d+)\s+Multiplicity", line)
        if match:
            return int(match.group(1))
    raise ParserError(f"molecular charge not found in '{output['path']}'")


def last_scf_energy(output):
    idx = find_last(output["lines"], SCF_ENERGY)
    if idx is None:
        raise ParserError(f"no SCF energy in '{output['path']}'")
    return float(output["lines"][idx].split("=")[1].split()[0])


def td_energy(output):
    """Total energy of the root, printed right after its 'Excited State' line."""
    idx = find_last(output["lines"], TD_ENERGY)
    if idx is None:
        raise ParserError(f"no TD total energy in '{output['path']}'")
    return float(output["lines"][idx].split("=")[1].split()[0])


def state_dipole(output):
    """Last total dipole printed (Debye) converted to atomic units."""
    lines = output["lines"]
    idx = find_last(lines, DIPOLE)
    if idx is None:
        raise ParserError(f"no dipole moment in '{output['path']}'")
    t = lines[idx + 1].split()          # X= x Y= y Z= z Tot= t
    return np.array([float(t[1]), float(t[3]), float(t[5])]) / AU_TO_DEBYE


def transition_dipoles(output, n_excited):
    lines = output["lines"]
    idx = find_last(lines, TRANS_DIPOLES)
    if idx is None:
        raise ParserError(f"no transition dipoles in '{output['path']}'")
    rows = {}
    for line in lines[idx + 2:]:          # skip the column header
        t = line.split()
        if not t or not is_int(t[0]):
            break
        rows[int(t[0])] = [float(x) for x in t[1:4]]
    missing = [k for k in range(1, n_excited + 1) if k not in rows]
    if missing:
        raise ParserError(f"transition dipoles missing for states {missing} in "
                          f"'{output['path']}' (increase NStates)")
    return np.array([rows[k] for k in range(1, n_excited + 1)])


def parse_esp(output, n_atoms):
    """ESP charges of the last dipole-constrained fit (Pop=(MK,dipole)).

    In an Opt+Freq output, the Freq job repeats the fit without the dipole
    constraint, so the last fit in the file is not always the right one.
    """
    lines = output["lines"]
    fits = [i for i, line in enumerate(lines) if ESP in line]
    if not fits:
        return None
    constrained = [i for i in fits
                   if any(ESP_CONSTRAINT in l for l in lines[max(0, i - 5):i])]
    idx = constrained[-1] if constrained else fits[-1]
    charges = []
    for line in lines[idx + 1:]:
        t = line.split()
        if len(t) == 1 and is_int(t[0]):          # column header
            continue
        if len(t) == 3 and is_int(t[0]) and not is_int(t[1]):
            charges.append(float(t[2]))
        else:
            break
    if len(charges) != n_atoms:
        raise ParserError(f"found {len(charges)} ESP charges in '{output['path']}', "
                          f"expected {n_atoms}")
    return charges


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
            f.write(f"{MASSES[ELEMENTS.index(s)]:.4f} {x:.6f} {y:.6f} {z:.6f}\n")
    written.append(path)

    return written


def print_summary(props, written):
    e = props["energies"]
    print("Files:")
    for k, path in enumerate(props["files"]):
        print(f"  state {k}: {path}")
    print(f"\nMethod / basis : {props['method']} / {props['basis']}")
    print(f"Charge         : {props['charge']}")
    print(f"Excited states : {len(e) - 1} ({'TDA' if props['tda'] else 'full TDDFT'})")
    print(f"Atoms          : {len(props['symbols'])}")
    print(f"ESP charges    : {'yes' if props['esp'] is not None else 'no'}\n")
    print("  State   Energy (Eh)        dE (eV)   |mu| (D)")
    for k in range(len(e)):
        mu = np.linalg.norm(props["dipoles"][k, k]) * AU_TO_DEBYE
        print(f"  {k:>5}   {e[k]:<16.8f}   {(e[k] - e[0]) * HARTREE_TO_EV:8.4f}   {mu:8.4f}")
    print("\nNote: excited-to-excited transition dipoles are not available "
          "from Gaussian and are written as zero.")
    print("\nWritten:")
    for path in written:
        print(f"  {path}")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Extract unperturbed properties for PMM/PyMM from Gaussian 16 "
                    "TDDFT outputs (one ground-state output + one per excited state).")
    ap.add_argument("outputs", nargs="+", help="Gaussian output files, in any order")
    ap.add_argument("-o", "--outdir", default=".",
                    help="directory for the .pmm files (default: current directory)")
    ap.add_argument("-q", "--quiet", action="store_true", help="only report errors")
    args = ap.parse_args(argv)

    try:
        props = parse(args.outputs)
        written = write_pmm(props, args.outdir)
    except (ParserError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    if not args.quiet:
        print_summary(props, written)
    return 0


if __name__ == "__main__":
    sys.exit(main())
