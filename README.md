# qm2pmm

**Extract unperturbed quantum-centre properties from electronic structure outputs for PMM and PyMM calculations.**

`qm2pmm` is a collection of scripts that parse the output of common quantum chemistry packages and convert the unperturbed (gas-phase) properties of a quantum centre into the input format required by the Perturbed Matrix Method (PMM) [[1]](#pmm) and its Python implementation, PyMM [[2]](#pymm-article).

---
 
## Background
 
In the Perturbed Matrix Method, the electronic properties of a quantum centre (QC) are first computed in isolation, the *unperturbed* reference. They are then perturbed by the electrostatic field of the surrounding environment, typically sampled from classical molecular dynamics. The unperturbed quantities usually include:

- **Geometry** of the QC reference structure (in Å), written to `geometry.pmm`. The first line contains the number of atoms, then each line is formatted as::
```
<atomic mass> <x> <y> <z> 
```
- **Electric dipole matrix** (in a.u.), written to `dipole_matrix.pmm`: permanent dipole moments of each state (diagonal) and transition dipole moments between states (off-diagonal). Each line is formatted as:
```
<state_n> <state_m> <x> <y> <z>
```
where <state_n> and <state_m> are the indices of the n-th and m-th electronic state.
- **Electronic state energies** (ground and excited states, in a.u.), written to `energies.pmm`. Each line is formatted as:
```
<state_n energy>
```
- **QC atomic charges** (ESP charges) for each electronic state, written to `esp_matrix.pmm`. It is formatted as follows:
```
<atom1 state1 charge> <atom2 state1 charge> ...
<atom1 state2 charge> <atom2 state2 charge> ...
```

All quantities are in the same Cartesian frame as `geometry.pmm`, which should therefore be used as the reference structure in PMM.
 
Collecting these quantities by hand from different QM codes is tedious and error-prone. `qm2pmm` automates the task and provides a consistent output regardless of which package produced the data.

## Supported packages

`qm2pmm` currently supports **DFT/TDDFT** calculations (full TDDFT or the Tamm-Dancoff approximation) with singlet excited states.

| Package | Script | Tested with |
|---------|--------|-------------|
| Q-Chem | [`tests/qchem/qchem.py`](tests/qchem/qchem.py) | Q-Chem 5.3 |
| Gaussian | [`tests/gaussian/gaussian.py`](tests/gaussian/gaussian.py) | Gaussian 16 Rev. B.01 |

## Requirements
 
- Python (>= 3.8)
- NumPy 

## Installation

No installation is needed. Clone the repository:
 
```bash
git clone git@github.com:AlessandroNNardi/qm2pmm.git
```
 
Each script is self-contained and can also be copied and used on its own.

## Usage
 
### Q-Chem
 
A single TDDFT output contains all the properties:
 
```bash
python qchem.py output.log
```

### Gaussian
 
Gaussian computes the properties of one excited state per calculation, so the script needs a ground-state output plus one TD output per excited state. The files can be given in any order, as each one is identified from its content:
 
```bash
python gaussian.py ground.log root1.log root2.log ...
```

### Options
 
The `.pmm` files are written to the current directory, or to another one with `-o`:
 
```bash
python /path/to/qm2pmm/tests/qchem/qchem.py output.log -o pmm_input/
```
 
| Option | Description |
|--------|-------------|
| `-o`, `--outdir DIR` | Directory for the output files (default: current directory) |
| `-q`, `--quiet` | Only report errors |

### Examples
 
The `tests/qchem/` and `tests/gaussian/` folders contain the QM outputs of guanine at the CAM-B3LYP/6-311++G(2d,2p) level and the `.pmm` files obtained from them. To reproduce them:
 
```bash
cd tests/qchem
python qchem.py gua_n_CAM-B3LYP_state_analysis.log
 
cd ../gaussian
python gaussian.py gua_n_CAM-B3LYP*.log
```

## Repository structure
 
```
qm2pmm/
├── src/
├── tests/
│   ├── qchem/       # qchem.py, an example Q-Chem output and its .pmm files
│   └── gaussian/    # gaussian.py, example Gaussian outputs and their .pmm files
├── LICENSE
└── README.md
```
 
## Glossary
 
* **Quantum center** (QC): portion of the system described at quantum-mechanical level.
<br></br>

## References
<a id="pmm">[1]</a>
Aschi, M., Spezia, R., Di Nola, A., and Amadei, A. (2001). A first-principles method to model perturbed electronic wavefunctions: the effect of an external homogeneous electric field. _Chemical physics letters_, 344(3-4), 374-380, https://doi.org/10.1016/S0009-2614(01)00638-8.

<a id="pymm-article">[2]</a>
Chen, C. G., Nardi, A. N., Amadei, A., and D’Abramo M. _Journal of Chemical Theory and Computation_, __Article ASAP__, https://doi.org/10.1021/acs.jctc.2c00767.
