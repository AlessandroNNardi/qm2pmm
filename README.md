# qm2pmm

**Extract unperturbed quantum-centre properties from electronic structure outputs for PMM and PyMM calculations.**

`qm2pmm` is a collection of scripts that parse the output of common quantum chemistry packages and convert the unperturbed (gas-phase) properties of a quantum centre into the input format required by the Perturbed Matrix Method (PMM) and its Python implementation, PyMM.

---
 
## Background
 
In the Perturbed Matrix Method, the electronic properties of a quantum centre (QC) are first computed in isolation, the *unperturbed* reference. They are then perturbed by the electrostatic field of the surrounding environment, typically sampled from classical molecular dynamics. The unperturbed quantities usually include:

- Geometry of the QC reference structure
- Electronic state energies (ground and excited states)
- Permanent dipole moments of each state
- Transition dipole moments between states

Collecting these quantities by hand from different QM codes is tedious and error-prone. `qm2pmm` automates the task and provides a consistent output regardless of which package produced the data.

## Supported packages

## Requirements
 
- Python <!-- TODO: confirm -->
- NumPy <!-- TODO: list actual dependencies -->

## Installation
