# qm2pmm

**Extract unperturbed quantum-centre properties from electronic structure outputs for PMM and PyMM calculations.**

`qm2pmm` is a collection of scripts that parse the output of common quantum chemistry packages and convert the unperturbed (gas-phase) properties of a quantum centre into the input format required by the Perturbed Matrix Method (PMM) [[1]](#pmm) and its Python implementation, PyMM [[2]](#pymm-article).

---
 
## Background
 
In the Perturbed Matrix Method, the electronic properties of a quantum centre (QC) are first computed in isolation, the *unperturbed* reference. They are then perturbed by the electrostatic field of the surrounding environment, typically sampled from classical molecular dynamics. The unperturbed quantities usually include:

- Geometry of the QC reference structure
* **-g** : QC reference geometry input file. Each line should be formatted as:
```
<atom symbol> <x> <y> <z> 
```
- Permanent dipole moments of each state and transition dipole moments between states
* **-dm** : electric dipole matrix file (in a.u.). Each line should be formatted as:

```
<state_n> <state_m> <x> <y> <z>
```
where <state_n> and <state_m> are the indices of the n-th and m-th electronic state
- Electronic state energies (ground and excited states)
* **-e** : electronic energies file (in a.u.). Each line should be formatted as:

```
<state_n energy>
```

* -ch : file containing QC atomic charges for each electronic state. It should be formatted as follows:

```
<atom1 state1 charge> <atom2 state1 charge> ...
<atom1 state2 charge> <atom2 state2 charge> ...
```


Collecting these quantities by hand from different QM codes is tedious and error-prone. `qm2pmm` automates the task and provides a consistent output regardless of which package produced the data.

## Supported packages

## Requirements
 
- Python <!-- TODO: confirm -->
- NumPy <!-- TODO: list actual dependencies -->

## Installation

## Glossary

* **Quantum center** (QC): portion of the system described at quantum-mechanical level.
<br></br>

## References
<a id="pmm">[1]</a>
Aschi, M., Spezia, R., Di Nola, A., and Amadei, A. (2001). A first-principles method to model perturbed electronic wavefunctions: the effect of an external homogeneous electric field. _Chemical physics letters_, 344(3-4), 374-380, https://doi.org/10.1016/S0009-2614(01)00638-8.

<a id="pymm-article">[2]</a>
Chen, C. G., Nardi, A. N., Amadei, A., and D’Abramo M. _Journal of Chemical Theory and Computation_, __Article ASAP__, https://doi.org/10.1021/acs.jctc.2c00767.
