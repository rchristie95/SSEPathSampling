# SSEPathSampling

This repository contains MATLAB scripts used for the study **"Transition Path and Interface Sampling of Stochastic Schrödinger Dynamics"**, published in *The Journal of Chemical Physics* **162**(11), 114120 (2025), [doi:10.1063/5.0246727](https://doi.org/10.1063/5.0246727). The [preprint](https://arxiv.org/abs/2411.00490) was first posted in 2024.
It implements quantum and classical trajectory sampling techniques such as
Transition Path Sampling (TPS) and Transition Interface Sampling (TIS).

## Requirements
- MATLAB
- [Optional] Parallel Computing Toolbox for running `parfor` loops

## Installation
Clone the repository and add its folder to your MATLAB path:

```bash
git clone https://github.com/rchristie95/SSEPathSampling.git
cd SSEPathSampling
```

If you lack the Parallel Computing Toolbox, replace `parfor` loops with `for` loops as noted in `ReadMe.txt`.

## Repository Structure
- **Dynamics Integrators** – Quantum integrators such as
  `SSEDynamicsTPS.m` evolve a wavefunction forward and backward around a
  chosen shooting index. Example initialization is shown in lines 1‑20 of
  the file.
- **Symmetry Transforms** – `PathTransformerSSE.m` applies time,
  parity or PT reversals, selectable via the `choice` parameter
  (see lines 1‑13).
- **Sampling Drivers** – Scripts including `SSETPSSymmetry.m` and
  `SSETransitionRateTIS.m` coordinate TPS/TIS procedures.
- **Post‑Processing** – Functions such as `adjust_histogram_weightsTIS.m`
  combine histogram data across interfaces.

Example usage of `SSETPSSymmetry` is outlined in the header comments of the file:

```matlab
SSETPSSymmetry(h2, h4, Gamma, kT, LambdaVals, windowMat, A, n_equilib, nPathsPerWindow)
```

The code saves trajectories and statistics to MAT-files that can be further
analyzed using the provided utilities.

## Running Example Scripts
To reproduce the paper's calculations, open MATLAB in this directory and run one
of the example scripts, e.g.:

```matlab
TISDetailedScript
```

Results are written to the current working directory.

## License
A license file is not provided. Please contact the authors of the original
paper for usage permissions.

## Citing

Robson Christie, Peter G. Bolhuis, and David T. Limmer, "Transition path and interface sampling of stochastic Schrödinger dynamics", *The Journal of Chemical Physics* **162**(11), 114120 (2025). https://doi.org/10.1063/5.0246727

```bibtex
@article{Christie2025Transition,
  author = {Christie, Robson and Bolhuis, Peter G. and Limmer, David T.},
  title = {Transition path and interface sampling of stochastic Schrödinger dynamics},
  journal = {The Journal of Chemical Physics},
  volume = {162},
  number = {11},
  pages = {114120},
  year = {2025},
  doi = {10.1063/5.0246727},
  url = {https://doi.org/10.1063/5.0246727}
}
```
