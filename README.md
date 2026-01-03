# Molecular Modelling: Markov-chain Monte Carlo, 2D Ising model

*Authors: Manel Díaz, Adrián Llamas, Itxaso Muñoz-Aldalur*

Equilibrium Markov-chain Monte Carlo (Metropolis single-spin flip / Glauber dynamics) simulation of the 2D Ising model with periodic boundary conditions, plus post-processing scripts to generate the figures for the MoMo computational project.

## Repository contents 
- Fortran source files for the Monte Carlo simulation (main program: `main.f90`) compiled into the executable `mdinamics`.  
- SLURM `0.run*.sh` scripts (for `sbatch`) to launch simulations on an HPC cluster.
- A posteriori binning code (`binning_program.f90`) used to compute statistical errors vs block size `m` from the time series files.
- `PLOTS.ipynb` notebook to generate all requested figures and store outputs in the `outputs/` folder.
- Final report is available as *.pdf at `report/` folder.

## Build (local)
Requirements: `gfortran` + `make`.

Compile (bash):
``make clean``
``make``

### Step by step:
1) Download `timeseries_*.dat` outputs of the main simulations:

https://ubarcelona-my.sharepoint.com/:u:/g/personal/imunozal74_alumnes_ub_edu/IQB-sG4y7rWQSp-5-eiCheMmARrglJA0uoZ-1TFRJTNL8CY?e=Ij63Td

2) Place the outputs in the same folder as the rest of the code.
3) Execute `binning_program.f90` and `PLOTS.ipynb` to get the results.
