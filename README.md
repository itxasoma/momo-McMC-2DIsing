# momo-McMC-2DIsing

Equilibrium Markov-chain Monte Carlo (Metropolis single-spin flip / Glauber dynamics) simulation of the 2D Ising model with periodic boundary conditions, plus post-processing scripts to generate the figures for the MoMo computational project.

## Repository contents 
- Fortran source files for the Monte Carlo simulation (main program: `main.f90`) compiled into the executable `mdinamics`.  
- SLURM `0.run*.sh` scripts (for `sbatch`) to launch simulations on an HPC cluster.
- A posteriori binning code (`binning_program.f90`) used to compute statistical errors vs block size `m` from the time series files (see notes below).
- `PLOTS.ipynb` notebook to generate all requested figures and store outputs in the `outputs/` folder.

## Build (local)
Requirements: `gfortran` + `make`.

Compile (bash):
``make clean``
``make``

Download `timeseries_*.dat` outputs of the main simulations:
https://ubarcelona-my.sharepoint.com/:u:/g/personal/imunozal74_alumnes_ub_edu/IQB-sG4y7rWQSp-5-eiCheMmARrglJA0uoZ-1TFRJTNL8CY?e=Ij63Td
