# Molecular Modelling: Markov-chain Monte Carlo, 2D Ising model

*Authors: Manel Díaz, Adrián Llamas, Itxaso Muñoz-Aldalur*

Equilibrium Markov-chain Monte Carlo (Metropolis single-spin flip / Glauber dynamics) simulation of the 2D Ising model with periodic boundary conditions, plus post-processing scripts to generate the figures for the MoMo computational project.

## Repository contents 
- Fortran source files for the Monte Carlo simulation (main program: `main.f90`) compiled into the executable `mdinamics`.
- SLURM `0.run*.sh` scripts (for `sbatch`) to launch simulations on an HPC cluster.
- A posteriori binning code (`binning_program.f90`) used to compute statistical errors vs block size `m` from the time series files.
- Jackknife resampling (see `jackknife/` folder) is used to compute the heat capacity and magnetic susceptibility, as an extra chapter.
- `PLOTS.ipynb` notebook to generate all requested figures and store outputs in the `outputs/` folder.
- Final report is available as a PDF in the `report/` folder.

## Build (local)
Requirements: `gfortran` + `make`.

Compile (bash):

`make clean`  

`make`

## Step by step
1) Download `timeseries_*.dat` outputs of the main simulations:

https://ubarcelona-my.sharepoint.com/:u:/g/personal/imunozal74_alumnes_ub_edu/IQB-sG4y7rWQSp-5-eiCheMmARrglJA0uoZ-1TFRJTNL8CY?e=Ij63Td

2) Place the outputs in the same folder as the rest of the code (or adjust paths inside the scripts/notebook accordingly).

3) Run the binning analysis and generate plots:
- Compile/run `binning_program.f90` to produce the `binning_program_L_100_T_*.dat` outputs.
- Run `PLOTS.ipynb` to generate all figures and store outputs in `outputs/`.

4) Optional (Jackknife analysis):
- Keep the `jackknife/` folder structure, place the required `timeseries_*.dat` files where the jackknife programs expect them, and run the provided `*.f90` and `*.gnu` scripts to obtain the heat capacity and susceptibility results.
- Note: `jack_Xi-C_vs_T.out` contains the final \(T, \chi, C\) values (and their errors).


