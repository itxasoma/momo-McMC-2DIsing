#!/bin/bash
#SBATCH -J ising_TAll
#SBATCH -o ising_TAll.out
#SBATCH -e ising_TAll.err
#SBATCH --mail-user=imunozal74@alumnes.ub.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -p std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04-00:00:00

set -e

TMIN=2.0
TMAX=3.0
STEP=0.1

TEMPLATE="DATA.in"

make clean && make

mkdir -p runs

pids=()

# Do 2.0 to 3.0 in 0.1 time steps 
for ti in {20..30}; do
  T=$(echo "scale=1; $ti/10" | bc)                       # 2.0, 2.1, ..., 3.0
  TFILE=$(printf "runs/DATA_T%0.1f.in" "$T")

  # Write a new input as DATA.in, but changing the 2nd line
  awk -v newT="${T}d0" 'NR==2{$0=newT} {print}' "$TEMPLATE" > "$TFILE"

  echo "Launching T=$T  ->  $TFILE"
  ./mdinamics "$TFILE" > "runs/out_T$(printf "%0.1f" "$T").log" 2>&1 &
  pids+=($!)
done

# Wait for all to finish
for pid in "${pids[@]}"; do
  wait "$pid"
done

echo "All runs finished."

