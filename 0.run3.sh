#!/bin/bash
#SBATCH -J ising_Tmore
#SBATCH -o ising_Tmore.out
#SBATCH -e ising_Tmore.err
#SBATCH --mail-user=imunozal74@alumnes.ub.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -p std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04-00:00:00

set -e

TEMPLATE="DATA.in"

make clean && make

mkdir -p runs

pids=()

# ----------------------------
# Sweep: T = 2.0..3.0 step 0.1
# ----------------------------
for ti in {20..30}; do
  T=$(echo "scale=1; $ti/10" | bc)
  TFILE=$(printf "runs/DATA_T%0.1f.in" "$T")

  awk -v newT="${T}d0" 'NR==2{$0=newT} {print}' "$TEMPLATE" > "$TFILE"

  echo "Launching sweep T=$T  ->  $TFILE"
  ./mdinamics "$TFILE" > "runs/out_T$(printf "%0.1f" "$T").log" 2>&1 &
  pids+=($!)
done

# ----------------------------
# Extra sim #1: L=100, T=2.27, nMCS=1e8
# ----------------------------
cat > runs/DATA_L100_T2p27_MC1e8.in <<'EOF'
100
2.27d0
100000000
10
1000
12345
4
EOF
echo "Launching extra: runs/DATA_L100_T2p27_MC1e8.in"
./mdinamics "runs/DATA_L100_T2p27_MC1e8.in" > "runs/out_L100_T2p27_MC1e8.log" 2>&1 &
pids+=($!)

# ----------------------------
# Extra sim #2: L=20, T=2.00, nMCS=1e6
# ----------------------------
cat > runs/DATA_L20_T2p00_MC1e6.in <<'EOF'
20
2.00d0
1000000
10
1000
12345
4
EOF
echo "Launching extra: runs/DATA_L20_T2p00_MC1e6.in"
./mdinamics "runs/DATA_L20_T2p00_MC1e6.in" > "runs/out_L20_T2p00_MC1e6.log" 2>&1 &
pids+=($!)

# Wait for all to finish (so the sbatch job doesn't exit early) [web:88]
for pid in "${pids[@]}"; do
  wait "$pid"
done

echo "All runs finished."


