#!/bin/bash
#$ -N ising_T2.60
#$ -pe smp 1
#$ -q cerqt02.q
#$ -S /bin/bash
#$ -cwd
#$ -o ising_T2.60.out
#$ -e ising_T2.60.err

module load gcc

echo "======================================"
echo "Simulation: L=100, T=2.60"
echo "Started: $(date)"
echo "======================================"

sed -i 's/T = [0-9.]*d0/T = 2.60d0/' main.f90

make clean
make

if [ $? -eq 0 ]; then
    ./mdinamics
    echo "Finished: $(date)"
else
    echo "Compilation failed"
    exit 1
fi
