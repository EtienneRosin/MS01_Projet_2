#!/bin/bash
#SBATCH --job-name=Jacobi_sequential_weak_scalability_measurement
#SBATCH --error=error_%j.txt
#SBATCH --partition=cpu_shared
#SBATCH --account=ams301
#SBATCH --output=weak_%j.txt
#SBATCH --time=00:30:00
#SBATCH --ntasks=40
## load modules
module load cmake
module load gcc
module load gmsh
module load openmpi
## execution

mpirun -display-map ${SLURM_SUBMIT_DIR}/build/3_weak_scalability
done