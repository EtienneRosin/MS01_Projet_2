#!/bin/bash
#SBATCH --job-name=Project_2_scalability_measurement
#SBATCH --error=error_%j.txt
#SBATCH --partition=cpu_shared
#SBATCH --account=ams301
#SBATCH --output=strong_%j.txt
#SBATCH --time=00:30:00
#SBATCH --ntasks=40
## load modules
module load cmake
module load gcc
module load gmsh
module load openmpi
## execution
gmsh -2 -part 40 ${SLURM_SUBMIT_DIR}/benchmark/mesh.geo
mpirun -display-map ${SLURM_SUBMIT_DIR}/solver
done