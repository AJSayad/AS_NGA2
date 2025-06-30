#!/bin/bash

#SBATCH --job-name=2Dshot3                # job name
#SBATCH --output=1kcpd2D.log              # log file
#SBATCH --nodes=8                         # number of nodes 
#SBATCH --ntasks-per-node=64              # number of tasks per node
#SBATCH --time=7-00:00:00                 # max time
#SBATCH --partition=compute-long          # partition
#SBATCH --mail-user=asayad1@stevens.edu   # user email address
#SBATCH --mail-type=ALL                   # type of notifications

export OMPI_MCA_btl_openib_allow_ib=true

# run NGA2
mpiexec -n 1152 ./nga.dp.gnu.opt.mpi.exe -i input -v 2 > 100cpd3D.out
