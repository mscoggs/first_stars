#!/bin/bash
#SBATCH -J JOB_NAME_RING
#SBATCH -o WORK_DIR/ICS_DIR/LOG/JOB_NAME_RING.o1
#SBATCH -e WORK_DIR/ICS_DIR/LOG/JOB_NAME_RING.e1 
#SBATCH -p PARTITION_RING
#SBATCH -N N_NODES_RING
#SBATCH -n N_MPI_RING
#SBATCH -t JOB_TIME_RING
#SBATCH --mail-user=mts2188@columbia.edu
#SBATCH --mail-type=NONE
#SBATCH -A TG-AST140041


module load hdf5
module list
pwd
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/.local/lib

date
LOG_FILE=WORK_DIR/ICS_DIR/LOG/JOB_NAME_RING
rm $LOG_FILE
mpirun -n N_MPI_RING WORK_DIR/execs/ring.exe pv ParticlePositions ParticleVelocities >> $LOG_FILE
date
