#!/bin/bash
#SBATCH -J JOB_NAME_ENZO
#SBATCH -o WORK_DIR/ICS_DIR/LOG/JOB_NAME_ENZO.o1
#SBATCH -e WORK_DIR/ICS_DIR/LOG/JOB_NAME_ENZO.e1
#SBATCH -p PARTITION_ENZO
#SBATCH -N N_NODES_ENZO
#SBATCH -n N_MPI_ENZO
#SBATCH -t JOB_TIME_ENZO
###SBTACH --exclusive
#SBATCH --mail-user=mts2188@columbia.edu
#SBATCH --mail-type=NONE
#SBATCH -A TG-AST140041

module load hdf5
module list
pwd
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/.local/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/.local
date

LOG_FILE=WORK_DIR/ICS_DIR/LOG/JOB_NAME_ENZO
rm $LOG_FILE 
for run in RUN_NUMS_BASH
do
	ibrun WORK_DIR/execs/enzo.exe -d params$run.txt >> $LOG_FILE
	#mpirun -n N_MPI_ENZO WORK_DIR/execs/enzo.exe -d params$run.txt >> $LOG_FILE

done
date
