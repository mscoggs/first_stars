#!/bin/bash
#SBATCH -J JOB_NAME_FOF
#SBATCH -o WORK_DIR/ICS_DIR/LOG/JOB_NAME_FOF.o1
#SBATCH -e WORK_DIR/ICS_DIR/LOG/JOB_NAME_FOF.e1
#SBATCH -p PARTITION_FOF
#SBATCH -N N_NODES_FOF
#SBATCH -n N_MPI_FOF
#SBATCH -t JOB_TIME_FOF
#SBATCH --mail-user=mts2188@columbia.edu
#SBATCH --mail-type=NONE
#SBATCH -A TG-AST140041

module load python
module list
pwd
date


LOG_FILE=WORK_DIR/ICS_DIR/LOG/JOB_NAME_FOF
rm $LOG_FILE
for run in RUN_NUMS_BASH
do
	ibrun python -u WORK_DIR/execs/fof.py OUT_DIR/sim$run >> $LOG_FILE
done
