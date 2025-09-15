#!/bin/bash
#SBATCH -J JOB_NAME_RS
#SBATCH -o WORK_DIR/ICS_DIR/LOG/JOB_NAME_RS.o1
#SBATCH -e WORK_DIR/ICS_DIR/LOG/JOB_NAME_RS.e1
#SBATCH -p PARTITION_RS
#SBATCH -N N_NODES_RS
#SBATCH -n N_MPI_RS
#SBATCH -t JOB_TIME_RS
#SBATCH --mail-user=mts2188@columbia.edu
#SBATCH --mail-type=NONE
#SBATCH -A TG-AST140041

module load gcc
module load hdf5
module list
pwd
date

LOG_FILE=WORK_DIR/ICS_DIR/LOG/JOB_NAME_RS
rm $LOG_FILE 
for run in RUN_NUMS_BASH
do
	DIR=OUT_DIR/sim$run/rockstar_halos
	WORK_DIR/execs/rockstar-galaxies -c WORK_DIR/ICS_DIR/rockstar$run.cfg  >& $LOG_FILE &

	until [ -e $DIR/auto-rockstar.cfg ]
	do
		sleep 5
	done

	ibrun WORK_DIR/execs/rockstar-galaxies -c $DIR/auto-rockstar.cfg  >> $LOG_FILE
done
jobs -p | xargs kill
date
