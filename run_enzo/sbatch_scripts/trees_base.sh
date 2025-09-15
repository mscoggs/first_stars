#!/bin/bash
#SBATCH -J JOB_NAME_TREES
#SBATCH -o WORK_DIR/ICS_DIR/LOG/JOB_NAME_TREES.o1
#SBATCH -e WORK_DIR/ICS_DIR/LOG/JOB_NAME_TREES.e1
#SBATCH -p PARTITION_TREES
#SBATCH -N N_NODES_TREES
#SBATCH -n N_MPI_TREES
#SBATCH -t JOB_TIME_TREES
#SBATCH --mail-user=mts2188@columbia.edu
#SBATCH --mail-type=NONE
#SBATCH -A TG-AST140041
#SBATCH -A TG-AST140041
module load python
module list
pwd
date
LOG_FILE=WORK_DIR/ICS_DIR/LOG/JOB_NAME_TREES
rm $LOG_FILE 

for run in RUN_NUMS_BASH
do
	cd $WORK/enzo_sims
	perl $WORK/rockstar-galaxies/scripts/gen_merger_cfg.pl OUT_DIR/sim$run/rockstar_halos/rockstar.cfg
	cd $WORK/consistent-trees
	make
	perl do_merger_tree.pl OUT_DIR/sim$run/rockstar_halos/outputs/merger_tree.cfg
done
