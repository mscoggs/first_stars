#!/bin/bash
#SBATCH -J trees
#SBATCH -o /work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/LOG/trees.o1
#SBATCH -e /work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/LOG/trees.e1
#SBATCH -p skx
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH --mail-user=mts2188@columbia.edu
#SBATCH --mail-type=NONE
#SBATCH -A TG-AST140041
#SBATCH -A TG-AST140041
module load python
module list
pwd
date
LOG_FILE=/work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/LOG/trees
rm $LOG_FILE 

for run in 7
do
	cd $WORK/enzo_sims
	perl $WORK/rockstar-galaxies/scripts/gen_merger_cfg.pl /scratch/08288/tg875874/enzo_sims/ICs_V0_N512_JLW_0.0/sim$run/rockstar_halos/rockstar.cfg
	cd $WORK/consistent-trees
	make
	perl do_merger_tree.pl /scratch/08288/tg875874/enzo_sims/ICs_V0_N512_JLW_0.0/sim$run/rockstar_halos/outputs/merger_tree.cfg
done
