#!/bin/bash
#SBATCH -J mcrit
#SBATCH -o /work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/LOG/mcrit.o1
#SBATCH -e /work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/LOG/mcrit.e1
#SBATCH -p icx
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH --mail-user=mts2188@columbia.edu
#SBATCH --mail-type=NONE
#SBATCH -A TG-AST140041
#SBATCH -A TG-AST140041
module load python
module list
pwd
date
LOG_FILE=/work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/LOG/mcrit
rm $LOG_FILE 

for run in 0 1 2 3 4
do
	for halo_type in rockstar_halos fof_halos
	do
		python -u /work2/08288/tg875874/stampede3/enzo_sims/execs/mcrit.py /scratch/08288/tg875874/enzo_sims/ICs_V0_N512_JLW_0.0/sim$run/$halo_type /work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/mcrit_files 50000.0 Round1 >> $LOG_FILE
	done
done
date
