#!/bin/bash
#SBATCH -J rockstar
#SBATCH -o /work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/LOG/rockstar.o1
#SBATCH -e /work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/LOG/rockstar.e1
#SBATCH -p icx
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 10:00:00
#SBATCH --mail-user=mts2188@columbia.edu
#SBATCH --mail-type=NONE
#SBATCH -A TG-AST140041

module load gcc
module load hdf5
module list
pwd
date

LOG_FILE=/work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/LOG/rockstar
rm $LOG_FILE 
for run in 0 1 2 3 4
do
	DIR=/scratch/08288/tg875874/enzo_sims/ICs_V0_N512_JLW_0.0/sim$run/rockstar_halos
	/work2/08288/tg875874/stampede3/enzo_sims/execs/rockstar-galaxies -c /work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/rockstar$run.cfg  >& $LOG_FILE &

	until [ -e $DIR/auto-rockstar.cfg ]
	do
		sleep 5
	done

	ibrun /work2/08288/tg875874/stampede3/enzo_sims/execs/rockstar-galaxies -c $DIR/auto-rockstar.cfg  >> $LOG_FILE
done
jobs -p | xargs kill
date
