#!/bin/bash
#SBATCH -J fof
#SBATCH -o /work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/LOG/fof.o1
#SBATCH -e /work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/LOG/fof.e1
#SBATCH -p skx
#SBATCH -N 16
#SBATCH -n 512
#SBATCH -t 4:00:00
#SBATCH --mail-user=mts2188@columbia.edu
#SBATCH --mail-type=NONE
#SBATCH -A TG-AST140041

module load python
module list
pwd
date


LOG_FILE=/work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/LOG/fof
rm $LOG_FILE
for run in 7
do
	ibrun python -u /work2/08288/tg875874/stampede3/enzo_sims/execs/fof.py /scratch/08288/tg875874/enzo_sims/ICs_V0_N512_JLW_0.0/sim$run >> $LOG_FILE
done
