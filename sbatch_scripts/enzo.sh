#!/bin/bash
#SBATCH -J enzo
#SBATCH -o /work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/LOG/enzo.o1
#SBATCH -e /work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/LOG/enzo.e1
#SBATCH -p skx
#SBATCH -N 16
#SBATCH -n 512
#SBATCH -t 12:00:00
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

LOG_FILE=/work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/LOG/enzo
rm $LOG_FILE 
for run in 7
do
	ibrun /work2/08288/tg875874/stampede3/enzo_sims/execs/enzo.exe -d params$run.txt >> $LOG_FILE
	#mpirun -n 512 /work2/08288/tg875874/stampede3/enzo_sims/execs/enzo.exe -d params$run.txt >> $LOG_FILE

done
date
