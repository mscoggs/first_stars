#!/bin/bash
#SBATCH -J gen_ICs
#SBATCH -o /work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/LOG/gen_ICs.o1
#SBATCH -e /work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/LOG/gen_ICs.e1 
#SBATCH -p icx
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:20:00
#SBATCH --mail-user=mts2188@columbia.edu
#SBATCH --mail-type=NONE
#SBATCH -A TG-AST140041


module load hdf5
module load gsl
module load fftw3
module list
pwd
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/.local/lib
export OMP_NUM_THREADS=4


date
LOG_FILE=/work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/LOG/gen_ICs
rm $LOG_FILE 
/work2/08288/tg875874/stampede3/enzo_sims/execs/transfer.x -B0.5 -N512 -J1 -V0 -Z100 -D1 -SinitSB_transfer_out >> $LOG_FILE
/work2/08288/tg875874/stampede3/enzo_sims/execs/genICs.x -H0.67 -O0.32 -B0.049 -L0.5 -V0 -N512 -G512 -Z100  -bICs_V0_N512 -g/work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/glass_128_usethis -o/work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/ >> $LOG_FILE
date

rm IO_Log
rm initSimCart*
mv Part* ICs_V0_N512_JLW_0.0
mv Grid* ICs_V0_N512_JLW_0.0
mv ICs_V0_N512.pk ICs_V0_N512_JLW_0.0/




