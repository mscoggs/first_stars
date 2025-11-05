#!/bin/bash
#SBATCH -J ring
#SBATCH -o /work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/LOG/ring.o1
#SBATCH -e /work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/LOG/ring.e1 
#SBATCH -p skx
#SBATCH -N 16
#SBATCH -n 512
#SBATCH -t 0:07:00
#SBATCH --mail-user=mts2188@columbia.edu
#SBATCH --mail-type=NONE
#SBATCH -A TG-AST140041


module load hdf5
module list
pwd
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/.local/lib

date
LOG_FILE=/work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0/LOG/ring
rm $LOG_FILE
mpirun -n 512 /work2/08288/tg875874/stampede3/enzo_sims/execs/ring.exe pv ParticlePositions ParticleVelocities >> $LOG_FILE
date
