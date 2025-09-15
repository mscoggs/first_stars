#!/bin/bash
#SBATCH -J JOB_NAME      
#SBATCH -o LOG/JOB_NAME.o1
#SBATCH -e LOG/JOB_NAME.e1 
#SBATCH -p PARTITION
#SBATCH -N N_NODES   
#SBATCH -n N_MPI      
#SBATCH -t JOB_TIME    
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
rm LOG/JOB_NAME
./execs/transfer.x -BBOX_SIZE -NPARTICLE_DIM -J1 -VV_BARYONIC -ZZ_START -D1 -SinitSB_transfer_out >> LOG/JOB_NAME
./execs/genICs.x -HBIG_H -OOMEGA_M -BOMEGA_B -LBOX_SIZE -VV_BARYONIC -NPARTICLE_DIM -GGRID_DIM -ZZ_START  -bTRANSFER_NAME -g./ICS_DIR/glass_128_usethis -o./ICS_DIR/ >> LOG/JOB_NAME
date

