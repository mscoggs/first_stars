#!/bin/bash

cd /work2/08288/tg875874/stampede3/enzo_sims/ICs_V0_N512_JLW_0.0
sbatch --dependency=singleton --job-name=AllJobs1757441357.0885122 /work2/08288/tg875874/stampede3/enzo_sims/sbatch_scripts/enzo.sh


sbatch --dependency=singleton --job-name=AllJobs1757441357.0885122 /work2/08288/tg875874/stampede3/enzo_sims/sbatch_scripts/rockstar.sh


sbatch --dependency=singleton --job-name=AllJobs1757441357.0885122 /work2/08288/tg875874/stampede3/enzo_sims/sbatch_scripts/trees.sh
