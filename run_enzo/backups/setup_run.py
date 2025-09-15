import os
import sys
import time
import copy
import glob
import numpy as np
from datetime import timedelta




class constants(object):
    #GENERAL
    BIG_H = 0.67
    OMEGA_M = 0.32
    OMEGA_B = 0.049
    OMEGA_CDM = OMEGA_M - OMEGA_B
    OMEGA_L = 1.0 - OMEGA_M
    OMEGA_K = 1.0-OMEGA_M-OMEGA_L
    MU_CONSTANT = 1.22 
    Z_START = 100
    Z_END = 15.0
    BOX_SIZE_1D = 0.5
    HALO_FOLDERS = " ".join(["rockstar_halos", "fof_halos"])


    '''
    PARTICLE_DIM = 1024
    RUN_NUMBERS = [0]
    REFINEMENT_LEVELS = [3] 
    REFINE_BYS = [2]
    JEANS_LENGTHS_CRITERIA = [True]
    '''

    '''
    PARTICLE_DIM = 512
    RUN_NUMBERS = [0,1,2,3,4]
    REFINEMENT_LEVELS = [3,4,5,6,5] 
    REFINE_BYS = [2,2,2,2,2]
    JEANS_LENGTHS_CRITERIA = [True,True,True,True,False]
    '''
    #BIG BATCH FOR NEW PAPER
    PARTICLE_DIM = 512
    RUN_NUMBERS = [0]
    REFINEMENT_LEVELS = [5]
    REFINE_BYS = [2]
    JEANS_LENGTHS_CRITERIA = [True]
    J_LW_VALUE = 1.0



    '''
    JEANS_LENGTHS_CRITERIA = [True]
    PARTICLE_DIM = 1024
    RUN_NUMBERS = [0]
    REFINEMENT_LEVELS = [3] 
    '''

    GRID_DIM = PARTICLE_DIM
    V_BARYONIC = 00
    #J_LW_VALUE = 0.0
    TRANSFER_NAME = "ICs_V"+str(int(V_BARYONIC))+"_N"+str(int(GRID_DIM))

    #DIRS
    WORK_DIR = "/work2/08288/tg875874/stampede3/enzo_sims"
    ICS_DIR = TRANSFER_NAME+"_JLW_"+str(J_LW_VALUE)
    OUT_DIR_BASE = "/scratch/08288/tg875874/enzo_sims"
    OUT_DIR = OUT_DIR_BASE + "/" + ICS_DIR

    #ICS
    JOB_NAME_IC = "gen_ICs"
    PARTITION_IC = "skx"
    N_NODES_IC = 1
    N_MPI_IC = 1
    JOB_TIME_IC_MIN = 7*(PARTICLE_DIM/512)**3
    JOB_TIME_IC = str(timedelta(minutes=round(JOB_TIME_IC_MIN+5)))
    #RING
    JOB_NAME_RING = "ring"
    PARTITION_RING = "skx"
    N_NODES_RING = int(PARTICLE_DIM/32)
    N_MPI_RING = PARTICLE_DIM
    JOB_TIME_RING_MIN = 2*(PARTICLE_DIM/512)**3
    JOB_TIME_RING = str(timedelta(minutes=round(JOB_TIME_RING_MIN+5)))


    #ENZO
    JOB_NAME_ENZO = "enzo"
    PARTITION_ENZO = "skx"
    #N_NODES_ENZO = int(PARTICLE_DIM/64)
    N_NODES_ENZO = int(PARTICLE_DIM/32)
    N_MPI_ENZO = PARTICLE_DIM
    JOB_TIME_ENZO_HR = 5.5*(PARTICLE_DIM/N_MPI_ENZO)**3
    JOB_TIME_ENZO = str(timedelta(hours=JOB_TIME_ENZO_HR*np.size(RUN_NUMBERS)+1))
    COSMOLOGY_OUTPUT_ZS = [30,27,25,24,23,22,21,20,19,18,17,16,15]
    COSMOLOGY_OUTPUT_LIST = ''.join(["CosmologyOutputRedshift["+str(num+1)+"] = "+str(z)+"\n" for num, z in enumerate(COSMOLOGY_OUTPUT_ZS)])


    #ROCKSTAR
    JOB_NAME_RS = "rockstar"
    SNAPSHOT_FILE_NAME = WORK_DIR+"/"+ICS_DIR+"/snapshots.dat"
    PARTITION_RS = "skx"
    N_NODES_RS = 1
    N_MPI_RS = 16
    JOB_TIME_RS_MIN = 25*((16/N_MPI_RS)*(PARTICLE_DIM/512))**3
    JOB_TIME_RS = str(timedelta(minutes=JOB_TIME_RS_MIN*np.size(RUN_NUMBERS)+30))
    RUN_NUMS_BASH = " ".join([str(RUN_IT) for RUN_IT in RUN_NUMBERS])

    #FOF
    JOB_NAME_FOF = "fof"
    PARTITION_FOF = "skx"
    N_NODES_FOF = int(PARTICLE_DIM/32)
    N_MPI_FOF = PARTICLE_DIM
    JOB_TIME_FOF_HR = 3.5
    JOB_TIME_FOF = str(timedelta(hours=JOB_TIME_FOF_HR*np.size(RUN_NUMBERS)+0.5))


    #MCRIT
    JOB_NAME_MCRIT = "mcrit"
    PARTITION_MCRIT = "skx"
    N_NODES_MCRIT = 1
    N_MPI_MCRIT = 1
    JOB_TIME_MCRIT_MIN = 30.0
    JOB_TIME_MCRIT = str(timedelta(minutes=JOB_TIME_MCRIT_MIN*np.size(RUN_NUMBERS)+30))
    MIN_MASS_MCRIT = 5e4
    SAVE_TAIL_MCRIT = "Round1"
    FMAX, FMIN = 0.75, 0.25
    MAX_STEPS, MIN_STEPS = 300,3



    #PARAMETERS_IN = sys.argv[-1]

    

def get_keys():
    keys = [name for name,value in vars(constants).items() if not name.startswith("__")]
    vals = [value for name,value in vars(constants).items() if not name.startswith("__")]
    return keys,vals

def write_sbatch(fname,work_dir,submit_dir =-1):
    submit_str = "#!/bin/bash\n\n"
    if(submit_dir != -1): submit_str += "cd "+submit_dir+"\n"
    submit_str += "sbatch "+work_dir+"/sbatch_scripts/"+fname+".sh\n"
    count=str(np.size(list(glob.iglob("submit*.sh"))))
    fname = count+"_"+fname
    with open("submit"+fname+".sh", "w") as file:
        file.write(submit_str)
    os.system("chmod +x submit"+fname+".sh")

def replace_text(fname,infile,outfile):
    consts = constants()
    keys,vals = get_keys()
    with open(infile, "r") as file:
        text = file.read()
        for key, val in zip(keys,vals):
            text = text.replace(key,str(val))
    
    if("param" in infile or "rockstar_base.cfg" in infile):
        for x in range(len(consts.RUN_NUMBERS)):
            temp_text = text
            outfile_temp = outfile+str(consts.RUN_NUMBERS[x])+"."+infile.split(".")[-1]
            temp_text = temp_text.replace("REFINEMENT_LEVEL", str(consts.REFINEMENT_LEVELS[x]))
            temp_text = temp_text.replace("REFINE_BY", str(consts.REFINE_BYS[x]))
            temp_text = temp_text.replace("RUN_NUMBER", str(consts.RUN_NUMBERS[x]))
            if(consts.JEANS_LENGTHS_CRITERIA[x]):
                temp_text = temp_text.replace("JEANS_LENGTH_ON", "6")
                temp_text = temp_text.replace("JEANS_LENGTH_VAL", "0")
            else:
                temp_text = temp_text.replace("JEANS_LENGTH_ON", "")
                temp_text = temp_text.replace("JEANS_LENGTH_VAL", "")
            with open(outfile_temp, "w") as file:
                file.write(temp_text)
    else:
        with open(outfile, "w") as file:
            file.write(text)

def make_snapshots_file():
    consts = constants()
    num_snaps = np.size(consts.COSMOLOGY_OUTPUT_ZS)
    fname = consts.SNAPSHOT_FILE_NAME
    text = "".join(["RD00"+"{:02d}".format(x+1)+"/RedshiftOutput00"+"{:02d}".format(x+1)+"\n" for x in range(num_snaps)])
    with open(fname, "w") as file:
        file.write(text)


def combine_all_jobs():
    consts = constants()
    all_files = np.sort(list(glob.iglob("submit*.sh")))
    text = ""
    for file0 in all_files:
        with open(file0, "r") as file:
            text += file.read()
    text = text.replace("#!/bin/bash",  "")
    text = text.replace("sbatch ", "sbatch --dependency=singleton --job-name=AllJobs"+str(time.time())+ " ")
    text = "#!/bin/bash\n"+text
    with open("submit_all_jobs.sh", "w") as file:
        file.write(text)
    os.system("chmod +x submit_all_jobs.sh")

def make_dir(dir_):
    if(not os.path.exists(dir_)): os.system("mkdir "+dir_)


def main():
    consts = constants()

    print("\nCREATING DIRECTORIES:")
    os.system("rm submit*")
    dirs_list = [consts.ICS_DIR, consts.OUT_DIR_BASE, consts.OUT_DIR, consts.ICS_DIR+"/mcrit_files"]
    for dir_ in dirs_list: make_dir(dir_)

    print("\nINITIAL CONDITIONS:")
    fname = consts.JOB_NAME_IC
    print("     Copying glass file")
    if(not os.path.exists(consts.ICS_DIR+"/glass_128_usethis")): os.system("cp /work2/08288/tg875874/stampede3/CICASS/glass/glass_128_usethis "+consts.ICS_DIR+"/glass_128_usethis")
    print("     Writing ICs batch file")
    replace_text(fname,"sbatch_scripts/"+fname+"_base.sh","sbatch_scripts/"+fname+".sh")
    write_sbatch(fname, consts.WORK_DIR)

    print("\nRING:")
    fname = consts.JOB_NAME_RING
    print("     Writing ring batch file")
    replace_text(fname,"sbatch_scripts/"+fname+"_base.sh","sbatch_scripts/"+fname+".sh")
    write_sbatch(fname,consts.WORK_DIR, submit_dir=consts.WORK_DIR+"/"+consts.ICS_DIR)

    print("\nENZO:")
    fname = consts.JOB_NAME_ENZO
    print("     Setting up LW_J21")
    os.system("cp inputs/* "+consts.ICS_DIR)
    replace_text(fname,"inputs/LW_J21.in", consts.ICS_DIR+"/LW_J21.in")
    print("     Setting up enzo parameter files")
    replace_text(fname,"inputs/params_base.txt", consts.ICS_DIR+"/params")
    print("     Creating enzo out dirs")
    for run in consts.RUN_NUMBERS: make_dir(consts.OUT_DIR+"/sim"+str(run))
    print("     Writing enzo batch file")
    replace_text(fname,"sbatch_scripts/"+fname+"_base.sh","sbatch_scripts/"+fname+".sh")
    write_sbatch(fname,consts.WORK_DIR, submit_dir=consts.WORK_DIR+"/"+consts.ICS_DIR)

    if("rockstar_halos" in consts.HALO_FOLDERS):
        print("\nROCKSTAR:")
        fname = consts.JOB_NAME_RS
        print("     Making pfs.dat file")
        for run in consts.RUN_NUMBERS: make_dir(consts.OUT_DIR+"/sim"+str(run)+"/rockstar_halos")
        make_snapshots_file()
        print("     Making rockstar.cfg files")
        replace_text(fname,"inputs/rockstar_base.cfg", consts.ICS_DIR+"/rockstar")
        print("     Writing rockstar batch file")
        replace_text(fname,"sbatch_scripts/"+fname+"_base.sh","sbatch_scripts/"+fname+".sh")
        write_sbatch(fname,consts.WORK_DIR)

    if("fof_halos" in consts.HALO_FOLDERS):
        print("\nFOF:")
        fname = consts.JOB_NAME_FOF
        for run in consts.RUN_NUMBERS: make_dir(consts.OUT_DIR+"/sim"+str(run)+"/fof_halos")
        print("     Writing fof batch file")
        replace_text(fname,"sbatch_scripts/"+fname+"_base.sh","sbatch_scripts/"+fname+".sh")
        write_sbatch(fname,consts.WORK_DIR)

    print("\nMCRIT:")
    fname = consts.JOB_NAME_MCRIT
    print("     Writing mcrit batch file")
    replace_text(fname,"sbatch_scripts/"+fname+"_base.sh","sbatch_scripts/"+fname+".sh")
    write_sbatch(fname,consts.WORK_DIR)

    print("\n")
    combine_all_jobs()
    make_dir(consts.WORK_DIR+"/"+consts.ICS_DIR+"/sbatch_saves")
    for submit_file in list(glob.iglob("submit*.sh")): os.system("cp "+submit_file+" "+consts.ICS_DIR+"/sbatch_saves/")
    os.system("cp sbatch_scripts/* "+consts.ICS_DIR+"/sbatch_saves/")


if __name__=="__main__":
    main()
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

	perl $WORK/rockstar-galaxies/scripts/gen_merger_cfg.pl OUT_DIR/sim$run/rockstar_halos/rockstar.cfg
done
date
