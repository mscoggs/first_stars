import os



class ICs_constants(object):
    ICS_SH_FILE = "gen_ICs"
    Z_START = 100
    BOX_SIZE = 0.5
    PARTICLE_DIM = 512
    GRID_DIM = PARTICLE_DIM
    V_BARYONIC = 00
    JOB_NAME = "gen_ICs"
    PARTITION = "icx"
    N_NODES = 1
    N_MPI = 1
    JOB_TIME = "01:00:00"
    BIG_H = 0.67
    OMEGA_M = 0.32
    OMEGA_B = 0.049
    TRANSFER_NAME = "ICs_V"+str(int(V_BARYONIC))+"_N"+str(int(GRID_DIM))
    ICS_DIR = TRANSFER_NAME



def get_IC_keys():
    keys = [name for name,value in vars(ICs_constants).items() if not name.startswith("__")]
    vals = [value for name,value in vars(ICs_constants).items() if not name.startswith("__")]
    return keys,vals


ics = ICs_constants()
keys,vals = get_IC_keys()
fname = ics.ICS_SH_FILE
submit_str = "#!/bin/bash\n\n"


if(not os.path.exists(ics.ICS_DIR)): os.system("mkdir "+ics.ICS_DIR)
if(not os.path.exists(ics.ICS_DIR+"/glass_128_usethis")): os.system("cp /work2/08288/tg875874/stampede3/CICASS/glass/glass_128_usethis "+ics.ICS_DIR+"/glass_128_usethis")

with open("sbatch_scripts/"+fname+"_base.sh", "r") as file:
    text = file.read()
    for key, val in zip(keys,vals):
        text = text.replace(key,str(val))

with open("sbatch_scripts/"+fname+".sh", "w") as file:
    file.write(text)

submit_str += "sbatch sbatch_scripts/"+fname+".sh\n"
with open("submit_"+fname+".sh", "w") as file:
    file.write(submit_str)
os.system("chmod +x submit_"+fname+".sh")





