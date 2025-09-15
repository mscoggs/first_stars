import numpy as np
import sys
import glob
import yt
import ytree
import os
import pandas as pd
from astropy.constants import G as G_consts
from setup_run import constants
from cosmo import *
consts = constants()




def is_sf(ds, halo_position,halo_radius,T_virs):
    star_forming= []

    left= (ds.index.grid_left_edge.to("Mpccm/h")).value
    right=(ds.index.grid_right_edge.to("Mpccm/h")).value
    levels= ds.index.grid_levels.flatten()

    for N in range(np.size(T_virs)):
        r = halo_radius[N]*2
        args = np.where(np.all(halo_position[N] >= left-r,axis=1) & np.all(halo_position[N] <= right+r, axis=1))
        halo_levels= levels[args]
        max_levels_args = np.where(halo_levels == np.max(halo_levels))[0]
        target_grids = ds.index.grids[args][max_levels_args]

        sf = 0
        for grid in target_grids:
            H2_fraction = grid["gas", "H2_p0_fraction"].value
            n_density = grid["gas", "number_density"].value
            temp = grid["gas", "temperature"].value
            density_args = np.where(n_density>=100)
            if(np.size(density_args) == 0): continue

            temp = temp[density_args]
            H2_fraction = H2_fraction[density_args]
            if(np.any(temp<=T_virs[N]*0.5)): sf = 1    #kulkarni conditions
        star_forming.append(sf)


    return np.array(star_forming)



def get_mcrit_kulkarni(x,y):

    bm,f,bs= iterate_bin_masses(x,y)
    if(np.size(bm) == 0): return x[np.where(y==1)[0][0]], -1,-1
    args = np.where((f<=consts.FMAX) & (f>=consts.FMIN))[0]
    bm_filt, f_filt = bm[args], f[args]
    mcrit  = np.interp(0.5, f_filt, bm_filt)
    err_up,err_dn = mcrit-bs/2.0, mcrit+bs/2.0
    mcrit, err_up, err_dn = 10**mcrit, 10**err_up, 10**err_dn
    return mcrit,err_up,err_dn




def iterate_bin_masses(x,y):
    x,y = np.log10(x),np.array(y)
    for steps in range(consts.MAX_STEPS, consts.MIN_STEPS-1, -1):
        bs = (np.max(x)-np.min(x))/steps
        bm,f = bin_mass(x,y,bs=bs)
        if(check_continuous_and_mono(bm,f)):
            return bm, f,bs
    return [], [], []

def check_continuous_and_mono(bm,f):
    args = np.where((f<=consts.FMAX) & (f>=consts.FMIN))[0]
    if(np.size(args)<2):return False
    cont = np.array(range(args[0],args[-1]+1,1))
    if(np.size(args) != np.size(cont)): return False
    if(not np.all(cont==args)): return False


    f_filt = f[args]
    f_dif = f_filt[1:]-f_filt[:-1]
    if(np.any(f_dif <= 0)): return False

    if(np.any(f[args[-1]+1:] <consts.FMAX )): return False
    if(np.any(f[:args[0]] >consts.FMIN )): return False
    return True


def bin_mass(x,y,bs=0.1):
    min_mass = np.min(x)
    max_mass = np.max(x)
    l,r = min_mass, min_mass+bs
    bm,f = [],[]
    while(l < max_mass):
        args = np.where((x>=l) & (x<r))
        cm = (l+r)/2.0
        l+=bs
        r+=bs
        if(np.size(args)==0): continue
        bm.append(cm)
        f.append(np.mean(y[args]))
    return np.array(bm), np.array(f)


def get_tree_data(t,min_mass):
    halo_mass = t["tree", "Mvir"].to("Msun").value
    halo_args = np.where(halo_mass >= min_mass)
    halo_mass = halo_mass[halo_args]
    redshift = (1/t["tree", "scale"]-1)[halo_args]
    halo_position = t["tree", "position"][halo_args].to("Mpc/h").value
    x,y,z = halo_position.T
    halo_radius = t["tree", "Rvir"][halo_args].to("Mpc/h").value
    T_virs = calc_T_vir(halo_mass,redshift)
    T_ids = t["tree", "Tree_root_ID"][halo_args]
    return halo_mass, redshift,x,y,z, halo_radius, T_virs, T_ids






def main():
    '''
    args = sys.argv[1:]
    if(len(args) != 5):
        print("Error, improper number of arguments")
        sys.exit(0)
    enzo_dir,halo_dir,out_dir, min_mass, save_tail = args

    min_mass = float(min_mass)
    sim_num = enzo_dir.split("sim")[-1]
    halo_type = (halo_dir.split("/")[-1]).split("_")[0]
    snap_nums = [snapstr[-4:] for snapstr in list(glob.glob(enzo_dir + "/RD*"))][::-1]

    names =['z','n_halos', 'mcrit_k', 'mcrit_euk', 'mcrit_edk', "mmin_k",  'mave_k', "mave_std_k", "n_sf_kulkani", 'mcrit_s', 'mcrit_eus', 'mcrit_eds',"mmin_s", 'mave_s', "mave_std_s", "n_sf_s"]
    df = pd.DataFrame(columns=names)
    df_row = 0
    
    '''

    enzo_dir = "/scratch/08288/tg875874/enzo_sims/ICs_V0_N512_JLW_1.0/sim0/"
    snap_nums = ["0010", "0011", "0012", "0013"][::-1]
    min_mass = 5e4






    print("loading trees")
    a = ytree.load("/scratch/08288/tg875874/enzo_sims/ICs_V0_N512_JLW_1.0/sim0/rockstar_halos/trees/tree_0_0_0.dat")

    halo_mass = a["Mvir"].to("Msun").value
    halo_args = np.where(halo_mass >= min_mass)
    halo_mass = halo_mass[halo_args]
    redshift = (1/a["scale"]-1)[halo_args]
    halo_position = a["position"][halo_args].to("Mpc/h").value
    halo_radius = a["Rvir"][halo_args].to("Mpc/h").value
    T_virs = calc_T_vir(halo_mass,redshift)

    print("\n loading enzo")
    ds  = yt.load(enzo_dir + "/RD" + snap_nums[0] + "/RedshiftOutput"+snap_nums[0])

    print("checking for SF")
    #sf = is_sf(ds, halo_position,halo_radius,T_virs)
    #sf_args = np.where(sf == 1)
    sf_args = np.array([0,1,2,3,5,7,8])
    sf_trees = np.array(list(a))[halo_args][sf_args]
    
    rows = [[],[],[],[],[],[],[],[]]
    for t in sf_trees:
        td = get_tree_data(t,min_mass)
        for i,info in enumerate(td):
            rows[i].extend(info)

    rows = np.array([np.array(x) for x in rows]).T
    names =["mass", "redshift", "x","y","z", "radius", "temp", "t_id"]
    df = pd.DataFrame(rows,columns=names)
    df["is_sf"] = 0
    print(df)

    for enzo_snap in snap_nums[1:]:
        print("iterating through trees across redshift")
        ds  = yt.load(enzo_dir + "/RD" + enzo_snap + "/RedshiftOutput"+enzo_snap)
        units = ds.units
        target_redshift = ds.current_redshift

        mask_z = (df["redshift"]-target_redshift).abs() < 0.1
        print(mask_z)
        df_target = df[mask_z]
        if(np.size(df_target) == 0):
            print("no halos for redshift", target_redshift)
            continue

        halo_position = np.array([df_target["x"], df_target["y"], df_target["z"]]).T
        halo_radius = np.array(df_target["radius"])
        T_virs = np.array(df_target["temp"])
        sf = is_sf(ds, halo_position,halo_radius,T_virs)
        df.loc[mask_z, "is_sf"] = sf

    print(df)
    sys.exit(0)







'''
        star_forming_k, star_forming_s  = apply_sf_condition(ds, halo_position,halo_radius,T_virs)
        if(np.size(np.where(star_forming_k == 1)) == 0):
            print("No star forming halos for z ",z)
            continue
        mcrit1,erru1,errd1 = get_mcrit_kulkarni(halo_mass,star_forming_k)
        mcrit2,erru2,errd2 = get_mcrit_kulkarni(halo_mass,star_forming_s)

        
        m_sf1 = halo_mass[np.where(star_forming_k == 1)]
        m_sf2 = halo_mass[np.where(star_forming_s == 1)]
        
        data = [z, np.size(halo_mass), mcrit1, erru1,errd1, np.min(m_sf1), np.mean(m_sf1), np.std(m_sf1), np.size(m_sf1),  mcrit2, erru2,errd2, np.min(m_sf2), np.mean(m_sf2), np.std(m_sf2), np.size(m_sf2)]
        df.loc[df_row] = np.array(data)
        df_row += 1

    df.to_csv(out_dir+"/mcrit_"+sim_num+"_"+halo_type+"_"+save_tail+".txt")
'''


if __name__=="__main__":
    main()
