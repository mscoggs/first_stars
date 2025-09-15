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
    print("checking for SF in", np.size(T_virs), "trees")

    for N in range(np.size(T_virs)):
        if(N%50 ==0): print("on", N," of", np.size(T_virs))
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
    ids = t["tree", "id"][halo_args]
    d_ids = t["tree", "desc_id"][halo_args]
    return halo_mass, redshift,x,y,z, halo_radius, T_virs, T_ids, ids, d_ids


def dynamical_heating(mass_list,dm_list, dt_list,z_list, temp_list):
    dmdt = dm_list/dt_list
    dedt =  temp_list/mass_list * K_CONST/(5.0/3.0-1.0) * dmdt/SEC_PER_MYR
    return dedt




def main():
    '''
    args = sys.argv[1:]
    if(len(args) != 5):
        print("Error, improper number of arguments")
        sys.exit(0)
    enzo_dir,out_dir, min_mass, save_tail = args

    min_mass = float(min_mass)
    '''

    '''

    enzo_dir = "/scratch/08288/tg875874/enzo_sims/ICs_V0_N512_JLW_1.0/sim0/"
    sim_num = enzo_dir.split("sim")[-1]

    snap_nums = [snapstr[-4:] for snapstr in list(glob.glob(enzo_dir + "/RD*"))][::-1]
    min_mass = 3e4


    print("snap nums:", snap_nums)






    print("loading trees")
    tree_file = enzo_dir + "/rockstar_halos/trees/tree_0_0_0.dat"
    a = ytree.load(tree_file)

    halo_mass = a["Mvir"].to("Msun").value
    halo_args = np.where(halo_mass >= min_mass)
    print(np.size(halo_args))

    print("loading trees")
    #candidate_trees = [np.array(list(a))[0]]
    candidate_trees = np.array(list(a))[halo_args]
    rows = [[],[],[],[],[],[],[],[],[],[]]
    for t in candidate_trees:
        td = get_tree_data(t,min_mass)
        for i,info in enumerate(td):
            rows[i].extend(info)

    rows = np.array([np.array(x) for x in rows]).T
    names =["mass", "redshift", "x","y","z", "radius", "temp", "t_id", "id", "d_id"]
    df = pd.DataFrame(rows,columns=names)
    df["is_sf"] = 0

    for enzo_snap in snap_nums:
        ds  = yt.load(enzo_dir + "/RD" + enzo_snap + "/RedshiftOutput"+enzo_snap)
        units = ds.units
        target_redshift = ds.current_redshift

        mask_z = (df["redshift"]-target_redshift).abs() < 0.1
        df_target = df[mask_z]
        if(np.size(df_target) == 0):
            print("no halos for redshift", target_redshift)
            continue

        halo_position = np.array([df_target["x"], df_target["y"], df_target["z"]]).T
        halo_radius = np.array(df_target["radius"])
        T_virs = np.array(df_target["temp"])
        sf = is_sf(ds, halo_position,halo_radius,T_virs)
        print(np.size(np.where(sf == 1.0)), "star forming")
        df.loc[mask_z, "is_sf"] = sf

    df.to_csv("test.txt")
    #df = pd.read_csv("test.txt", index_col = 0)
    df["is_sf_first"] = 0
    sf_already_marked = []
    print(df)

    for tree_id, tree_df in df.groupby("t_id"):
        leaves = set(tree_df["id"]) - set(tree_df["d_id"][tree_df["d_id"] != -1])
        for leaf_id in leaves:
            current_id = leaf_id
            branch_halos = []
            while current_id != -1:
                if(current_id not in np.array(tree_df["id"])):
                        print("skipping", current_id)
                        break
                halo = tree_df.loc[tree_df["id"] == current_id].iloc[0]
                branch_halos.append(halo["id"])
                current_id = halo["d_id"]
            for hid in branch_halos:
                if df.loc[df["id"] == hid, "is_sf"].values[0] == 1 and hid not in sf_already_marked:
                    df.loc[df["id"] == hid, "is_sf_first"] = 1
                    sf_already_marked.extend(branch_halos)

            




    print(df)
    print(np.size(np.where(df["is_sf_first"] == 1)))
    df.to_csv("test2.txt")
    '''

    df = pd.read_csv("test2.txt", index_col = 0)

    #getting the mm_prog_ids
    idx = df.loc[df["d_id"] != -1].groupby("d_id")["mass"].idxmax()
    mm_map = df.loc[idx].set_index("d_id")["id"]
    df["mm_prog_id"] = df["id"].map(mm_map).fillna(-1)
    print(df)

    #calculating dm
    dm = df["mass"] - df["mm_prog_id"].replace(-1,pd.NA).map(df.set_index("id")["mass"])
    df["time"] = z_to_time(np.array(df["redshift"]))
    dt = df["time"] - df["mm_prog_id"].replace(-1,pd.NA).map(df.set_index("id")["time"])
    print(dm)
    print(dt)
    
    


    df["dm"] = dm
    df["dt"] = dt
    dmdt = dm/dt
    df["dmdt"] = dmdt
    df["dyn_heating"] =  (df["temp"]/df["mass"] * consts.K_CONST/(5.0/3.0-1.0) * dmdt/consts.SEC_PER_MYR).clip(lower=0).fillna(0)

    df.to_csv("mcrit.txt")

    




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

'''


if __name__=="__main__":
    main()
