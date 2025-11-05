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



def apply_sf_condition(ds, halo_position,halo_radius,T_virs):
    star_forming_k, star_forming_s = [], []

    left= (ds.index.grid_left_edge.to("Mpccm/h")).value
    right=(ds.index.grid_right_edge.to("Mpccm/h")).value
    levels= ds.index.grid_levels.flatten()
    #print("Working on halo", N,"out of",np.size(T_virs), "| Total sf:", np.sum(star_forming_k))

    for N in range(np.size(T_virs)):
        r = halo_radius[N]*2
        args = np.where(np.all(halo_position[N] >= left-r,axis=1) & np.all(halo_position[N] <= right+r, axis=1))
        halo_levels= levels[args]
        max_levels_args = np.where(halo_levels == np.max(halo_levels))[0]
        target_grids = ds.index.grids[args][max_levels_args]

        sfk,sfs = 0,0
        for grid in target_grids:
            H2_fraction = grid["gas", "H2_p0_fraction"].value
            n_density = grid["gas", "number_density"].value
            temp = grid["gas", "temperature"].value
            density_args = np.where(n_density>=100)
            if(np.size(density_args) == 0): continue

            temp = temp[density_args]
            H2_fraction = H2_fraction[density_args]
            if(np.any(temp<=T_virs[N]*0.5)): sfk = 1    #kulkarni conditions
            if(np.any((temp<500) & (H2_fraction>= 1e-4))): sfs = 1 #schauer conditions

        star_forming_k.append(sfk)
        star_forming_s.append(sfs)

    return np.array(star_forming_k), np.array(star_forming_s)



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








def main():
    args = sys.argv[1:]
    if(len(args) != 5):
        print("Error, improper number of arguments")
        sys.exit(0)
    enzo_dir,halo_dir,out_dir, min_mass, save_tail = args
    min_mass = float(min_mass)
    sim_num = enzo_dir.split("sim")[-1]
    halo_type = (halo_dir.split("/")[-1]).split("_")[0]
    snap_nums = [snapstr[-4:] for snapstr in list(glob.glob(enzo_dir + "/RD*"))]

    names =['z','n_halos', 'mcrit_k', 'mcrit_euk', 'mcrit_edk', "mmin_k",  'mave_k', "mave_std_k", "n_sf_kulkani", 'mcrit_s', 'mcrit_eus', 'mcrit_eds',"mmin_s", 'mave_s', "mave_std_s", "n_sf_s"]
    df = pd.DataFrame(columns=names)
    df_row = 0

    for enzo_snap in snap_nums:
        if("13" not in enzo_snap): continue 
        print(halo_type)
        print("Working on snap", enzo_snap,"\n\n")
        ds  = yt.load(enzo_dir + "/RD" + enzo_snap + "/RedshiftOutput"+enzo_snap)
        u = ds.units
        z = ds.current_redshift
        if("rockstar" in halo_dir): 
            #halo_file = halo_dir + "/halos_RD"+enzo_snap+".0.bin"
            out_file = "out_"+str(int(enzo_snap[-2:])-1)+".list"
            print("this it out:", out_file)
            os.system("cp "+halo_dir+"/"+out_file + " " + enzo_dir +"/RD" + enzo_snap)
            halo_file = enzo_dir+"/RD"+enzo_snap+"/"+out_file
            halos = ytree.load(halo_file)
            halo_vcirc = (np.sqrt(G_consts*halos["Mvir"].to("Msun").value*u.Msun/((np.array(halos["Rvir"])*u.kpccm/u.h).to("kpc")))).to("km/s")
            print(halo_vcirc)
            sys.exit(0)


            halo_mass = np.array(halos["M200c"].to("Msun"))
            halo_x = np.array(halos["X"])
            halo_y = np.array(halos["Y"])
            halo_z = np.array(halos["Z"])
            halo_position = np.array([halo_x,halo_y,halo_z]).T
            halo_radius = (np.array(halos["Rvir"])*u.kpccm/u.h).to("Mpccm/h").value


        elif("fof" in halo_dir): 
            halo_file = halo_dir + "/RedshiftOutput"+enzo_snap+"/RedshiftOutput"+enzo_snap+".0.h5"
            halos = yt.load(halo_file).all_data()
            halo_mass = halos['particle_mass'].to("Msun").value
            halo_position = halos['particle_position'].to("Mpccm/h").value
            halo_radius  = halos['virial_radius'].to("Mpccm/h").value
        halo_args = np.where(halo_mass >= min_mass)

        if(np.size(halo_args) == 0):
            print("no halos for redshift", z)
            continue

        halo_mass = halo_mass[halo_args]
        halo_position = halo_position[halo_args]
        halo_radius = halo_radius[halo_args]


        T_virs = (calc_T_vir(halo_mass,z)*u.K).value
        

        star_forming_k, star_forming_s  = apply_sf_condition(ds, halo_position,halo_radius,T_virs)
        #print("saving data")
        #df = pd.DataFrame(np.array([halo_mass, star_forming_k, star_forming_s]).T,columns=["mass", "sfk", "sfs"])
        #df.to_csv("test"+enzo_snap+".txt")
        #sys.exit(0)
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


if __name__=="__main__":
    main()
