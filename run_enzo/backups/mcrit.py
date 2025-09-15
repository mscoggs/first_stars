import numpy as np
import sys
import glob
import yt
import ytree
import os
import pandas as pd
from setup_run import constants
from cosmo import *
consts = constants()



def apply_sf_condition(ds, halo_position,halo_radius,T_virs):
    star_forming_k, star_forming_s = [], []

    left= (ds.index.grid_left_edge.to("Mpccm/h")).value
    right=(ds.index.grid_right_edge.to("Mpccm/h")).value
    levels= ds.index.grid_levels.flatten()

    for N in range(np.size(T_virs)):
        print("Working on halo", N,"out of",np.size(T_virs), "| Total sf:", np.sum(star_forming_k))
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

    print(np.sum(star_forming_k))
    print(np.sum(star_forming_s))

    return star_forming_k, star_forming_s



def get_mcrit_kulkarni(x,y):
   #CHANGE MIIN MASS
    bs_start = 0.05*consts.MIN_MASS_MCRIT
    min_num_bins = 3
    max_bs = np.max(x)/min_num_bins
    max_its = int(np.ceil(max_bs/bs_start))
    LOW_START, UP_START = 0.25,0.75
    OFFSET_NUMS=10

    for bound_offset in np.arange(0,0.25,0.05):
        LOW, UP = LOW_START-bound_offset, UP_START + bound_offset

        for bs_factor in range(1, max_its+1):
            for offset in range(OFFSET_NUMS):
                bs = bs_start*bs_factor
                l = consts.MIN_MASS_MCRIT - bs*offset/OFFSET_NUMS
                r = l+bs
                f=[]
                m=[]
                while(l<np.max(x)):
                    args = np.where((l<=x) & (x<r))
                    if(np.size(args)==0):
                        l=r
                        r+=bs
                        continue
                    else:
                        f.append(np.sum(y[args])*1.0/np.size(args))
                        m.append((l+r)/2.0)
                        l=r
                        r+=bs
                f,m = np.array(f), np.array(m)
                filt = np.where((f>= LOW)&(f<= UP))
                reduced_f = f[filt]
                reduced_m = m[filt]
                if(np.size(reduced_f) <=1): #too few points
                    continue
                maxm,minm = np.max(reduced_m), np.min(reduced_m)
                filt2 = np.where((m >=minm) & (m <= maxm))
                reduced_f = f[filt2]
                reduced_m = m[filt2]
                mono = True
                for i in range(np.size(reduced_f)-1):
                    if(reduced_f[i+1]<= reduced_f[i]):mono = False
                if(mono):
                    filt = np.where((f>= LOW)&(f<= UP))
                    if(np.min(filt)-1 >= 0):
                        filt = np.append(filt, np.array([np.min(filt)-1]))
                    if(np.max(filt)+1 < np.size(f)):
                        filt = np.append(filt, np.array([np.max(filt)+1]))
                    reduced_f = f[filt]
                    reduced_m = m[filt]
                    mcrit = np.interp(0.5, reduced_f, reduced_m)
                    return mcrit, bs/2, f, m
    return np.min(np.array(x)[np.where(np.array(y)==1)]), -1, 0,0




def main():
    args = sys.argv[1:]
    print(args)
    if(len(args) != 5):
        print("Error, improper number of arguments")
        sys.exit(0)
    enzo_dir,halo_dir,out_dir, min_mass, save_tail = args
    min_mass = float(min_mass)
    sim_num = enzo_dir.split("sim")[-1]
    halo_type = (halo_dir.split("/")[-1]).split("_")[0]
    snap_nums = [snapstr[-4:] for snapstr in list(glob.glob(enzo_dir + "/RD*"))]

    names =['z','n_halos', 'mcrit_kulkarni', 'mcrit_std_kulkarni', "mmin_kulkarni"  'mave_kulkarni', "mave_std_kulkarni", "n_sf_kulkani", 'mcrit_schauer', 'mcrit_std_schauer', "mmin_schauer", 'mave_schauer', "mave_std_schauer", "n_sf_schauer"]
    df = pd.DataFrame(columns=names)
    

    for enzo_snap in snap_nums:
        enzo_snap = "0013"
        print("Working on snap", enzo_snap,"\n\n")
        ds  = yt.load(enzo_dir + "/RD" + enzo_snap + "/RedshiftOutput"+enzo_snap)
        u = ds.units
        z = ds.current_redshift
        if("rockstar" in halo_dir): 
            #halo_file = halo_dir + "/halos_RD"+enzo_snap+".0.bin"
            out_file = "out_"+str(int(enzo_snap.split("0")[-1])-1)+".list"
            print(out_file)
            os.system("cp "+halo_dir+"/"+out_file + " " + enzo_dir +"/RD" + enzo_snap)
            halo_file = enzo_dir+"/RD"+enzo_snap+"/"+out_file
            halos = ytree.load(halo_file)
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
        print("saving data")
        df = pd.DataFrame(np.array([halo_mass, star_forming_k, star_forming_s]).T,columns=["mass", "sfk", "sfs"])
        df.to_csv("test"+enzo_snap+".txt")
        sys.exit(0)
        mcrit1, std1, f1, m1 = get_mcrit_kulkarni(np.array(halo_mass),np.array(star_forming_k))
        mcrit2, std2, f2, m2 = get_mcrit_kulkarni(np.array(halo_mass),np.array(star_forming_s))

        
        m_sf1 = halo_mass[np.where(np.array(star_forming_k) == 1)]
        m_sf2 = halo_mass[np.where(np.array(star_forming_s) == 1)]
        data = [z, np.size(halo_mass), mcrit1, std1, np.min(m_sf1), np.mean(m_sf1), np.std(m_sf1), np.size(m_sf1),  mcrit2, std2, np.min(m_sf2), np.mean(m_sf2), np.std(m_sf2), np.size(m_sf2)]
        df =pd.concat([df, pd.DataFrame(data)], ignore_index=True)

    df.to_csv(out_dir+"/mcrit_"+sim_num+"_"+halo_type+"_"+save_tail+".txt")


if __name__=="__main__":
    main()
