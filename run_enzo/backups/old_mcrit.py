import numpy as np
import math
import sys
import glob
import yt
from yt.units import kpc
import ytree
from astropy.cosmology import WMAP7   # WMAP 7-year cosmology
from astropy import units as u
import astropy.cosmology.units as cu
import os
import shutil
import pandas as pd






def get_enzo_dir(tr):
    candidates = list(glob.iglob(BASE_DIR+ ENZO_DIR_FORMAT))
    word = 'CosmologyCurrentRedshift'
    for file in candidates:
        with open(file, 'r') as fp:
            lines = fp.readlines()
            for row in lines:
                if word in row:
                    this_redshift = float(row.split()[-1])
                    if(np.abs(this_redshift-tr) < REDSHIFT_TOLERANCE):
                        return file
    print("error, no dir for redshift", tr)
    return -1


def T_vir(m,z):
    h=H
    u = 1.22
    #eq 26 from https://arxiv.org/pdf/astro-ph/0010468.pdf
    #returns virial temp in units of Kelvin
    OMEGA_K = 1-OMEGA_M-OMEGA_L
    OMEGA_M_Z=OMEGA_M*np.float_power(1+z,3)/ (OMEGA_M*np.float_power(1+z,3) + OMEGA_L + OMEGA_K*np.float_power(1+z,2))
    d = OMEGA_M_Z-1

    delta_c = 18*np.pi*np.pi + 82*d -39*d*d
    #a = 1.98e4*(u/0.6)*np.float_power(m*h/1e8, 2.0/3)*np.float_power(OMEGA_M/OMEGA_M_Z * delta_c/(18*np.pi*np.pi),1.0/3) *(1+z)/10.0
    a = 1.98e4*(u/0.6)*np.float_power(m*h/1e8, 2.0/3)
    b = np.float_power(OMEGA_M/OMEGA_M_Z * delta_c/(18*np.pi*np.pi),1.0/3) 
    c = (1+z)/10.0
    d = a*b*c
    return d



def r_vir(m_h,z):
    #returns kpc/h
    h=H
    OMEGA_K = 1-OMEGA_M-OMEGA_L
    OMEGA_M_Z=OMEGA_M*np.float_power(1+z,3)/ (OMEGA_M*np.float_power(1+z,3) + OMEGA_L + OMEGA_K*np.float_power(1+z,2))
    d = OMEGA_M_Z-1
    delta_c = 18*np.pi*np.pi + 82*d -39*d*d
    return 0.784* np.float_power(m_h*h/1e8, 1.0/3)*np.float_power(OMEGA_M/OMEGA_M_Z * (delta_c/(18*np.pi*np.pi)), -1.0/3) * (10/(1+z))


def apply_sf_condition(ds, halo_position,halo_radius,T_virs):
    max_density = []
    max_temp = []
    star_forming_k = []
    star_forming_s = []
    star_forming_cm = []
    for N in range(np.size(T_virs)):
        r = halo_radius[N]*2
        left = (ds.index.grid_left_edge.to("Mpccm/h")).value-r
        right = (ds.index.grid_right_edge.to("Mpccm/h")).value+r
        args = np.where(np.all(halo_position[N] >= left,axis=1) & np.all(halo_position[N] <= right, axis=1))
        levels= ds.index.grid_levels[args]
        n_grids = np.size(levels)
        H2_fraction = ds.index.grids[args][np.argmax(levels)]["gas", "H2_p0_fraction"]
        n_density = ds.index.grids[args][np.argmax(levels)]["gas", "number_density"]
        temp = ds.index.grids[args][np.argmax(levels)]["gas", "temperature"]
        cm = ds.index.grids[args][np.argmax(levels)]["gas", "cell_mass"].to("Msun/h")
        #ds.fields
        density_args = np.where(n_density>=100)
        temp_at_density_args = temp[density_args]
        if(np.size(density_args) == 0): 
            star_forming_k.append(0)
            star_forming_s.append(0)
        else:
            temp_at_density_args = temp[density_args]
            H2_fraction_at_density_args = H2_fraction[density_args]
            cm_at_density = cm[density_args]
            if(np.any(temp_at_density_args<=T_virs[N]*0.5)):
                star_forming_k.append(1)
                print("MEAN CM: \n", np.min(cm_at_density[np.where(temp_at_density_args<=T_virs[N]*0.5)]))
                print("\n\n")
                star_forming_cm.append(np.min(cm_at_density[np.where(temp_at_density_args<=T_virs[N]*0.5)]))
            else:
                star_forming_k.append(0)

            if((np.any((temp_at_density_args<500) & (H2_fraction_at_density_args>= 1e-4)))):
                star_forming_s.append(1)
            else: 
                star_forming_s.append(0)

    return star_forming_k, star_forming_s, star_forming_cm



def get_mcrit_kulkarni(x,y):
    print(np.size(x), np.size(y))
    print(list(x))
    print(list(y))
    bs_start = 0.05*MIN_MASS
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
                l = MIN_MASS- bs*offset/OFFSET_NUMS
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
    #print("\n NO SUCCESS IN MCRIT, USING MMIN")
    return np.min(np.array(x)[np.where(np.array(y)==1)]), -1, 0,0







global MIN_MASS, HALO_DEF, BASE_DIR, ENZO_DIR_FORMAT, ROCKSTAR_DIR_FORMAT, OMEGA_M, OMEGA_L, H, REDSHIFT_TOLERANCE
BASE_DIR = "/scratch/08288/tg875874/enzo_sims/sim"
ENZO_DIR_FORMAT = "RD00??/RedshiftOutput????"
ROCKSTAR_DIR_FORMAT = "rockstar_halos/out*" #"rockstar_halos/out*"
HALO_DEF = "RS" #RS for rockstar, FOF for fof
MIN_MASS = 1e5 #note that anna used 5e4
OMEGA_M, OMEGA_L,H = 0.32, 0.68, 0.67    # H = ds.hubble_constant           #     OMEGA_M = ds.omega_matter #     OMEGA_L = ds.omega_lambda
REDSHIFT_TOLERANCE = 0.01


def main():
    global MIN_MASS, HALO_DEF, BASE_DIR

    args = sys.argv[1:]
    SIMNUM = args[0]
    BASE_DIR += SIMNUM + "/"
    if(len(args)==4):
        HALO_DEF = args[1]
        MIN_MASS = float(args[2])
        SAVENAME = str(args[3])
        #if(MASS_DEF not in ["M200c", "Mvir", "Mvir_all"]):
        #    print("Error, mass definition", MASS_DEF, "not found")
        #    sys.exit(0)
    elif(len(args) != 0):
        print("Error, improper number of arguments")#. Pass 3 (sf criteria, mass def, and min mass) or 0 (assuming defaults)")
        sys.exit(0)

    

    target_redshift = [30,27,25,24,23,22,21,20,19,18,17,16,15]
    #target_redshift = [17,16,15]
    redshift = []

    mcrit_list1,mcrit_std_list1 = [],[]
    mcrit_list2,mcrit_std_list2 = [],[]
    mmin_list1 = []
    mmin_list2 = []
    mave_list1, mave_std_list1 = [],[]
    mave_list2, mave_std_list2 = [],[]
    n_halos = []
    n_sf1, n_sf2 = [],[]

    print("Running mcrit for target redshift", target_redshift, "with halo def, and min mass:", HALO_DEF, MIN_MASS,"\n\n")
    for tr in target_redshift:
        print("Working on redshift", tr,"\n\n")
        enzo_dir = get_enzo_dir(tr)
        enzo_snap = enzo_dir[-4:]
        ds  = yt.load(enzo_dir)
        if(HALO_DEF == "RS"): halo_file = BASE_DIR +  "rockstar_halos/halos_RD"+enzo_snap+".0.bin"
        if(HALO_DEF == "FOF"): halo_file = BASE_DIR +  "fof_halos/RedshiftOutput"+enzo_snap+"/RedshiftOutput"+enzo_snap+".0.h5"

        halos = yt.load(halo_file).all_data()
        u = ds.units
        z = ds.current_redshift
        
        halo_position = halos['particle_position'].to("Mpccm/h").value
        halo_radius = halos['virial_radius'].to("Mpccm/h").value
        halo_mass = halos['particle_mass'].to("Msun").value

        
        halo_args = np.where(halo_mass >= MIN_MASS)
        halo_mass = halo_mass[halo_args]
        halo_position = halo_position[halo_args]
        print("total number of haloes above min mass:", np.size(halo_mass), "\n\n")
        if(np.size(halo_mass) == 0):
            print("no halos for redshift", tr)
            continue

        T_virs =  (T_vir(halo_mass,z)*u.K).value

        print(np.size(T_virs), np.size(halo_mass))
        star_forming_k, star_forming_s, star_forming_cm = apply_sf_condition(ds, halo_position,halo_radius,T_virs)
        print("total number of star forming k:", np.sum(star_forming_k), "\n\n")
        print("total number of star forming s:", np.sum(star_forming_s), "\n\n")
        print(np.mean(star_forming_cm))
        if(np.sum(star_forming_k)+np.sum(star_forming_s) == 0):
            print("no star forming halos", tr)
            continue
        mcrit1, std1, f1, m1 = get_mcrit_kulkarni(np.array(halo_mass),np.array(star_forming_k))
        mcrit2, std2, f2, m2 = get_mcrit_kulkarni(np.array(halo_mass),np.array(star_forming_s))

        mcrit_list1.append(mcrit1)
        mcrit_std_list1.append(std1)
        m_sf1 = halo_mass[np.where(np.array(star_forming_k) == 1)]
        mmin_list1.append(np.min(m_sf1))
        mave_list1.append(np.mean(m_sf1))
        mave_std_list1.append(np.std(m_sf1))


        mcrit_list2.append(mcrit2)
        mcrit_std_list2.append(std2)
        m_sf2 = halo_mass[np.where(np.array(star_forming_s) == 1)]
        mmin_list2.append(np.min(m_sf2))
        mave_list2.append(np.mean(m_sf2))
        mave_std_list2.append(np.std(m_sf2))

        redshift.append(tr)
        n_halos.append(np.size(halo_mass))
        n_sf1.append(np.sum(star_forming_k))
        n_sf2.append(np.sum(star_forming_s))
        

    df = pd.DataFrame({'z': redshift, 
        'mcrit_kulkarni': mcrit_list1, 'mcrit_std_kulkarni': mcrit_std_list1, "mmin_kulkarni" : mmin_list1, 'mave_kulkarni': mave_list1, "mave_std_kulkarni": mave_std_list1, 
        'mcrit_schauer': mcrit_list2, 'mcrit_std_schauer': mcrit_std_list2, "mmin_schauer" : mmin_list2, 'mave_schauer': mave_list2, "mave_std_schauer": mave_std_list2, 
        "n_halos": n_halos, "n_sf_kulkani": n_sf1, "n_sf_schauer": n_sf2, "m_ave_cell_sf": np.mean(star_forming_cm)})
    df.to_csv("mcrit_sim_"+str(SIMNUM)+"_halo_def_"+HALO_DEF+"_"+str(MIN_MASS)+"_"+str(SAVENAME)+".txt")


if __name__=="__main__":
    main()
