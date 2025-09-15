import numpy as np
from setup_run import constants
consts = constants()


def calc_omega_m_z(z):
    return consts.OMEGA_M*np.float_power(1+z,3)/ (consts.OMEGA_M*np.float_power(1+z,3) + consts.OMEGA_L + consts.OMEGA_K*np.float_power(1+z,2))

def calc_T_vir(m_halo,z): #Returns in units of Kelvin, eq 26 from https://arxiv.org/pdf/astro-ph/0010468.pdf, this function takes in units of m_halo (in Msun, not Msun/h) and redshift z
    omz=calc_omega_m_z(z)
    d = omz-1
    delta_c = 18*np.pi*np.pi + 82*d -39*d*d
    a = 1.98e4*(consts.MU_CONSTANT/0.6)*np.float_power(m_halo*consts.BIG_H/1e8, 2.0/3)
    b = np.float_power(consts.OMEGA_M/omz * delta_c/(18*np.pi*np.pi),1.0/3) 
    c = (1+z)/10.0
    return a*b*c



def calc_R_vir(m_halo,z): #returns kpc/h
    omz=calc_omega_m_z(z)
    d = omz-1
    delta_c = 18*np.pi*np.pi + 82*d -39*d*d
    return 0.784* np.float_power(m_halo*consts.BIG_H/1e8, 1.0/3)*np.float_power(consts.OMEGA_M/omz * (delta_c/(18*np.pi*np.pi)), -1.0/3) * (10/(1+z))


def z_to_time(z):
        a=2./(3.*np.sqrt(consts.OMEGA_L))
        b = np.sqrt(consts.OMEGA_L/consts.OMEGA_M)*np.float_power(1+z, -1.5)
        c = 71*1e5*3.154e+7/3.086e18
        return a/c*np.log(b+np.sqrt(1+b*b))

