from astropy import constants as const
from plasmapy.physics import parameters
from astropy import units as u
from sympy import pi,oo,simplify
from sympy import integrate, symbols
import disper
import numpy as np
import matplotlib.pyplot as plt
import particle
from astropy import constants as const
from sympy import sin,cos,exp

def wave_growth_para(name,density,B,T_perp,T_para,waveperOmega,distribution):
    if name == 'e':
        n =density[0]
    else:
        n = density[1]

    electron = particle.particles(name,n,B,T_perp,T_para)
    species = ['e','p']
    theta_input = 0
    theta_deg = theta_input * pi / 180

    gyro_frequency = parameters.gyrofrequency(B, name)
    wave_frequency = waveperOmega * gyro_frequency

    refraction_index_para = disper.solve_dispersion(B, species, density, wave_frequency, theta_input)
    print("nn",refraction_index_para )
    vx = symbols('vx')
    vz = symbols('vz')
    vth_x = symbols('vth_x')
    vth_z = symbols('vth_z')
    # Tx = symbols('Tx')
    # Tz = symbols('Tz')
    f = electron.F_0()
    # result_x = integrate(f, (vz, -oo, oo))
    # result = integrate(2 * pi * vx * result_x, (vx, 0, oo))

    # set w and Omega and k
    Omega = gyro_frequency.value
    w = wave_frequency.value

    c = const.c.value
    k = refraction_index_para * w / c
    k_para = k * cos(theta_deg)
    k_para_c = k_para*c


    #print(result)
    G1_term1, G1_term2 = electron.G_1(distribution)
    print(G1_term1,G1_term2,f)
   # print(G1_term1, G1_term2)
    #print(refraction_index_para)
    G1 = G1_term1 - (k_para / w) * G1_term2

    G1_para =  - (k_para / w) * G1_term2
    #
    vz_res = (w - Omega) / (k_para)
    G1_new = (G1.subs(vz,vz_res))
    print("test!!")
    print(G1_new)
    # G1_res = G1_para.subs(vz, (w - Omega) / (k_para))
    #print(G1_res)
    # wave_growth_para_coeff= pi**2 * Omega*c/n
    wave_para_int =(pi**2 * Omega * w / k)* integrate(vx ** 2 * G1_new, (vx, 0, oo))
    wave_para_int_new = wave_para_int.subs([(vth_x,electron.vth_x),(vth_z,electron.vth_z)])
    print(wave_para_int_new)
    #wave_growth_eva = wave_para_int.subs(vz, vz_res)
    #print(wave_para_int)
    return wave_para_int_new, vz_res, G1_new

