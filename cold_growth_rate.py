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

def anisotropy(name,density,B,T_perp,T_para,waveperOmega,distribution):
    if name == 'e':
        n =density[0]
    else:
        n = density[1]

    electron = particle.particles(name,n,B,T_perp,T_para)
    species = ['e','p']
    theta_input = 0
    theta_deg = theta_input * pi / 180
    vx = symbols('vx')
    vz = symbols('vz')
    vth_x = symbols('vth_x')
    vth_z = symbols('vth_z')
    G1_term1, G1_term2 = electron.G_1(distribution)
    if distribution == 'drift':
        f = electron.F_1()
    else:
        f = electron.F_0()

    gyro_frequency = parameters.gyrofrequency(B, name)
    wave_frequency = waveperOmega * gyro_frequency

    refraction_index_para = disper.solve_dispersion(B, species, density, wave_frequency, theta_input)

    # result_x = integrate(f, (vz, -oo, oo))
    # result = integrate(2 * pi * vx * result_x, (vx, 0, oo))

    # set w and Omega and k
    Omega = gyro_frequency.value
    w = wave_frequency.value

    c = const.c.value
    k = refraction_index_para * w / c
    k_para = k * cos(theta_deg)
    k_para_c = k_para * c

    vz_res = (w - Omega) / (k_para)
    A_bottom = 2 * integrate(f*vx, (vx, 0, oo))
    A_top = vz**-1 *integrate(vx**2 * G1_term2,(vx,0,oo))
    A_res = simplify((A_top/A_bottom).subs(vz,vz_res))
    A_res_value = A_res.subs([(vth_x,electron.vth_x),(vth_z,electron.vth_z)])
    Eta_res = 2*pi *(-vz)*integrate(vx*f,(vx,0,oo))
    Eta_res_v= Eta_res.subs(vz,vz_res)
    Eta_res_value = Eta_res_v.subs([(vth_x, electron.vth_x), (vth_z, electron.vth_z)])
    A_c = 1/((Omega - w ) -1)
    wave_growth = pi * Omega *(1 - w/Omega)**2 *Eta_res_value*(A_res_value - A_c)
    return A_res_value, Eta_res_value , wave_growth.evalf()

def wave_growth_oblique_m(name,density,B,T_perp,T_para,waveperOmega,distribution,theta_input, m ):
    """
    Calculate the m part of growth_rate
    Based on Kennel 1966
    """
    if name == 'e':
        n =density[0]
    else:
        n = density[1]

    electron = particle.particles(name,n,B,T_perp,T_para)
    species = ['e','p']
    theta_deg = theta_input * pi / 180
    # Here is positive signed = False
    gyro_frequency = parameters.gyrofrequency(B, name)
    wave_frequency = waveperOmega * gyro_frequency
    # First need to


    # Fist solve the refraction_index
    refraction_index_para = disper.solve_dispersion(B, species, density, wave_frequency, theta_input)

    Omega = electron.gyro_f.value
    ww = wave_frequency.value

    c = const.c.value
    kk = refraction_index_para * ww / c
    kk_para = kk * cos(theta_deg)
    kk_perp = kk * sin(theta_deg)

    # set the symbol
    vx = symbols('vx')
    vz = symbols('vz')
    vth_x = symbols('vth_x')
    vth_z = symbols('vth_z')
    r = symbols('r')
    w = symbols('w')
    k_para = symbols('k_para')
    theta = symbols('theta')


    # First integrate the v_parallel
    # we need four term
    # G1, weight function, Dirac Function, vx**2
    # First integrate v_parallel
    G1_term1, G1_term2 = electron.G_1(distribution)
    G1 = G1_term1 - (kk_para/ww) *G1_term2
    v_para_int = integrate(G1 * electron.delta_function(m).subs([(w,ww),(k_para,kk_para)]),(vz,-oo,oo))
    print(v_para_int)

    if electron.name == 'e':
        flag = -1
    else:
        flag = 1
    # flag determin the r integrate range -1 means -oo to 0
    # 1 means 0 to oo

    ### The reason here need to change that is because the Jm integration can not with parameters
    v_para_int_value = (v_para_int.subs([(vth_x,electron.vth_x),(vth_z,electron.vth_z)]))
    #v_perp_int_value = v_para_int_value * vx**2 * electron.weight_function_p(m)
    # # change the vx to (k *v / Omega) and input the vth_z and vth_x
    r_vx = r *Omega/kk_perp
    v_para_int_value_new = v_para_int_value.subs(vx,r_vx)


    v_perp_int = v_para_int_value_new*((Omega/kk_perp)**3 *r**2*electron.weight_function_p(m)).subs(theta,theta_deg)
    # Dothe simplify?
    #v_perp_int_test = simplify(v_perp_int)
    print("vy success")
   #Now do the final integration!
    if flag == -1:

        v_perp_int_value = integrate(v_perp_int,(r,-oo,0))
    else:
        v_perp_int_value = integrate(v_perp_int, (r, 0, oo))
    growth_rate_oblique = (- pi**2 * Omega *ww/kk) * v_perp_int_value
    return growth_rate_oblique.evalf()
    #return v_perp_int_test
