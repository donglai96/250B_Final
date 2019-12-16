from plasmapy.physics import parameters
from astropy import units as u
from sympy import pi,oo
from sympy import integrate, symbols
import disper
import numpy as np
import matplotlib.pyplot as plt
import cold_growth_rate

from sympy import ln
import particle


B = 1e-6 * u.T
n = [300e6* u.m ** -3, 300e6 * u.m ** -3]

wave_over_Omega = 0.38
T_perp = 400
T_para = 100
distribution = 'drft'

x = np.arange(0.1,0.5,0.005)

growth_rate_array_4 = []
vz_res_array = []

#
# for w in x:
#     growth_rate, vz_res, G1 =  cold_growth_rate.wave_growth_para('e',n,B,T_perp,T_para,w,distribution)
#
#     print('v_res',vz_res)
#     growth_rate_array_4 = np.append(growth_rate_array_4,growth_rate)
#     vz_res_array= np.append(vz_res_array,vz_res)
#
# print(growth_rate_array_4)
# print(G1)
m = 0
landau_res_array = []
x = np.arange(1,50,3)
for w in x:
    theta_input = w
    int_vz = cold_growth_rate.wave_growth_oblique_m('e',n,B,T_perp,T_para,wave_over_Omega,distribution,theta_input,m)
    landau_res_array = np.append(landau_res_array,int_vz)
    print(w)
print(landau_res_array)
plt.plot(x,landau_res_array)
plt.show()
#int_vz = cold_growth_rate.wave_growth_oblique_m(name,density,B,T_perp,T_para,waveperOmega,distribution,theta_input, m ):

# Draw a pic of m = 0 part with degree