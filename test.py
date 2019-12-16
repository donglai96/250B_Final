#%%
from astropy import constants as const
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
from scipy.optimize import fsolve
from astropy import constants as const
from sympy import sin,cos,exp

species = ['e','p']

B = 1e-6 * u.T
n = [300e6* u.m ** -3, 300e6 * u.m ** -3]
theta_input = 0
theta_deg = theta_input *pi/180
gyro_frequency = parameters.gyrofrequency(B, 'e')
wave_frequency = 0.05 * gyro_frequency
print(wave_frequency)

#%%

# test dispersion  solve
refraction_index = disper.solve_dispersion(B,species,n,wave_frequency,theta_input)
print(refraction_index)
# Draw a picture of n and w/Omega

x = np.arange(0.1,0.5,0.05)
refraction_index_array = np.zeros(len(x))
print(refraction_index_array)
refraction_index_array = []
for w in x:
    refraction_index = disper.solve_dispersion(B,species,n,w *gyro_frequency,theta_input)
    refraction_index_array = np.append(refraction_index_array,refraction_index)
print(refraction_index_array)

plt.plot(x,refraction_index_array)
plt.show()




#%%

B = 1e-6 * u.T
n = [300e6* u.m ** -3, 300e6 * u.m ** -3]
theta_input = 0
wave_over_Omega = 0.05
T_perp = 400
T_para = 100
distribution = 'drift'
electron = particle.particles('e',n,B,T_perp,T_para)
print("thermal",electron.vth_x,electron.vth_z)

# vz = symbols('vz')
growth_rate_base,vz_res_base,G0 = cold_growth_rate.wave_growth_para('e',n,B,T_perp,T_para,wave_over_Omega,distribution = 'drift')

print("The base of growth rate",growth_rate_base)

x = np.arange(0.1,0.5,0.05)

growth_rate_array_4 = []
vz_res_array = []
growth_value_base = growth_rate_base.evalf()

for w in x:
    growth_rate, vz_res, G1 =  cold_growth_rate.wave_growth_para('e',n,B,T_perp,T_para,w,distribution)

    print('v_res',vz_res)
    growth_rate_array_4 = np.append(growth_rate_array_4,growth_rate)
    vz_res_array= np.append(vz_res_array,vz_res)

print(growth_rate_array_4)
print(G1)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(x,growth_rate_array_4)
ax1.set_xlabel('w/Omega')
ax1.set_ylabel('growth_rate')
ax1.set_title('Growth rate and V_res')
ax2 = ax1.twinx()
ax2.plot(x,vz_res_array,linestyle = '-',label = 'v_res', color = 'red')
ax2.plot(x,-x *4e7/x,linestyle = '--',label = 'v_drift')
ax2.set_ylabel('Resonant Velocity')
fig.legend()
plt.savefig('whistler_growth.eps')
plt.show()
#
# growth_rate_array_2 = []

# Test for Langmiurt wave
# k_min = 0.1
# dk1 = 0.1
# k_mid = 1
# dk2 = 1
# kmax = 10
# k = np.append(np.arange(k_min,k_mid,dk1),np.arange(k_mid,kmax,dk2))
# c = 1
# print(k)
# w = disper.wave_growth_Langmuir(k,c)
# wre = np.real(w)
# wie = np.imag(w)
# print(wre,wie)