# test for bi_maxwellian distribution
# electron = particle.particles('e',n[0],B,400,100)
# print(electron.F_0())
# theta_input = 0
# theta_deg = theta_input*pi/180
# wave_frequency = 0.01 *gyro_frequency
# refraction_index_para = disper.solve_dispersion(B,species,n,wave_frequency, theta_input)
# vx = symbols('vx')
# vz = symbols('vz')
# Tx = symbols('Tx')
# Tz = symbols('Tz')
# f = electron.F_0()
# result_x= integrate(f,(vz,-oo,oo))
# result = integrate(2 * pi*vx *result_x,(vx,0,oo))

# set w and Omega and k
# Omega = gyro_frequency.value
# w = wave_frequency.value
# Omega_dot = Omega/Omega
# c = const.c.value
# print("light.speed",c)
# k = refraction_index_para*w /c
# k_para = k*cos(theta_deg)
# k_para_c = k_para*c
# nn = refraction_index_para
#
# print(result)
# G1_term1, G1_term2 = electron.G_1()
# print(G1_term1,G1_term2)
# print(refraction_index_para)
# G1 = G1_term1 - (k_para/w) * G1_term2
# print(G1)
# print("vz.c",(w-Omega)/k_para)
# G1_res = G1.subs(vz,(w-Omega)/(k_para))
# print(G1_res)
# # wave_growth_para_coeff= pi**2 * Omega*c/n
# wave_para_int = integrate(vx**2 * G1,(vx,0,oo))
# print(wave_para_int)
#%%
# Test of the parallel n the course reader
#
# A_res_top = integrate(vx**2/vz * G1_term2,(vx,0,oo))
# print(A_res_top)
# A_res_down = 2 *integrate(vx*electron.F_0(),(vx,0,oo))
# print(A_res_down)
# Eta_res = pi*(-vz) *A_res_down
# print(k_para,"K")
# Eta_value = Eta_res.subs(vz,(w-Omega)/(k_para))
# print(Eta_res,Eta_value)
# A_c = 1/(Omega/w - 1)
# A_res = A_res_top/A_res_down
# growth_rate_course = (k/w)*pi**-1 *Eta_res*(A_res_top/A_res_down - A_c)
# print("test")
# from sympy import simplify
# print((growth_rate_course/wave_para_int).subs(vz,(w-Omega)/(k_para)))
#%%
# G1 term test
# G1_term_11 = integrate(vx**2 * G1_term1,(vx,0,oo))
# G1_term_22 = integrate(vx * electron.F_0(),(vx,0,oo))
# print(G1_term_11,G1_term_22)
#
# print("test for normalization",growth_rate_course.subs(vz,(w-Omega)/(k_para_c)))