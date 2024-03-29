from sympy import symbols,pi
from astropy import units as u
from plasmapy.physics import parameters
from sympy import Matrix,simplify,solve,root
from sympy import sin, cos
from scipy import special






def dispersion_matrix():
    L = symbols('L')
    R = symbols('R')
    P = symbols('P')
    theta_rad = symbols('theta')
    # Only solve the positive refraction
    nn = symbols('nn',positive = True)
    m_11 = 2 * (R - nn ** 2 + 0.5 * nn ** 2 * (sin(theta_rad)) ** 2)
    m_12 = nn ** 2 * sin(theta_rad) ** 2
    m_13 = nn ** 2 * cos(theta_rad) * sin(theta_rad)
    m_21 = m_12
    m_22 = 2 * (L - nn ** 2 + 0.5 * nn ** 2 * (sin(theta_rad)) ** 2)
    m_23 = m_13
    m_31 = m_13
    m_32 = m_13
    m_33 = P - nn ** 2 * (sin(theta_rad) ** 2)
    D_0 = Matrix([[m_11, m_12, m_13], [m_21, m_22, m_23], [m_31, m_32, m_33]])
    return D_0

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u


import math
from astropy import constants as const
from scipy.optimize import fsolve
import scipy.special as spl

def cold_plasma_LRP(B, species, n, omega ,flag= True ):
    """

    @param B: Magnetic field in z direction
    @param species: ['e','p'] means electrons and protons
    @param n: The density of different species
    @param omega: wave_frequency
    @param flag: if flag == 1 then omega means to be a number, if flag is not True the omega means a symbol
    @return: (l, R, P)

    """

    if flag:
        wave_frequency = omega
        w = wave_frequency.value
    else:
        wave_frequency = symbols('w')
        w = wave_frequency
    L, R, P = 1, 1, 1

    for s, n_s in zip(species, n):
        omega_c = parameters.gyrofrequency(B=B, particle=s, signed=True)
        omega_p = parameters.plasma_frequency(n=n_s, particle=s)

        L += - omega_p.value ** 2 / (w * (w - omega_c.value))
        R += - omega_p.value ** 2 / (w * (w + omega_c.value))
        P += - omega_p.value ** 2 / w ** 2
    return L, R, P


def solve_dispersion(B, species, n, wave_frequency, theta_input):
    L = symbols('L')
    R = symbols('R')
    P = symbols('P')
    theta = symbols('theta')



    # Get the function
    theta_deg = theta_input *pi/180
    D_0 = dispersion_matrix()
    L_wave, R_wave, P_wave = cold_plasma_LRP(B, species, n, wave_frequency)
    #print("L,R,P is ", L_wave, R_wave, P_wave)
    f = simplify(D_0.det())
   # print("D_0.det is test")
    #print(f.subs([(L, L_wave), (P, P_wave), (R, R_wave), (theta, theta_deg)]))
    D_0_sim = f.subs([(L, L_wave), (P, P_wave), (R, R_wave), (theta, theta_deg)])
    refractive_index = solve(D_0_sim)
    print("nn",refractive_index)
    return refractive_index[0]

def faddeeva(x):
    """
    Faddeeva function is defined as

    w(z) = exp(-z**2) * erfc(-i * z)
    Computation of the Complex Error Function J.A.C. Weideman

    The faddeeva function is wrapped in C++ in scipy


    """
    return special.wofz(x)

def zeta(x):
    zeta = faddeeva(x) * 1j * np.sqrt(np.pi)
    return zeta
def wave_growth_Langmuir(k, c):
    """
    This function is for calculating the dispersion relationship
    f is maxwellian distribution
    c is initial value
    """
    w = []
    for kk in k:

        def f(x):
            f = 1 + kk ** 2 + x * zeta(x)
            return f
        print(f(1+1j))
        print(c)

        b = fsolve(f,c)
        print(b)
        x = c*np.sqrt(2)*kk
        w = np.append(w,x)


    return w
