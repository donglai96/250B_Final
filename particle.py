from plasmapy.physics import parameters
from sympy import symbols, diff, exp, pi, cos, sin
from astropy import units as u
from sympy import besselj,DiracDelta


class particles():
    def __init__(self, name, density, magnetic_field, Tx, Tz):
        self.name = name
        self.density = density
        self.gyro_f =  parameters.gyrofrequency(B=magnetic_field, particle=name, signed=True)
        self.plasma_f = parameters.plasma_frequency(density, particle=name)
        self.Tx = Tx
        self.Tz = Tz
        self.vth_x = parameters.thermal_speed(Tx*u.K,name).value
        self.vth_z = parameters.thermal_speed(Tz*u.K,name).value

    def F_0(self):

        """
        Bi-Maxwellian
        """
        #x means perpendicular and z means parallel

        vx = symbols('vx')
        vz = symbols('vz')
        vth_x = symbols('vth_x')
        vth_z = symbols('vth_z')
        exp_term = exp((-vx**2/(vth_x)**2) +(-(vz)**2/(vth_z)**2))
        coeff_term = pi**(-3/2)*vth_x**-2 * vth_z**(-1)
        # exp_term = exp((-vx ** 2 / 8) + (-(vz) ** 2 /  2))
        # coeff_term = (2*pi) ** (-3 / 2) *(128)**(-1/2)

        return exp_term* coeff_term
    def F_1(self):

        """
        Bi-Maxwellian-Drift
        """
        #x means perpendicular and z means parallel

        vx = symbols('vx')
        vz = symbols('vz')
        vth_x = symbols('vth_x')
        vth_z = symbols('vth_z')
        exp_term = exp((-vx**2/(vth_x)**2) +(-(vz+4e7)**2/(vth_z)**2))
        coeff_term = pi**(-3/2)*vth_x**-2 * vth_z**(-1)
        return exp_term* coeff_term
    def G_1(self,distribution):
        """
        Return to 2 term of G1,
        The first is df/dvx
        The second is vz*df/dvx - vx*df/dz

        """
        if distribution =='drift':
            F = self.F_1()
        else:
            F = self.F_0()
        vx = symbols('vx')
        vz = symbols('vz')
        g1_term1 = diff(F,vx)
        g1_term2 =( vz * diff(F,vx) - vx*diff(F,vz))
        return g1_term1, g1_term2

    def Jm_p(self,m):
        """
        The argument of Jm change to r
        r = k_perp*v_perp/Omega
        @param m: the order of Bessel function
        @return:
        """
        # The resaon use r here is because Jm(r) is easy to integrate
        r = symbols('r')
        argument = r
        return besselj(m,argument)

    def delta_function(self, m):
        """
        The argument of Dirac is v_para - (w -m*Omega/k_parallel)
        @param m: The order
        @return: Dirac function of argument
        """
        vz = symbols('vz')
        k_para = symbols('k_para')
        w = symbols('w')
        argument = vz - (w - m * self.gyro_f.value) / k_para
        return DiracDelta(argument)

    def weight_function_p(self,m):
        theta = symbols('theta')
        gamma_m = (((1 + cos(theta))*self.Jm_p(m+1) + (1-cos(theta))*self.Jm_p(m-1))/(2*cos(theta)))**2

        return gamma_m