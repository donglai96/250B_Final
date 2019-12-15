from plasmapy.physics import parameters
from sympy import symbols, diff, exp, pi
from astropy import units as u

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
        Tx = symbols('Tx')
        Tz = symbols('Tz')
        exp_term = exp((-vx**2/(vth_x)**2) +(-(vz+4e7)**2/(vth_z)**2))
        coeff_term = pi**(-3/2)*vth_x**-2 * vth_z**(-1)
        return exp_term* coeff_term
    def F_1(self):

        """
        Drift_Bi-Maxwellian
        """
        #x means perpendicular and z means parallel

        vx = symbols('vx')
        vz = symbols('vz')
        v_drift = symbols('v_drift')
        Tz = symbols('Tz')
        exp_term = exp((-vx**2/self.Tx) +(-vz**2/self.Tz))
        coeff_term = pi**(-3/2)*self.Tx**-1 * self.Tz**(-1/2)
        return exp_term* coeff_term

    def G_1(self):
        """
        Return to 2 term of G1,
        The first is df/dvx
        The second is vz*df/dvx - vx*df/dz

        """
        vx = symbols('vx')
        vz = symbols('vz')
        g1_term1 = diff(self.F_0(),vx)
        g1_term2 =( vz * diff(self.F_0(),vx) - vx*diff(self.F_0(),vz))
        return g1_term1, g1_term2
