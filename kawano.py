# -*- coding: utf-8 -*-

import os
import itertools
import numpy
import argparse
from scipy.integrate import simps, dblquad
from scipy.interpolate import interp1d
from subprocess import Popen, PIPE

from collections import namedtuple
from common import UNITS, CONST, utils
from common.integrators import integrate_1D
from library import SM

T_kawano = 10 * UNITS.MeV

q = (SM.particles.hadrons.neutron['mass'] - SM.particles.hadrons.proton['mass'])
# q = 1.2933 * UNITS.MeV

GF = CONST.G_F
m_n = SM.particles.hadrons.neutron['mass']
m_p = SM.particles.hadrons.proton['mass']
m_e = SM.particles.leptons.electron['mass']
gA = 1.27
Vud = 0.974
mun = -1.913
mup = 2.793
mA = 1.026 * UNITS.GeV
mV = 0.843 * UNITS.GeV
mPi = 0.13957 * UNITS.GeV
aa = 0.942
bb = 4.61

nu_transition = 20. * UNITS.MeV
a = None

Particles = namedtuple("Particles", "electron neutrino")
particles = None

sigma_nnu_to_pe_lowE = None
sigma_pnu_to_ne_lowE = None
sigma_nnu_to_pe_highE = None
sigma_pnu_to_ne_highE = None
sigma_pe_to_nnu_lowE = None
sigma_ne_to_pnu_lowE = None

def init_kawano(electron=None, neutrino=None):
    global particles
    particles = Particles(electron=electron, neutrino=neutrino)
    # import_interpolation_tables()


def run(data_folder, input="s4.dat", output="kawano_output.dat"):
    p = Popen(utils.getenv('KAWANO', 'KAWANO/kawano_noneq'), stdin=PIPE, env={
        "INPUT": os.path.join(data_folder, input),
        "OUTPUT": os.path.join(data_folder, output)
    })
    p.communicate(bytes(os.linesep.join([
        # ...
        "",
        # Run
        "4",
        # Go
        "2",
        # ...
        "",
        # Exit
        "4",
        # Output
        "5",
        # Request output file
        "1",
        # ...
        "", "", "", ""
    ]), 'utf-8'))
    with open(os.path.join(data_folder, output), "r") as kawano_output:
        return kawano_output.read()

def matrix_element(PePp, PnuPn, PePn, PnuPp, PePnu):
    return 16. * GF**2 * Vud**2 * ((1+gA)**2 * PePp * PnuPn + (1-gA)**2 * PePn * PnuPp + (gA**2-1)*m_n*m_p*PePnu)

def coulomb_correction(pe):

    beta = numpy.sqrt(1. - m_e**2/(pe**2 + m_e**2))
    Fcor = 2.*numpy.pi*CONST.alpha / beta / (1. - numpy.exp(-2.*numpy.pi*CONST.alpha/beta))

    return Fcor

# def radiative_correction(pnu, pe):
#     """ zero-temperature radiative correction """

#     y = pnu / m_e
#     eps = numpy.sqrt(pe**2 + m_e**2) / m_e
#     beta = numpy.sqrt(1. - 1/eps**2)
#     Rad = numpy.arctanh(beta) / beta

#     Ccor = 4.11542 + 4. * (Rad-1.) * (y/(3.*eps) -3./2. + numpy.log(2.*y)) + Rad * (2.*(1+beta**2) + (y**2)/(6.*eps**2) - 4.*beta*Rad) - 4. * (2. + 11.*beta + (224./9.)*beta**2 + (89./3.)*beta**3 + (1496./75.)*beta**4 + (596./75.)*beta**5 + (128./49.)*beta**6) / (1.+beta)**6

#     return CONST.rad_corr_prefac * (1. + CONST.alpha*Ccor / (2.*numpy.pi))

def dSigmadEe_highE(Enu, Ee, sign):
    Qsq = -(m_n**2 - m_p**2 - 2. * m_p * (Enu - Ee))
    tau = Qsq / (4. * m_p**2)
    FA = gA / (1 + Qsq/mA**2)**2
    FP = 2 * m_p**2 * FA / (Qsq + mPi**2)
    GD = 1 / (1 + Qsq/mV**2)**2
    FOneV = (1. + mun * (aa * tau) / (1. + bb * tau) + tau * (mup - mun)) * GD / (1. + tau)
    FTwoV = (mup - mun - 1. - mun * (aa * tau) / (1 + bb * tau)) * GD / (1. + tau)
    sMinu = 4. * m_p * Enu - Qsq - m_e**2
    partX = (m_e**2 + Qsq) * ((1. + tau) * FA**2 - (1. - tau) * FOneV**2 + tau * (1. - tau) * FTwoV**2 + 4. * tau * FOneV * FTwoV - (m_e**2) / (4. * m_p**2) * ((FOneV + FTwoV)**2 + (FA + 2. * FP)**2 - (Qsq / m_p**2 + 4.) * FP**2)) / m_p**2
    partY = Qsq * FA * (FOneV + FTwoV) / m_p**2
    partZ = 0.25 * (FA**2 + FOneV**2 + tau * FTwoV**2)
    return m_p**3 * GF**2 * Vud**2 * (partX + sign * sMinu * partY / m_p**2 + sMinu**2 * partZ / m_p**4) / (4. * numpy.pi * Enu**2)

def n_to_penu():
    """ n ⟶ p + ν_e' + e """

    def Integrand(Enu, Ee):
        return -8. * GF**2 * m_n * Vud**2 * (2. * Ee**2 * (gA-1)**2 * m_n - Ee * (gA-1) * (m_p * (2. * Enu * (gA+1) - gA * m_p + m_p) + (gA-1) * m_e**2 + (gA-1) * m_n**2) + Enu * (gA+1)**2 * (2. * Enu * m_n + m_e**2 - m_n**2 + m_p**2)) * (1. - particles.neutrino.distribution(Enu * a)) * (1. - particles.electron.distribution(numpy.sqrt(Ee**2-m_e**2) * a)) * coulomb_correction(numpy.sqrt(Ee**2 - m_e**2))# * radiative_correction(Enu, numpy.sqrt(Ee**2 - m_e**2))

    # return dblquad(Integrand, m_e, (m_n**2-m_p**2+m_e**2)/(2.*m_n), lambda x: (0.5 * (m_n**2 - m_p**2 + m_e**2) - m_n * x)/(m_n - x + numpy.sqrt(x**2 - m_e**2)), lambda x: (0.5 * (m_n**2 - m_p**2 + m_e**2) - m_n * x)/(m_n - x - numpy.sqrt(x**2 - m_e**2)))[0] / (4. * numpy.pi)**3 / m_n

    return integrate_1D(numpy.vectorize(lambda u1: integrate_1D(lambda u2: Integrand(u2, u1),
        bounds=((0.5 * (m_n**2 - m_p**2 + m_e**2) - m_n * u1)/(m_n - u1 + numpy.sqrt(u1**2 - m_e**2)),
        (0.5 * (m_n**2 - m_p**2 + m_e**2) - m_n * u1)/(m_n - u1 - numpy.sqrt(u1**2 - m_e**2))))[0]),
        bounds=(m_e, (m_n**2-m_p**2+m_e**2)/(2.*m_n)))[0] / (4. * numpy.pi)**3 / m_n

def penu_to_n():
    """ p + e + ν_e' ⟶  n """

    def Integrand(Enu, Ee):
        return  Ee * Enu**2 * numpy.sqrt(Ee**2 - m_e**2) * particles.neutrino.distribution(Enu * a) * particles.electron.distribution(numpy.sqrt(Ee**2-m_e**2) * a) * coulomb_correction(numpy.sqrt(Ee**2 - m_e**2))# * radiative_correction(Enu, numpy.sqrt(Ee**2 - m_e**2))

    return GF**2 * Vud**2 * (1 + 3. * gA**2) * integrate_1D(lambda x: Integrand(x, m_n - m_p - x), bounds=(0., m_n - m_p - m_e))[0] / 2. / numpy.pi**3

def pnu_to_ne():
    """ p + ν_e' ⟶  n + e """

    def dSigmadEe_lowE(Enu, Ee):
        return matrix_element(
            m_p * Ee,
            (m_n**2 - m_p**2 - m_e**2 + 2*m_p*Ee) / 2.,
            (m_p**2 + 2*m_p*Enu - m_n**2 - m_e**2) / 2.,
            m_p * Enu,
            (m_e**2 - (m_n**2 + m_p**2 - 2*m_p*(m_p + Enu - Ee))) / 2.
        ) / (32. * numpy.pi * m_p * Enu**2)

    def bounds_e_lowE(Enu):
        return ((((m_n**2-m_p**2-m_e**2) / (2. * numpy.sqrt(m_p**2 + 2.*m_p*Enu)))**2 - ((m_p*Enu)/numpy.sqrt(m_p**2+2.*m_p*Enu)+numpy.sqrt(((m_p**2+2*m_p*Enu+m_e**2-m_n**2)/(2.*numpy.sqrt(m_p**2+2*m_p*Enu)))**2-m_e**2))**2 + m_p**2 + 2*m_p*Enu - m_n**2) / (2*m_p),
        (((m_n**2-m_p**2-m_e**2) / (2. * numpy.sqrt(m_p**2 + 2.*m_p*Enu)))**2 - ((m_p*Enu)/numpy.sqrt(m_p**2+2.*m_p*Enu)-numpy.sqrt(((m_p**2+2*m_p*Enu+m_e**2-m_n**2)/(2.*numpy.sqrt(m_p**2+2*m_p*Enu)))**2-m_e**2))**2 + m_p**2 + 2*m_p*Enu - m_n**2) / (2*m_p)
        )

    def bounds_e_highE(Enu):
        return ((-(2. * Enu**2 * m_p - m_e**2 * m_p - Enu * m_e**2 + Enu * numpy.sqrt((m_p**2 + 2. * m_p * Enu - m_e**2)**2 - 2. * (m_p**2 + 2. * m_p * Enu + m_e**2) * m_p**2 + m_p**4)) / (2. * Enu + m_p) - m_n**2 + m_p**2 + 2. * m_p * Enu) / (2. * m_p),
        (-(2. * Enu**2 * m_p - m_e**2 * m_p - Enu * m_e**2 - Enu * numpy.sqrt((m_p**2 + 2. * m_p * Enu - m_e**2)**2 - 2. * (m_p**2 + 2. * m_p * Enu + m_e**2) * m_p**2 + m_p**4)) / (2. * Enu + m_p) - m_n**2 + m_p**2 + 2. * m_p * Enu) / (2. * m_p)
        )

    def bounds_nu_lowE():
        return (((m_n + m_e)**2 - m_p**2) / (2. * m_p), nu_transition)

    def bounds_nu_highE():
        return (nu_transition, particles.neutrino.grid.MAX_MOMENTUM / a)

    def Integrand_low_energy(Enu, Ee):
        return Enu**2 * dSigmadEe_lowE(Enu, Ee) * particles.neutrino.distribution(Enu * a) * (1. - particles.electron.distribution(numpy.sqrt(Ee**2-m_e**2) * a))# * radiative_correction(Enu, numpy.sqrt(Ee**2-m_e**2))

    def Integrand_high_energy(Enu, Ee):
        return Enu**2 * dSigmadEe_highE(Enu, Ee, -1.) * particles.neutrino.distribution(Enu * a) * (1. - particles.electron.distribution(numpy.sqrt(Ee**2-m_e**2) * a))

    if particles.neutrino.grid.MAX_MOMENTUM / a <= nu_transition:
        return integrate_1D(numpy.vectorize(lambda Enu: integrate_1D(lambda Ee: Integrand_low_energy(Enu, Ee), bounds=bounds_e_lowE(Enu))[0]), bounds=bounds_nu_lowE())[0] / 2. / numpy.pi**2

    return (integrate_1D(numpy.vectorize(lambda Enu: integrate_1D(lambda Ee: Integrand_low_energy(Enu, Ee), bounds=bounds_e_lowE(Enu))[0]), bounds=bounds_nu_lowE())[0] + integrate_1D(numpy.vectorize(lambda Enu: integrate_1D(lambda Ee: Integrand_high_energy(Enu, Ee), bounds=bounds_e_highE(Enu))[0]), bounds=bounds_nu_highE()))[0] / 2. / numpy.pi**2

def ne_to_pnu():
    """ n + e ⟶  p + ν_e' """

    def dSigmadEnu(Enu, Ee):
        return matrix_element(
            (m_p**2 + m_e**2 + 2*m_n*Enu - m_n**2) / 2.,
            m_n * Enu,
            m_n * Ee,
            (m_n**2 + m_e**2 + 2*m_n*Ee - m_p**2) / 2.,
            (m_e**2 - (m_n**2 + m_p**2 - 2*m_n*(m_n - Enu + Ee))) / 2.
        ) / (32. * numpy.pi * m_n * (Ee**2-m_e**2))

    def bounds_el():
        return (0., particles.electron.grid.MAX_MOMENTUM / a)

    def bounds_nu(pe):
        return ((((m_p**2-m_n**2+m_e**2) / (2. * numpy.sqrt(m_n**2 + m_e**2 + 2.*m_n*numpy.sqrt(pe**2+m_e**2))))**2 - ((m_n*pe)/numpy.sqrt(m_n**2 + m_e**2 + 2.*m_n*numpy.sqrt(pe**2+m_e**2))+(m_n**2 + m_e**2 + 2.*m_n*numpy.sqrt(pe**2+m_e**2)-m_p**2)/(2.*numpy.sqrt(m_n**2 + m_e**2 + 2.*m_n*numpy.sqrt(pe**2+m_e**2))))**2 + m_n**2 + 2*m_n*numpy.sqrt(pe**2+m_e**2) - m_p**2) / (2*m_n),
        (((m_p**2-m_n**2+m_e**2) / (2. * numpy.sqrt(m_n**2 + m_e**2 + 2.*m_n*numpy.sqrt(pe**2+m_e**2))))**2 - ((m_n*pe)/numpy.sqrt(m_n**2 + m_e**2 + 2.*m_n*numpy.sqrt(pe**2+m_e**2))-(m_n**2 + m_e**2 + 2.*m_n*numpy.sqrt(pe**2+m_e**2)-m_p**2)/(2.*numpy.sqrt(m_n**2 + m_e**2 + 2.*m_n*numpy.sqrt(pe**2+m_e**2))))**2 + m_n**2 + 2*m_n*numpy.sqrt(pe**2+m_e**2) - m_p**2) / (2*m_n)
        )

    def Integrand(Enu, pe):
        return pe**2 * numpy.sqrt(1. - m_e**2 / (pe**2+m_e**2)) * dSigmadEnu(Enu, numpy.sqrt(pe**2+m_e**2)) * particles.electron.distribution(pe * a) * (1. - particles.neutrino.distribution(Enu * a))# * radiative_correction(Enu, pe)

    return integrate_1D(numpy.vectorize(lambda pe: integrate_1D(lambda Enu: Integrand(Enu, pe), bounds=bounds_nu(pe))[0]), bounds=bounds_el())[0] / 2. / numpy.pi**2

def nnu_to_pe():
    """ n + ν_e ⟶  p + e """

    def dSigmadEe_lowE(Enu, Ee):
        return matrix_element(
            (m_n**2 + 2*m_n*Enu - m_p**2 - m_e**2) / 2.,
            m_n*Enu,
            m_n * Ee,
            (m_p**2 + 2*m_n*Ee - m_n**2 - m_e**2) / 2.,
            (m_e**2 - (m_n**2 + m_p**2 - 2*m_n*(m_n + Enu - Ee))) / 2.
        ) / (32. * numpy.pi * m_n * Enu**2)

    def bounds_e_lowE(Enu):
        return ((((m_p**2-m_n**2-m_e**2) / (2. * numpy.sqrt(m_n**2 + 2.*m_n*Enu)))**2 - ((m_n*Enu)/numpy.sqrt(m_n**2+2.*m_n*Enu)+numpy.sqrt(((m_n**2+2*m_n*Enu+m_e**2-m_p**2)/(2.*numpy.sqrt(m_n**2+2*m_n*Enu)))**2-m_e**2))**2 + m_n**2 + 2*m_n*Enu - m_p**2) / (2*m_n),
        (((m_p**2-m_n**2-m_e**2) / (2. * numpy.sqrt(m_n**2 + 2.*m_n*Enu)))**2 - ((m_n*Enu)/numpy.sqrt(m_n**2+2.*m_n*Enu)-numpy.sqrt(((m_n**2+2*m_n*Enu+m_e**2-m_p**2)/(2.*numpy.sqrt(m_n**2+2*m_n*Enu)))**2-m_e**2))**2 + m_n**2 + 2*m_n*Enu - m_p**2) / (2*m_n)
        )

    def bounds_e_highE(Enu):
        return ((-(2. * Enu**2 * m_p - m_e**2 * m_p - Enu * m_e**2 + Enu * numpy.sqrt((m_p**2 + 2. * m_p * Enu - m_e**2)**2 - 2. * (m_p**2 + 2. * m_p * Enu + m_e**2) * m_p**2 + m_p**4)) / (2. * Enu + m_p) - m_n**2 + m_p**2 + 2. * m_p * Enu) / (2. * m_p),
        (-(2. * Enu**2 * m_p - m_e**2 * m_p - Enu * m_e**2 - Enu * numpy.sqrt((m_p**2 + 2. * m_p * Enu - m_e**2)**2 - 2. * (m_p**2 + 2. * m_p * Enu + m_e**2) * m_p**2 + m_p**4)) / (2. * Enu + m_p) - m_n**2 + m_p**2 + 2. * m_p * Enu) / (2. * m_p)
        )

    def bounds_nu_lowE():
        return (0., nu_transition)

    def bounds_nu_highE():
        return (nu_transition, particles.neutrino.grid.MAX_MOMENTUM / a)

    def Integrand_low_energy(Enu, Ee):
        return Enu**2 * dSigmadEe_lowE(Enu, Ee) * particles.neutrino.distribution(Enu * a) * (1. - particles.electron.distribution(numpy.sqrt(Ee**2-m_e**2) * a)) * coulomb_correction(numpy.sqrt(Ee**2-m_e**2))# * radiative_correction(Enu, numpy.sqrt(Ee**2-m_e**2))

    def Integrand_high_energy(Enu, Ee):
        return Enu**2 * dSigmadEe_highE(Enu, Ee, 1.) * particles.neutrino.distribution(Enu * a) * (1. - particles.electron.distribution(numpy.sqrt(Ee**2 - m_e**2) * a)) * coulomb_correction(numpy.sqrt(Ee**2 - m_e**2))

    if particles.neutrino.grid.MAX_MOMENTUM / a <= nu_transition:
        return integrate_1D(numpy.vectorize(lambda Enu: integrate_1D(lambda Ee: Integrand_low_energy(Enu, Ee), bounds=bounds_e_lowE(Enu))[0]), bounds=bounds_nu_lowE())[0] / 2. / numpy.pi**2

    return (integrate_1D(numpy.vectorize(lambda Enu: integrate_1D(lambda Ee: Integrand_low_energy(Enu, Ee), bounds=bounds_e_lowE(Enu))[0]), bounds=bounds_nu_lowE())[0] + integrate_1D(numpy.vectorize(lambda Enu: integrate_1D(lambda Ee: Integrand_high_energy(Enu, Ee), bounds=bounds_e_highE(Enu))[0]), bounds=bounds_nu_highE())[0]) / 2. / numpy.pi**2

def pe_to_nnu():
    """ p + e ⟶  n + ν_e """

    def dSigmadEnu(Enu, Ee):
        return matrix_element(
            m_p * Ee,
            (m_p**2 + m_e**2 + 2*m_p*Ee - m_n**2) / 2.,
            (m_n**2 + m_e**2 + 2*m_p*Enu - m_p**2) / 2.,
            m_p * Enu,
            (m_e**2 - (m_n**2 + m_p**2 - 2*m_p*(m_p - Enu + Ee))) / 2.
        ) / (32. * numpy.pi * m_p * (Ee**2-m_e**2))

    def bounds_el():
        return (numpy.sqrt(((m_n**2 - m_e**2 - m_p**2) / (2. * m_p))**2 - m_e**2), particles.electron.grid.MAX_MOMENTUM / a)

    def bounds_nu(pe):
        return ((((m_n**2-m_p**2+m_e**2) / (2. * numpy.sqrt(m_p**2 + m_e**2 + 2.*m_p*numpy.sqrt(pe**2+m_e**2))))**2 - ((m_p*pe)/numpy.sqrt(m_p**2 + m_e**2 + 2.*m_p*numpy.sqrt(pe**2+m_e**2))+(m_p**2 + m_e**2 + 2.*m_p*numpy.sqrt(pe**2+m_e**2)-m_n**2)/(2.*numpy.sqrt(m_p**2 + m_e**2 + 2.*m_p*numpy.sqrt(pe**2+m_e**2))))**2 + m_p**2 + 2*m_p*numpy.sqrt(pe**2+m_e**2) - m_n**2) / (2*m_p),
        (((m_n**2-m_p**2+m_e**2) / (2. * numpy.sqrt(m_p**2 + m_e**2 + 2.*m_p*numpy.sqrt(pe**2+m_e**2))))**2 - ((m_p*pe)/numpy.sqrt(m_p**2 + m_e**2 + 2.*m_p*numpy.sqrt(pe**2+m_e**2))-(m_p**2 + m_e**2 + 2.*m_p*numpy.sqrt(pe**2+m_e**2)-m_n**2)/(2.*numpy.sqrt(m_p**2 + m_e**2 + 2.*m_p*numpy.sqrt(pe**2+m_e**2))))**2 + m_p**2 + 2*m_p*numpy.sqrt(pe**2+m_e**2) - m_n**2) / (2*m_p)
        )

    def Integrand(Enu, pe):
        return pe**2 * numpy.sqrt(1. - m_e**2 / (pe**2+m_e**2)) * dSigmadEnu(Enu, numpy.sqrt(pe**2+m_e**2)) * particles.electron.distribution(pe * a) * (1. - particles.neutrino.distribution(Enu * a)) * coulomb_correction(pe)# * radiative_correction(Enu, pe)

    if particles.electron.grid.MAX_MOMENTUM / a <= numpy.sqrt(((m_n**2 - m_e**2 - m_p**2) / (2. * m_p))**2 - m_e**2):
        return 0.

    return integrate_1D(numpy.vectorize(lambda pe: integrate_1D(lambda Enu: Integrand(Enu, pe), bounds=bounds_nu(pe))[0]), bounds=bounds_el())[0] / 2. / numpy.pi**2

def baryonic_rates(_a):
    global a
    a = _a

    data = numpy.array([
        nnu_to_pe(),
        pe_to_nnu(),
        n_to_penu(),
        penu_to_n(),
        ne_to_pnu(),
        pnu_to_ne()
    ])

    data /= CONST.n_decay_rate_vac

    return data

heading = [
    ["t", 's', UNITS.s],
    ["x", 'MeV', UNITS.MeV],
    ["Tg", '10^9K', UNITS.K9],
    ["dTg/dt", '10^9K/s', UNITS.K9 / UNITS.s],
    ["rho_tot", 'g/cm^3', UNITS.g_cm3],
    ["H", '1/s', 1 / UNITS.s],
    ["n nue->p e", 'dimensionless', 1.],
    ["p e->n nue", 'dimensionless', 1.],
    ["n->p e nue", 'dimensionless', 1.],
    ["p e nue->n", 'dimensionless', 1.],
    ["n e->p nue", 'dimensionless', 1.],
    ["p nue->n e", 'dimensionless', 1.]
]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run KAWANO program for the given input file')
    parser.add_argument('--folder', required=True)
    parser.add_argument('--input', default='s4.dat')
    parser.add_argument('--output', default='kawano_output.dat')
    args = parser.parse_args()
    print(run(args.folder, input=args.input, output=args.output))
