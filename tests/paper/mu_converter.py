import numpy as np
import argparse


def Gamma_nu_and_e(mass, theta):
    x = 0.511 / mass
    LL = np.log((1 - 3 * x**2 - (1 - x**2) * np.sqrt(1 - 4 * x**2))/(x**2 * (1 + np.sqrt(1 - 4 * x**2))))
    if not np.isfinite(LL):
        LL = 0.
    Gamma_ee = (
        (1.166e-11 * theta)**2 * mass**5 / (192 * np.pi**3)
        * (
            0.25 * (1 - 4 * 0.2312 + 8 * (0.2312)**2) * ((1 - 14 * x**2 - 2 * x**4 - 12 * x**6) * np.sqrt(1 - 4 * x**2) + 12 * x**4 * (x**4 - 1) * LL)
            + 4 * 0.5 * 0.2312 * (2 * 0.2312 - 1) * (x**2 * (2 + 10 * x**2 - 12 * x**4) * np.sqrt(1 - 4 * x**2) + 6 * x**4 * (1 - 2 * x**2 + 2 * x**4) * LL)
        )
    )
    return ((1.166e-11 * theta)**2 * mass**5) / (192 * np.pi**3) + Gamma_ee

# def Gamma_nu_and_e(mass, theta):
#     return ((1.166e-11 * theta)**2 * mass**5 \
#             * (1 + 0.25 * (1 - 4 * 0.2312 + 8 * (0.2312)**2))) / (192 * np.pi**3)


def Gamma_muon(mass, theta):
    x = 105.658 / mass
    return ((1.166e-11 * theta)**2 * mass**5
            * (1 - 8 * x**2 + 8 * x**6 - x**8 - 12 * x**4 * np.log(x**2))) / (192 * np.pi**3)


def Gamma_pi0(mass, theta):
    x = 134.98 / mass
    return ((1.166e-11 * theta * 130.2)**2 * mass**3 * (1 - x**2)**2) / (32 * np.pi)


def Gamma_muon_muon(mass, theta):
    x = 105.658 / mass
    LL = np.log((1 - 3 * x**2 - (1 - x**2) * np.sqrt(1 - 4 * x**2))/(x**2 * (1 + np.sqrt(1 - 4 * x**2))))
    Gamma = (
        (1.166e-11 * theta)**2 * mass**5 / (192 * np.pi**3)
        * (
            0.25 * (1 + 4 * 0.2312 + 8 * (0.2312)**2) * ((1 - 14 * x**2 - 2 * x**4 - 12 * x**6) * np.sqrt(1 - 4 * x**2) + 12 * x**4 * (x**4 - 1) * LL)
            + 4 * 0.5 * 0.2312 * (2 * 0.2312 + 1) * (x**2 * (2 + 10 * x**2 - 12 * x**4) * np.sqrt(1 - 4 * x**2) + 6 * x**4 * (1 - 2 * x**2 + 2 * x**4) * LL)
        )
    )
    return Gamma


def Gamma_pi_charged(mass, theta):
    xpi = 139.57 / mass
    xmu = 105.658 / mass
    return (
        (1.166e-11 * theta * 130.2 * 0.974)**2 * mass**3
        * ((1 - xmu**2)**2 - xpi**2 * (1 + xmu**2))
        * np.sqrt(1 + xmu**4 + xpi**4 - 2 * xmu**2 - 2 * xpi**2 - 2 * xmu**2 * xpi**2)
    ) / (16 * np.pi)


def Gamma_eta(mass, theta):
    x = 547.86 / mass
    return ((1.166e-11 * theta * 81.7)**2 * mass**3 * (1 - x**2)**2) / (32 * np.pi)


def Gamma_neutral_rho(mass, theta):
    x = 775.49 / mass
    Gamma = (1.166e-11 * theta * 208.9 * (1 - 2 * 0.2312))**2 * mass**3 * (1 + 2 * x**2) * (1 - x**2)**2 / (32 * np.pi)
    return Gamma


def Gamma_omega(mass, theta):
    x = 782.65 / mass
    Gamma = (1.166e-11 * theta * 195.5 * 4 * 0.2312 / 3)**2 * mass**3 * (1 + 2 * x**2) * (1 - x**2)**2 / (32 * np.pi)
    return Gamma


def Gamma_charged_rho(mass, theta):
    xrho = 775.11 / mass
    xmu = 105.658 / mass
    Gamma = (
        (1.166e-11 * theta * 0.974 * 209)**2 * mass**3 * ((1 - xmu**2)**2 + xrho**2 * (1 + xmu**2) -2 * xrho**4)
        * np.sqrt(1 + xmu**4 + xrho**4 - 2 * xmu**2 - 2 * xrho**2 - 2 * xmu**2 * xrho**2)
    ) / (16 * np.pi)
    return Gamma


def Gamma_eta_prime(mass, theta):
    x = 957.78 / mass
    Gamma = ((1.166e-11 * theta * 94.7)**2 * mass**3 * (1 - x**2)**2) / (32 * np.pi)
    return Gamma


def lifetime(mass, theta):
    Gamma = Gamma_nu_and_e(mass, theta)
    if mass >= 105.658 + 0.511:
        Gamma += Gamma_muon(mass, theta)
    if mass >= 134.98:
        Gamma += Gamma_pi0(mass, theta)
    if mass >= 2 * 105.658:
        Gamma += Gamma_muon_muon(mass, theta)
    if mass >= 139.57 + 105.658:
        Gamma += Gamma_pi_charged(mass, theta)
    if mass >= 547.86:
        Gamma += Gamma_eta(mass, theta)
    if mass >= 775.49:
        Gamma += Gamma_neutral_rho(mass, theta)
    if mass >= 782.65:
        Gamma += Gamma_omega(mass, theta)
    if mass >= 775.11 + 105.658:
        Gamma += Gamma_charged_rho(mass, theta)
    if mass >= 957.78:
        Gamma += Gamma_eta_prime(mass, theta)
    if mass >= 1019.46:
        raise Exception("Masses above 1019.46 MeV are not supported")

    return 6.58278e-22 / Gamma


def mixing_angle(mass, tau):
    Gamma = Gamma_nu_and_e(mass, 1.)
    if mass >= 105.658 + 0.511:
        Gamma += Gamma_muon(mass, 1)
    if mass >= 134.98:
        Gamma += Gamma_pi0(mass, 1)
    if mass >= 2 * 105.658:
        Gamma += Gamma_muon_muon(mass, 1)
    if mass >= 139.57 + 105.658:
        Gamma += Gamma_pi_charged(mass, 1)
    if mass >= 547.86:
        Gamma += Gamma_eta(mass, 1)
    if mass >= 775.49:
        Gamma += Gamma_neutral_rho(mass, 1)
    if mass >= 782.65:
        Gamma += Gamma_omega(mass, 1)
    if mass >= 775.11 + 105.658:
        Gamma += Gamma_charged_rho(mass, 1)
    if mass >= 957.78:
        Gamma += Gamma_eta_prime(mass, 1)
    if mass >= 1019.46:
        raise Exception("Masses above 1019.46 MeV are not supported")

    return np.sqrt(6.58278e-22 / Gamma / tau)


parser = argparse.ArgumentParser(description='Return mixing angle for given mass [MeV] and lifetime [s]')
parser.add_argument('--mass', required=True)
parser.add_argument('--tau', required=True)
args = parser.parse_args()

print(mixing_angle(float(args.mass), float(args.tau)))

