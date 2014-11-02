from common import UNITS
from particles import STATISTICS


class StandardModelParticles:
    """ A collection of Standard Model particles templates to be used as a reference and to avoid\
        typical mistakes such as wrong degree of freedom for the neutrinos (2, not 4 - there are\
        no right-handed neutrinos and left-handed antineutrinos). """

    photon = {
        'name': 'Photon',
        'statistics': STATISTICS.BOSON,
        'dof': 2
    }
    neutrino_e = {
        'name': 'Electron neutrino',
        'statistics': STATISTICS.FERMION,
        'dof': 2,
        'decoupling_temperature': 3 * UNITS.MeV
    }
    neutrino_mu = {
        'name': 'Muon neutrino',
        'statistics': STATISTICS.FERMION,
        'dof': 2,
        'decoupling_temperature': 3 * UNITS.MeV
    }
    neutrino_tau = {
        'name': 'Tau neutrino',
        'statistics': STATISTICS.FERMION,
        'dof': 2,
        'decoupling_temperature': 3 * UNITS.MeV
    }
    neutron = {
        'name': 'Neutron',
        'statistics': STATISTICS.FERMION,
        'mass': 0.939 * UNITS.GeV,
        'dof': 4
    }
    proton = {
        'name': 'Proton',
        'statistics': STATISTICS.FERMION,
        'mass': 0.938 * UNITS.GeV,
        'dof': 4
    }
    electron = {
        'name': 'Electron',
        'mass': 0.511 * UNITS.MeV,
        'statistics': STATISTICS.FERMION,
        'dof': 4
    }
