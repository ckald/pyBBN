# -*- coding: utf-8 -*-
"""
# Standard Model particles and interactions
"""

from __future__ import division

import itertools

from common import UNITS, CONST
from particles import STATISTICS
from interactions import Interaction
from interactions.four_particle import FourParticleM


class WeakM(FourParticleM):

    """ ## Weak interactions matrix element
        Weak processes usually include a common factor of $32 G_F^2$ """

    def __init__(self, *args, **kwargs):
        super(WeakM, self).__init__(*args, **kwargs)

        self.K1 *= 32 * CONST.G_F**2
        self.K2 *= 32 * CONST.G_F**2


class particles(object):

    """ A collection of Standard Model particles templates to be used as a reference and \
        to avoid typical mistakes such as wrong degree of freedom for the neutrinos \
        (2, not 4 - there are no right-handed neutrinos and left-handed antineutrinos). """

    photon = {
        'name': 'Photon',
        'symbol': 'γ',
        'statistics': STATISTICS.BOSON,
        'dof': 2,
        'majorana': True
    }

    class hadrons(object):
        neutron = {
            'name': 'Neutron',
            'symbol': 'n',
            'statistics': STATISTICS.FERMION,
            'mass': 0.939 * UNITS.GeV,
            'dof': 4,
            'majorana': False
        }
        proton = {
            'name': 'Proton',
            'symbol': 'p',
            'statistics': STATISTICS.FERMION,
            'mass': 0.938 * UNITS.GeV,
            'dof': 4,
            'majorana': False
        }

    class leptons(object):
        neutrino_e = {
            'name': 'Electron neutrino',
            'symbol': 'ν_e',
            'statistics': STATISTICS.FERMION,
            'mass': 0,
            'dof': 2,
            'decoupling_temperature': 5 * UNITS.MeV,
            'majorana': False,
            'flavour': 'electron'
        }
        neutrino_mu = {
            'name': 'Muon neutrino',
            'symbol': 'ν_μ',
            'statistics': STATISTICS.FERMION,
            'mass': 0,
            'dof': 2,
            'decoupling_temperature': 5 * UNITS.MeV,
            'majorana': False,
            'flavour': 'muon'
        }
        neutrino_tau = {
            'name': 'Tau neutrino',
            'symbol': 'ν_τ',
            'statistics': STATISTICS.FERMION,
            'mass': 0,
            'dof': 2,
            'decoupling_temperature': 5 * UNITS.MeV,
            'majorana': False,
            'flavour': 'tau'
        }

        electron = {
            'name': 'Electron',
            'symbol': 'e',
            'statistics': STATISTICS.FERMION,
            'mass': 0.511 * UNITS.MeV,
            'dof': 4,
            'majorana': False,
            'flavour': 'electron'
        }
        muon = {
            'name': 'Muon',
            'symbol': 'μ',
            'statistics': STATISTICS.FERMION,
            'mass': 105.7 * UNITS.MeV,
            'dof': 4,
            'majorana': False,
            'flavour': 'muon'
        }
        tau = {
            'name': 'Tau',
            'symbol': 'τ',
            'statistics': STATISTICS.FERMION,
            'mass': 1777 * UNITS.MeV,
            'dof': 4,
            'majorana': False,
            'flavour': 'tau'
        }

    class quarks(object):

        CKM = {
            (1, 1): 0.974,
            (1, 2): 0.225,
            (1, 3): 0.003,
            (2, 1): 0.225,
            (2, 2): 0.973,
            (2, 3): 0.041,
            (3, 1): 0.009,
            (3, 2): 0.040,
            (3, 3): 0.999,
        }

        up = {
            'name': 'Up quark',
            'symbol': 'u',
            'statistics': STATISTICS.FERMION,
            'mass': 3 * UNITS.MeV,
            'dof': 2,
            'majorana': False,
            'family': 1,
            'Q': 2./3.
        }
        down = {
            'name': 'Down quark',
            'symbol': 'd',
            'statistics': STATISTICS.FERMION,
            'mass': 5 * UNITS.MeV,
            'dof': 2,
            'majorana': False,
            'family': 1,
            'Q': -1./3.
        }
        charm = {
            'name': 'Charm quark',
            'symbol': 'c',
            'statistics': STATISTICS.FERMION,
            'mass': 1300 * UNITS.MeV,
            'dof': 2,
            'majorana': False,
            'family': 2,
            'Q': 2./3.
        }
        strange = {
            'name': 'Strange quark',
            'symbol': 's',
            'statistics': STATISTICS.FERMION,
            'mass': 200 * UNITS.MeV,
            'dof': 2,
            'majorana': False,
            'family': 2,
            'Q': -1./3.
        }
        top = {
            'name': 'Top quark',
            'symbol': 't',
            'statistics': STATISTICS.FERMION,
            'mass': 180 * UNITS.GeV,
            'dof': 2,
            'majorana': False,
            'family': 3,
            'Q': 2./3.
        }
        bottom = {
            'name': 'Bottom quark',
            'symbol': 'b',
            'statistics': STATISTICS.FERMION,
            'mass': 4.3 * UNITS.GeV,
            'dof': 2,
            'majorana': False,
            'family': 3,
            'Q': -1./3.
        }


class interactions(object):
    """ Collection of Standard Model interactions double-triple-checked and used \
        in all tests for consistency. """

    @staticmethod
    def neutrino_scattering(neutrino_a, neutrino_b):
        """ \begin{align}
                \nu_\alpha + \nu_\beta &\to \nu_\alpha + \nu_\beta
            \end{align}

            \begin{equation}
                |\mathcal{M}|^2 = 32 G_F^2 (p_0 \cdot p_1) (p_2 \cdot p_3)
            \end{equation}
        """

        return Interaction(
            name="Neutrino species scattering",
            particles=((neutrino_a, neutrino_b), (neutrino_a, neutrino_b)),
            antiparticles=((False, False), (False, False)),
            decoupling_temperature=0 * UNITS.MeV,
            Ms=(WeakM(K1=1., order=(0, 1, 2, 3)),)
        )

    @staticmethod
    def neutrinos_to_leptons(neutrino=None, lepton=None, g_L=CONST.g_R + 0.5):
        """ \begin{align}
                \nu_\alpha + \overline{\nu_\alpha} &\to e^- + e^+
            \end{align}

            \begin{align}
                |\mathcal{M}|^2 = 32 G_F^2 \left(
                \\\\ 4 \, g_L^2 \, (p_0 \cdot p_3) (p_1 \cdot p_2) +
                \\\\ +4 \, g_R^2 \, (p_0 \cdot p_2) (p_1 \cdot p_3) +
                \\\\ +4 \, g_L g_R \, m_2 m_3 (p_0 \cdot p_1)
                \\\\ \right)
            \end{align}

            Depending of the neutrino generations involved, $g_L$ can either be equal to \
            $g_R + \frac12$ (for $\nu_e$) or $g_R - \frac12$ (for others).
        """
        return Interaction(
            name="Neutrino pair annihilation into lepton pair",
            particles=((neutrino, neutrino), (lepton, lepton)),
            antiparticles=((False, True), (False, True)),
            decoupling_temperature=0 * UNITS.MeV,
            Ms=(
                WeakM(K1=4 * g_L**2, order=(0, 3, 1, 2)),
                WeakM(K1=4 * CONST.g_R**2, order=(0, 2, 1, 3)),
                WeakM(K2=4 * g_L * CONST.g_R, order=(2, 3, 0, 1)),
            )
        )

    @classmethod
    def neutrino_interactions(cls, leptons=None, neutrinos=None):

        g_R = CONST.g_R
        inters = []

        # Neutrinos scatterings
        for neutrino_a, neutrino_b in itertools.combinations_with_replacement(neutrinos, 2):
            inters.append(cls.neutrino_scattering(neutrino_a, neutrino_b))

        # Interactions of neutrinos and leptons
        for lepton in leptons:
            for neutrino in neutrinos:
                g_L = g_R + 0.5 if lepton.flavour == neutrino.flavour else g_R - 0.5
                inters.append(cls.neutrinos_to_leptons(g_L=g_L, lepton=lepton, neutrino=neutrino))

        return inters