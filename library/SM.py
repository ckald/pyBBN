# -*- coding: utf-8 -*-
"""
# Standard Model particles and interactions
"""

from __future__ import division

import numpy as np
import itertools
from math import sin, cos

from common import UNITS, CONST, statistics as STATISTICS
from interactions import CrossGeneratingInteraction
from interactions.four_particle import FourParticleM, FourParticleIntegral


class WeakM(FourParticleM):

    """ ## Weak interactions matrix element
        Weak processes usually include a common factor of $32 G_F^2$ """

    def __init__(self, *args, **kwargs):
        super(WeakM, self).__init__(*args, **kwargs)

        self.const = 32 * CONST.G_F**2
        self.K1 *= self.const
        self.K2 *= self.const

    def __str__(self):
        """ String-like representation of the matrix element """
        ret = ""
        if self.K1:
            ret += "K1=32 G_F^2 {: .2e} ".format(self.K1 / self.const)
        if self.K2:
            ret += "K2=32 G_F^2 {: .2e} ".format(self.K2 / self.const)
        return ret + self.order_format()


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
            'mass': 939.565378 * UNITS.MeV,
            'dof': 4,
            'majorana': False
        }
        proton = {
            'name': 'Proton',
            'symbol': 'p',
            'statistics': STATISTICS.FERMION,
            'mass': 938.272046 * UNITS.MeV,
            'dof': 4,
            'majorana': False
        }
        neutral_pion = {
            'name': 'Neutral pion',
            'symbol': 'π0',
            'statistics': STATISTICS.BOSON,
            'mass': 134.98 * UNITS.MeV,
            'dof': 1,
            'majorana': True,
            'Q': 0,
            'decay_constant': 130. / np.sqrt(2) * UNITS.MeV,
            'type': 'scalar'
        }
        charged_pion = {
            'name': 'Charged pion',
            'symbol': 'π',
            'statistics': STATISTICS.BOSON,
            'mass': 139.57 * UNITS.MeV,
            'dof': 2,
            'majorana': False,
            'Q': -1,
            'decay_constant': 130. * UNITS.MeV,
            'type': 'scalar'
        }
        neutral_rho = {
            'name': 'Neutral rho',
            'symbol': 'ρ0',
            'statistics': STATISTICS.BOSON,
            'mass': 775.49 * UNITS.MeV,
            'dof': 3,
            'majorana': True,
            'Q': 0,
            'decay_constant': 209. * UNITS.MeV,
            'type': 'vector'
        }
        charged_rho = {
            'name': 'Charged rho',
            'symbol': 'ρ',
            'statistics': STATISTICS.BOSON,
            'mass': 775.11 * UNITS.MeV,
            'dof': 3,
            'majorana': False,
            'Q': -1,
            'decay_constant': 209. * UNITS.MeV,
            'type': 'vector'
        }
        eta = {
            'name': 'Eta',
            'symbol': 'η',
            'statistics': STATISTICS.BOSON,
            'mass': 547.86 * UNITS.MeV,
            'dof': 1,
            'majorana': True,
            'Q': 0,
            'decay_constant': 156. * UNITS.MeV,
            'type': 'scalar'
        }
        eta_prime = {
            'name': 'Eta prime',
            'symbol': 'η*',
            'statistics': STATISTICS.BOSON,
            'mass': 957.78 * UNITS.MeV,
            'dof': 1,
            'majorana': True,
            'Q': 0,
            'decay_constant': 152. * UNITS.MeV,
            'type': 'scalar'
        }
        omega = {
            'name': 'Omega',
            'symbol': 'ω',
            'statistics': STATISTICS.BOSON,
            'mass': 782.65 * UNITS.MeV,
            'dof': 3,
            'majorana': True,
            'Q': 0,
            'decay_constant': 195. * UNITS.MeV,
            'type': 'vector'
        }
        phi = {
            'name': 'Phi',
            'symbol': 'ϕ',
            'statistics': STATISTICS.BOSON,
            'mass': 1019.46 * UNITS.MeV,
            'dof': 3,
            'majorana': True,
            'Q': 0,
            'decay_constant': 229. * UNITS.MeV,
            'type': 'vector'
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

        @staticmethod
        def oscillations_map(angles=(0.5905, 0.805404, 0.152346)):

            oscillations = {
                ('electron', 'electron'):
                    1 - .5 * (sin(2*angles[2])**2 + cos(angles[2])**4 * sin(2*angles[0])**2),
                ('electron', 'muon'):
                    .5 * cos(angles[2])**2 * sin(2*angles[0])**2,
                ('electron', 'tau'):
                    sin(angles[2])**2 * cos(angles[2])**2 * (2 - .5 * sin(2*angles[0])**2),
                ('muon', 'muon'):
                    1 - .5 * sin(2*angles[0])**2,
                ('muon', 'tau'):
                    .5 * sin(angles[2])**2 * sin(2*angles[0])**2,
                ('tau', 'tau'):
                    1 - sin(angles[2])**2 * (2 * cos(angles[2])**2 +
                                             .5 * sin(angles[2])**2 * sin(2*angles[0])**2)
            }
            oscillations[('muon', 'electron')] = oscillations[('electron', 'muon')]
            oscillations[('tau', 'electron')] = oscillations[('electron', 'tau')]
            oscillations[('tau', 'muon')] = oscillations[('muon', 'tau')]
            return oscillations

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

        return CrossGeneratingInteraction(
            name="Neutrino species scattering",
            particles=((neutrino_a, neutrino_b), (neutrino_a, neutrino_b)),
            antiparticles=((False, False), (False, False)),
            washout_temperature=0 * UNITS.MeV,
            Ms=(WeakM(K1=1., order=(0, 1, 2, 3)),),
            integral_type=FourParticleIntegral
        )

    @staticmethod
    def neutrinos_to_leptons(neutrino=None, lepton=None, g_L=CONST.g_R + 0.5):
        """ \begin{align}
                \nu_\alpha + \overline{\nu_\alpha} &\to e^- + e^+
            \end{align}

            \begin{align}
                |\mathcal{M}|^2 = 32 G_F^2 \left(
                \\ 4 \, g_L^2 \, (p_0 \cdot p_3) (p_1 \cdot p_2) +
                \\ +4 \, g_R^2 \, (p_0 \cdot p_2) (p_1 \cdot p_3) +
                \\ +4 \, g_L g_R \, m_2 m_3 (p_0 \cdot p_1)
                \\ \right)
            \end{align}

            Depending of the neutrino generations involved, $g_L$ can either be equal to \
            $g_R + \frac12$ (for $\nu_e$) or $g_R - \frac12$ (for others).
        """
        return CrossGeneratingInteraction(
            name="Neutrino pair annihilation into lepton pair",
            particles=((neutrino, neutrino), (lepton, lepton)),
            antiparticles=((False, True), (False, True)),
            washout_temperature=0 * UNITS.MeV,
            Ms=(
                WeakM(K1=4 * g_L**2, order=(0, 3, 1, 2)),
                WeakM(K1=4 * CONST.g_R**2, order=(0, 2, 1, 3)),
                WeakM(K2=4 * g_L * CONST.g_R, order=(2, 3, 0, 1)),
            ),
            integral_type=FourParticleIntegral
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

    @classmethod
    def baryons_interaction(cls, neutron=None, proton=None, neutrino=None, electron=None):

        return CrossGeneratingInteraction(
            name="Baryons interaction",
            particles=((neutron,), (proton, electron, neutrino)),
            antiparticles=((False,), (False, False, True)),
            washout_temperature=0 * UNITS.MeV,
            Ms=(WeakM(K1=2. * particles.quarks.CKM[(1, 1)]**2, order=(0, 1, 2, 3)),),
            integral_type=FourParticleIntegral
        )
