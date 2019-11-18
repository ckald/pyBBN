# -*- coding: utf-8 -*-
"""
# Standard Model particles and interactions
"""

from __future__ import division

import numpy as np
import itertools
from collections import Counter
from math import sin, cos, atan

import environment
from common import UNITS, CONST, statistics as STATISTICS
from interactions import CrossGeneratingInteraction
from interactions.three_particle import ThreeParticleM, ThreeParticleIntegral
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
        if self.K:
            ret += "|M|²={: .4e}".format(self.K)
            return ret
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
        'majorana': True,
        'Q': 0
    }

    gluon = {
        'name': 'Gluon',
        'symbol': 'g',
        'statistics': STATISTICS.BOSON,
        'dof': 16,
        'majorana': True,
        'Q': 0
    }

    class hadrons(object):
        neutron = {
            'name': 'Neutron',
            'symbol': 'n',
            'statistics': STATISTICS.FERMION,
            'mass': 939.565378 * UNITS.MeV,
            'dof': 4,
            'majorana': False,
            'Q': 0
        }
        proton = {
            'name': 'Proton',
            'symbol': 'p',
            'statistics': STATISTICS.FERMION,
            'mass': 938.272046 * UNITS.MeV,
            'dof': 4,
            'majorana': False,
            'Q': -1
        }
        neutral_pion = {
            'name': 'Neutral pion',
            'symbol': 'π0',
            'statistics': STATISTICS.BOSON,
            'mass': 134.98 * UNITS.MeV,
            'dof': 1,
            'decoupling_temperature': 6.8 * UNITS.MeV, # T @ Gamma_decay > Gamma_annihilation
            'majorana': True,
            'Q': 0,
            'decay_constant': 130.2 * UNITS.MeV,
            'type': 'scalar',
            'fast_decay': True,
            'lifetime': 8.52e-17 * UNITS.s,
            'BR': {'γγ': 0.98823}
        }
        charged_pion = {
            'name': 'Charged pion',
            'symbol': 'π',
            'statistics': STATISTICS.BOSON,
            'mass': 139.57 * UNITS.MeV,
            'dof': 2,
            'decoupling_temperature': 6.8 * UNITS.MeV, # T @ Gamma_decay > Gamma_annihilation
            'majorana': False,
            'Q': -1,
            'decay_constant': 130.2 * UNITS.MeV,
            'type': 'scalar',
            'fast_decay': True,
            'thermalization': True,
            'lifetime': 2.603e-8 * UNITS.s,
            'BR': {'μν_μ': 0.999877}
        }
        charged_kaon = {
            'name': 'Charged kaon',
            'symbol': 'K',
            'statistics': STATISTICS.BOSON,
            'mass': 493.677 * UNITS.MeV,
            'dof': 2,
            'decoupling_temperature': 25.1 * UNITS.MeV, # T @ Gamma_decay > Gamma_annihilation
            'majorana': False,
            'Q': -1,
            'decay_constant': 155.6 * UNITS.MeV,
            'type': 'scalar',
            'fast_decay': True,
            'thermalization': True, #TODO: Is this true?
            'lifetime': 1.238e-8 * UNITS.s,
            'BR': {
                'π0eν_e': 0.0507,
                'πππ': 0.05583,
                'μν_μ': 0.6356,
                'ππ0': 0.2067
            }
        }
        kaon_long = {
            'name': 'Kaon long',
            'symbol': 'K0L',
            'statistics': STATISTICS.BOSON,
            'mass': 497.611 * UNITS.MeV,
            'dof': 1,
            'majorana': True,
            'Q': 0,
            'decay_constant': 155.6 * UNITS.MeV, #TODO: Is this true?
            'type': 'scalar',
            'fast_decay': True,
            'lifetime': 5.116e-8 * UNITS.s,
            'BR': {
                'πeν_e': 0.4055,
                'πμν_μ': 0.2704,
                'π0π0π0': 0.1952,
                'πππ0': 0.1254
            }
        }
        kaon_short = {
            'name': 'Kaon short',
            'symbol': 'K0S',
            'statistics': STATISTICS.BOSON,
            'mass': 497.611 * UNITS.MeV,
            'dof': 1,
            'majorana': True,
            'Q': 0,
            'decay_constant': 155.6 * UNITS.MeV, #TODO: Is this true?
            'type': 'scalar',
            'fast_decay': True,
            'lifetime': 8.954e-11 * UNITS.s,
            'BR': {
                'ππ': 0.6920,
                'π0π0': 0.3069
            }
        }
        eta = {
            'name': 'Eta',
            'symbol': 'η',
            'statistics': STATISTICS.BOSON,
            'mass': 547.86 * UNITS.MeV,
            'dof': 1,
            'majorana': True,
            'Q': 0,
            'decay_constant': 81.7 * UNITS.MeV,
            'type': 'scalar',
            'fast_decay': True,
            'lifetime': 5.025e-19 * UNITS.s,
            'BR': {
                'γγ': 0.3941,
                'π0π0π0': 0.3268,
                'πππ0': 0.2292,
                'ππγ': 0.0422
            }
        }
        charged_rho = {
            'name': 'Charged rho',
            'symbol': 'ρ',
            'statistics': STATISTICS.BOSON,
            'mass': 775.11 * UNITS.MeV,
            'dof': 6,
            'majorana': False,
            'Q': -1,
            'decay_constant': 209. * UNITS.MeV,
            'type': 'vector',
            'fast_decay': True,
            'lifetime': 4.415e-24 * UNITS.s,
            'BR': {'ππ0': 1.}
        }
        neutral_rho = {
            'name': 'Neutral rho',
            'symbol': 'ρ0',
            'statistics': STATISTICS.BOSON,
            'mass': 775.49 * UNITS.MeV,
            'dof': 3,
            'majorana': True,
            'Q': 0,
            'decay_constant': 208.9 * UNITS.MeV,
            'type': 'vector',
            'fast_decay': True,
            'lifetime': 4.415e-24 * UNITS.s,
            'BR': {'ππ': 1.}
        }
        omega = {
            'name': 'Omega',
            'symbol': 'ω',
            'statistics': STATISTICS.BOSON,
            'mass': 782.65 * UNITS.MeV,
            'dof': 3,
            'majorana': True,
            'Q': 0,
            'decay_constant': 195.5 * UNITS.MeV,
            'type': 'vector',
            'fast_decay': True,
            'lifetime': 7.753e-23 * UNITS.s,
            'BR': {
                'πππ0': 0.892,
                'π0γ': 0.0840
            }
        }
        eta_prime = {
            'name': 'Eta prime',
            'symbol': 'η*',
            'statistics': STATISTICS.BOSON,
            'mass': 957.78 * UNITS.MeV,
            'dof': 1,
            'majorana': True,
            'Q': 0,
            'decay_constant': -94.7 * UNITS.MeV,
            'type': 'scalar',
            'fast_decay': True,
            'lifetime': 3.358e-21 * UNITS.s,
            'BR': {
                'ππη': 0.426,
                'ρ0γ': 0.289,
                'π0π0η': 0.228
            }
        }
        phi = {
            'name': 'Phi',
            'symbol': 'ϕ',
            'statistics': STATISTICS.BOSON,
            'mass': 1019.46 * UNITS.MeV,
            'dof': 3,
            'majorana': True,
            'Q': 0,
            'decay_constant': 229.5 * UNITS.MeV,
            'type': 'vector',
            'fast_decay': True,
            'lifetime': 1.549e-22 * UNITS.s,
            'BR': {
                'KK': 0.489,
                'K0LK0S': 0.342,
                'ρ0π0': 0.1532
            }
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
            'Q': 0,
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
            'Q': 0,
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
            'Q': 0,
            'flavour': 'tau'
        }
        electron = {
            'name': 'Electron',
            'symbol': 'e',
            'statistics': STATISTICS.FERMION,
            'mass': 0.511 * UNITS.MeV,
            'dof': 4,
            'majorana': False,
            'Q': -1,
            'flavour': 'electron'
        }
        muon = {
            'name': 'Muon',
            'symbol': 'μ',
            'statistics': STATISTICS.FERMION,
            'mass': 105.658 * UNITS.MeV,
            'dof': 4,
            'decoupling_temperature': 4.6 * UNITS.MeV, # T @ Gamma_decay > Gamma_annihilation
            'majorana': False,
            'Q': -1,
            'flavour': 'muon',
            'fast_decay': True,
            'thermalization': True,
            'lifetime': 2.197e-6 * UNITS.s,
            'BR': {'eν_eν_μ': 1.}
        }
        tau = {
            'name': 'Tau',
            'symbol': 'τ',
            'statistics': STATISTICS.FERMION,
            'mass': 1776.82 * UNITS.MeV,
            'dof': 4,
            'majorana': False,
            'Q': -1,
            'flavour': 'tau',
            'fast_decay': True,
            'lifetime': 2.903e-13 * UNITS.s
        }

        @staticmethod
        def oscillations_map(MSW_12=None, MSW_13=None, matter_effects=True):

            # angles=(0.5905, 0.71011, 0.152346) # theta_12, theta_23, theta_13
            theta_12 = 0.57636
            theta_23 = 0.71011 if environment.get('NORMAL_HIERARCHY_NEUTRINOS') else 0.87487
            theta_13 = 0.14715

            if matter_effects:
                theta_12 = 0.5 * np.arctan(np.sin(2.*theta_12) / (np.cos(2.*theta_12) + MSW_12))
                theta_23 = 0.5 * np.arctan(np.sin(2.*theta_23) / (np.cos(2.*theta_23) + MSW_13)) # MSW_23 \approx MSW_13
                theta_13 = 0.5 * np.arctan(np.sin(2.*theta_13) / (np.cos(2.*theta_13) + MSW_13))

            PMNS = {
                (1, 1): np.cos(theta_12) * np.cos(theta_13),
                (1, 2): np.cos(theta_13) * np.sin(theta_12),
                (1, 3): np.sin(theta_13),
                (2, 1): -np.cos(theta_23) * np.sin(theta_12) - np.cos(theta_12) * np.sin(theta_13) * np.sin(theta_23),
                (2, 2): np.cos(theta_12) * np.cos(theta_23) - np.sin(theta_12) * np.sin(theta_13) * np.sin(theta_23),
                (2, 3): np.cos(theta_13) * np.sin(theta_23),
                (3, 1): np.sin(theta_23) * np.sin(theta_12) - np.cos(theta_12) * np.cos(theta_23) * np.sin(theta_13),
                (3, 2): -np.cos(theta_12) * np.sin(theta_23) - np.cos(theta_23) * np.sin(theta_12) * np.sin(theta_13),
                (3, 3): np.cos(theta_13) * np.cos(theta_23),
            }

            # oscillations = {
            #     ('electron', 'electron'):
            #         1 - .5 * (np.sin(2*theta_13)**2 + np.cos(theta_13)**4 * np.sin(2*theta_12)**2),
            #     ('electron', 'muon'):
            #         .5 * np.cos(theta_13)**2 * np.sin(2*theta_12)**2,
            #     ('electron', 'tau'):
            #         np.sin(theta_13)**2 * np.cos(theta_13)**2 * (2 - .5 * np.sin(2*theta_12)**2),
            #     ('muon', 'muon'):
            #         1 - .5 * np.sin(2*theta_12)**2,
            #     ('muon', 'tau'):
            #         .5 * np.sin(theta_13)**2 * np.sin(2*theta_12)**2,
            #     ('tau', 'tau'):
            #         1 - np.sin(theta_13)**2 * (2 * np.cos(theta_13)**2 +
            #                                  .5 * np.sin(theta_13)**2 * np.sin(2*theta_12)**2)
            # }

            oscillations = {
                ('electron', 'electron'):
                    1 - 2 * (PMNS[(1, 1)])**2 * (PMNS[(1, 2)])**2 - 2 * PMNS[(1, 3)]**2 * (1 - PMNS[(1, 3)]**2),
                ('electron', 'muon'):
                    -2 * PMNS[(1, 1)] * PMNS[(1, 2)] * PMNS[(2, 1)] * PMNS[(2, 2)] + 2 * PMNS[(1, 3)]**2 * PMNS[(2, 3)]**2,
                ('electron', 'tau'):
                    -2 * PMNS[(1, 1)] * PMNS[(1, 2)] * PMNS[(3, 1)] * PMNS[(3, 2)] + 2 * PMNS[(1, 3)]**2 * PMNS[(3, 3)]**2,
                ('muon', 'muon'):
                    1 - 2 * (PMNS[(2, 1)])**2 * (PMNS[(2, 2)])**2 - 2 * PMNS[(2, 3)]**2 * (1 - PMNS[(2, 3)]**2),
                ('muon', 'tau'):
                    -2 * PMNS[(2, 1)] * PMNS[(2, 2)] * PMNS[(3, 1)] * PMNS[(3, 2)] + 2 * PMNS[(2, 3)]**2 * PMNS[(3, 3)]**2,
                ('tau', 'tau'):
                    1 - 2 * (PMNS[(3, 1)])**2 * (PMNS[(3, 2)])**2 - 2 * PMNS[(3, 3)]**2 * (1 - PMNS[(3, 3)]**2)
            }

            oscillations[('muon', 'electron')] = oscillations[('electron', 'muon')]
            oscillations[('tau', 'electron')] = oscillations[('electron', 'tau')]
            oscillations[('tau', 'muon')] = oscillations[('muon', 'tau')]
            return oscillations

    class quarks(object):

        CKM = {
            (1, 1): 0.974,
            (1, 2): 0.225,
            (1, 3): 0.004,
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
            'mass': 2.2 * UNITS.MeV,
            'dof': 12,
            'majorana': False,
            'family': 1,
            'Q': 2./3.
        }
        down = {
            'name': 'Down quark',
            'symbol': 'd',
            'statistics': STATISTICS.FERMION,
            'mass': 4.7 * UNITS.MeV,
            'dof': 12,
            'majorana': False,
            'family': 1,
            'Q': -1./3.
        }
        charm = {
            'name': 'Charm quark',
            'symbol': 'c',
            'statistics': STATISTICS.FERMION,
            'mass': 1.28 * UNITS.GeV,
            'dof': 12,
            'majorana': False,
            'family': 2,
            'Q': 2./3.
        }
        strange = {
            'name': 'Strange quark',
            'symbol': 's',
            'statistics': STATISTICS.FERMION,
            'mass': 96 * UNITS.MeV,
            'dof': 12,
            'majorana': False,
            'family': 2,
            'Q': -1./3.
        }
        top = {
            'name': 'Top quark',
            'symbol': 't',
            'statistics': STATISTICS.FERMION,
            'mass': 173.1 * UNITS.GeV,
            'dof': 12,
            'majorana': False,
            'family': 3,
            'Q': 2./3.
        }
        bottom = {
            'name': 'Bottom quark',
            'symbol': 'b',
            'statistics': STATISTICS.FERMION,
            'mass': 4.18 * UNITS.GeV,
            'dof': 12,
            'majorana': False,
            'family': 3,
            'Q': -1./3.
        }


class interactions(object):
    """ Collection of Standard Model interactions double-triple-checked and used \
        in all tests for consistency. """

    @staticmethod
    def neutrino_scattering(neutrino_a, neutrino_b, kind=None):
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
            integral_type=FourParticleIntegral,
            kind=kind
        )

    @staticmethod
    def neutrinos_to_leptons(neutrino=None, lepton=None, g_L=CONST.g_R + 0.5, kind=None):
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
                WeakM(K1=4 * CONST.g_R**2, order=(0, 3, 1, 2)),
                WeakM(K1=4 * g_L**2, order=(0, 2, 1, 3)),
                WeakM(K2=4 * g_L * CONST.g_R, order=(2, 3, 0, 1)),
            ),
            integral_type=FourParticleIntegral,
            kind=kind
        )

    @staticmethod
    def leptons_CC(active1=None, active2=None, lepton1=None, lepton2=None, kind=None):

        return CrossGeneratingInteraction(
            name="Charged current lepton-neutrino interactions",
            particles=((lepton1, lepton2), (active1, active2)),
            antiparticles=((False, True), (False, True)),
            Ms=(
                WeakM(K1=4., order=(0, 3, 1, 2)),
            ),
            integral_type=FourParticleIntegral,
            kind=kind
        )

    @staticmethod
    def leptons_leptons_NC(lepton1=None, lepton2=None, g_L=CONST.g_R - 0.5, kind=None):

        return CrossGeneratingInteraction(
            name="Neutral current lepton-lepton interactions",
            particles=((lepton1, lepton2), (lepton1, lepton2)),
            antiparticles=((False, False), (False, False)),
            Ms=(
                WeakM(K1=4 * (g_L**4 + CONST.g_R**4), order=(0, 1, 2, 3)),
                WeakM(K1=4 * 2 * g_L**2 * CONST.g_R**2, order=(0, 3, 1, 2)),
                WeakM(K2=4 * -g_L * CONST.g_R * (g_L**2 + CONST.g_R**2), order=(0, 2, 1, 3)),
            ),
            integral_type=FourParticleIntegral,
            kind=kind
        )

    @classmethod
    def neutrino_interactions(cls, leptons=None, neutrinos=None, kind=None):

        g_R = CONST.g_R
        inters = []

        # Neutrinos scatterings
        for neutrino_a, neutrino_b in itertools.combinations_with_replacement(neutrinos, 2):
            inters.append(cls.neutrino_scattering(neutrino_a, neutrino_b, kind))

        # Interactions of neutrinos and leptons
        for lepton in leptons:
            for neutrino in neutrinos:
                g_L = g_R + 0.5 if lepton.flavour == neutrino.flavour else g_R - 0.5
                inters.append(cls.neutrinos_to_leptons(g_L=g_L, lepton=lepton, neutrino=neutrino, kind=kind))

        return inters

    @classmethod
    def lepton_interactions(cls, leptons=None, neutrinos=None, SM_inters=False, kind=None):

        g_R = CONST.g_R
        inters = []

        def reaction_type(reaction):
            return sum(reactant.side for reactant in reaction)

        # Interactions of neutrinos and leptons
        for lepton in leptons:
            for neutrino in neutrinos:
                # if lepton.name != 'Electron' and SM_inters == False:
                #     g_L = g_R + 0.5 if lepton.flavour == neutrino.flavour else g_R - 0.5
                #     inters.append(cls.neutrinos_to_leptons(g_L=g_L, lepton=lepton, neutrino=neutrino, kind=kind))

                if len(leptons) > 1 and len(neutrinos) > 1:
                    other_lepton = [par for par in leptons if par != lepton]
                    other_neutrino = [par for par in neutrinos if par != neutrino]
                    for other_lep in other_lepton:
                        for other_neut in other_neutrino:
                            if neutrino.flavour == lepton.flavour and other_neut.flavour == other_lep.flavour:
                                inter = cls.leptons_CC(
                                    active1=neutrino,
                                    active2=other_neut,
                                    lepton1=lepton,
                                    lepton2=other_lep,
                                    kind=kind
                                )

                                inters.append(inter)
        # if len(leptons) > 1:
        #     lepton = leptons[0]
        #     other_lepton = [par for par in leptons if par != lepton]
        #     for other_lep in other_lepton:
        #         inters.append(cls.leptons_leptons_NC(
        #             lepton1=lepton,
        #             lepton2=other_lep,
        #             kind=kind
        #         ))
        return inters

    @staticmethod
    def decay_charged_pion(meson=None, lepton=None, neutrino=None, kind=None):

        CKM = particles.quarks.CKM[(1, 1)]

        return [CrossGeneratingInteraction(
            name="Charged pion decay into neutrino and charged lepton",
            particles=((meson, ), (lepton, neutrino)),
            antiparticles=antiparticles,
            Ms=(ThreeParticleM(
                K=2 * (CONST.G_F * meson.decay_constant * CKM)**2 * lepton.mass**4 * (
                    ((meson.mass / lepton.mass)**2 - 1)
                )
            ),),
            integral_type=ThreeParticleIntegral,
            kind=kind
        ) for antiparticles in [
            ((False, ), (False, True)),
            ((True, ), (True, False))
        ]]

    @staticmethod
    def decay_neutral_pion(meson=None, photon=None, kind=None):

        return [CrossGeneratingInteraction(
            name="Neutral pion decay into two photons",
            particles=((meson, ), (photon, photon)),
            antiparticles=((False, ), (False, False)),
            Ms=(ThreeParticleM(
                K=CONST.alpha**2 * meson.mass**4 / (2 * np.pi**2 * meson.decay_constant**2)
            ),),
            integral_type=ThreeParticleIntegral,
            kind=kind
        )]

    @staticmethod
    def four_particle_meson_decay(parts=None, antiparts=None, ME=None, kind=None):

        return [CrossGeneratingInteraction(
            name="Meson decay",
            particles=((parts[0], parts[1]), (parts[2], parts[3])),
            antiparticles=((antiparts[0], antiparts[1]), (antiparts[2], antiparts[3])),
            Ms=(WeakM(K=ME),),
            integral_type=FourParticleIntegral,
            kind=kind
        )]

    @staticmethod
    def three_particle_meson_decay(parts=None, antiparts=None, ME=None, kind=None):

        return [CrossGeneratingInteraction(
            name="Meson decay",
            particles=((parts[0], ), (parts[1], parts[2])),
            antiparticles=((antiparts[0], ), (antiparts[1], antiparts[2])),
            Ms=(ThreeParticleM(K=ME),),
            integral_type=ThreeParticleIntegral,
            kind=kind
        )]

    @classmethod
    def meson_interactions(cls, primary_mesons=None, mesons=None, leptons=None, neutrinos=None,
                        photon=None, muon_tau=None, kind=None):

        inters = []

        photon = photon[0] if photon else []
        for lepton in leptons or []:
            if lepton.name == 'Electron':
                electron = lepton
            if lepton.name == 'Muon':
                muon = lepton
            if lepton.name == 'Tau':
                tau = lepton
        for neutrino in neutrinos or []:
            if neutrino.name == 'Electron neutrino':
                neutrino_e = neutrino
            if neutrino.name == 'Muon neutrino':
                neutrino_mu = neutrino
            if neutrino.name == 'Tau neutrino':
                neutrino_tau = neutrino
        for meson in mesons or []:
            if meson.name == 'Charged pion':
                charged_pion = meson
            if meson.name == 'Neutral pion':
                neutral_pion = meson
            if meson.name == 'Charged kaon':
                charged_kaon = meson
            if meson.name == 'Kaon long':
                kaon_long = meson
            if meson.name == 'Kaon short':
                kaon_short = meson
            if meson.name == 'Eta':
                eta = meson
            if meson.name == 'Neutral rho':
                neutral_rho = meson
            if meson.name == 'Charged rho':
                charged_rho = meson
            if meson.name == 'Omega':
                omega = meson
            if meson.name == 'Eta prime':
                eta_prime = meson
            if meson.name == 'Phi':
                phi = meson

        matrix_elements = {
            # Four-particle
            'EtaThreePi0': 0.0870984,
            'EtaPiPlusPiMinPi0': 0.0690629,
            'EtaPiPlusPiMinGamma': 0.0046653,
            'OmegaPiPlusPiMinPi0': 1145.69,
            'EtaPrimePiPlusPiMinEta': 43.888,
            'EtaPrimeTwoPi0Eta': 20.0986,
            'KaonMinPi0ElecNue': 1.42906e-13,
            'KaonMinPiPlusTwoPiMin': 1.85537e-12,
            'KaonLPiPlusElecNue': 2.80345e-13,
            'KaonLPiPlusMuonNumu': 3.03627e-13,
            'KaonLThreePi0': 1.05573e-12,
            'KaonLPiPlusPiMinPi0': 8.26989e-13,

            # Three-particle
            'RhoPlusPiPlusPi0': 1.8639e7 * UNITS.MeV**2,
            'Rho0PiPlusPiMin': 1.86839e7 * UNITS.MeV**2,
            'OmegaPi0Gamma': 85508.6 * UNITS.MeV**2,
            'EtaPrimeRho0Gamma': 8044.63 * UNITS.MeV**2,
            'PhiKPlusKMin': 1.28798e6 * UNITS.MeV**2,
            'PhiK0LK0S': 1.03471e6 * UNITS.MeV**2,
            'PhiRho0Pi0': 286706. * UNITS.MeV**2,
            'EtaTwoGamma': 14.2174 * UNITS.MeV**2,
            'KaonMinMuonNumu': 8.78918e-10 * UNITS.MeV**2,
            'KaonMinPiMinPi0': 3.28177e-10 * UNITS.MeV**2,
            'KaonSTwoPi0': 6.718e-8 * UNITS.MeV**2,
            'KaonSPiPlusPiMin': 1.53713e-7 * UNITS.MeV**2
        }

        def already_there(interaction, name):
            if not interaction:
                return False

            for inter in interaction:
                for integral in inter.integrals:
                    if Counter(item.specie.name for item in integral.reaction if item.specie.name == name):
                        return True

            return False

        def neutral_pion_decay(inters=None):
            inters += cls.decay_neutral_pion(meson=neutral_pion,
                                            photon=photon,
                                            kind=kind)

        def charged_pion_decay(inters=None): # Check if this only includes muon decay
            if not (muon_tau[0] or already_there(inters, 'Muon')):
                inters += cls.lepton_interactions(leptons=leptons,
                                                neutrinos=neutrinos,
                                                kind=kind)

            inters += cls.decay_charged_pion(meson=charged_pion,
                                            lepton=muon,
                                            neutrino=neutrino_mu,
                                            kind=kind)

        def charged_kaon_decay(inters=None):
            if not already_there(inters, 'Charged pion'):
                charged_pion_decay(inters)

            if not already_there(inters, 'Neutral pion'):
                neutral_pion_decay(inters)

            if not muon_tau[0] or already_there(inters, 'Muon'):
                interactions.lepton_interactions(leptons=leptons,
                                                neutrinos=neutrinos,
                                                kind=kind)

            inters += cls.four_particle_meson_decay(parts=[charged_kaon, neutral_pion, electron, neutrino_e],
                                                    antiparts=[False, False, False, True],
                                                    ME=matrix_elements['KaonMinPi0ElecNue'],
                                                    kind=kind)
            inters += cls.four_particle_meson_decay(parts=[charged_kaon, charged_pion, charged_pion, charged_pion],
                                                    antiparts=[False, True, False, False],
                                                    ME=matrix_elements['KaonMinPiPlusTwoPiMin'],
                                                    kind=kind)
            inters += cls.three_particle_meson_decay(parts=[charged_kaon, muon, neutrino_mu],
                                                    antiparts=[False, False, True],
                                                    ME=matrix_elements['KaonMinMuonNumu'],
                                                    kind=kind)
            inters += cls.three_particle_meson_decay(parts=[charged_kaon, charged_pion, neutral_pion],
                                                    antiparts=[False, False, False],
                                                    ME=matrix_elements['KaonMinPiMinPi0'],
                                                    kind=kind)

        def kaon_long_decay(inters=None):
            if not already_there(inters, 'Charged pion'):
                charged_pion_decay(inters)

            if not already_there(inters, 'Neutral pion'):
                neutral_pion_decay(inters)

            if not muon_tau[0] or already_there(inters, 'Muon'):
                interactions.lepton_interactions(leptons=leptons,
                                                neutrinos=neutrinos,
                                                kind=kind)

            inters += cls.four_particle_meson_decay(parts=[kaon_long, charged_pion, electron, neutrino_e],
                                                    antiparts=[False, True, False, True],
                                                    ME=matrix_elements['KaonLPiPlusElecNue'],
                                                    kind=kind)
            inters += cls.four_particle_meson_decay(parts=[kaon_long, charged_pion, muon, neutrino_mu],
                                                    antiparts=[False, True, False, True],
                                                    ME=matrix_elements['KaonLPiPlusMuonNumu'],
                                                    kind=kind)
            inters += cls.four_particle_meson_decay(parts=[kaon_long, neutral_pion, neutral_pion, neutral_pion],
                                                    antiparts=[False, False, False, False],
                                                    ME=matrix_elements['KaonLThreePi0'],
                                                    kind=kind)
            inters += cls.four_particle_meson_decay(parts=[kaon_long, neutral_pion, charged_pion, charged_pion],
                                                    antiparts=[False, False, False, True],
                                                    ME=matrix_elements['KaonLPiPlusPiMinPi0'],
                                                    kind=kind)

        def kaon_short_decay(inters=None):
            if not already_there(inters, 'Charged pion'):
                charged_pion_decay(inters)

            if not already_there(inters, 'Neutral pion'):
                neutral_pion_decay(inters)

            inters += cls.three_particle_meson_decay(parts=[kaon_short, neutral_pion, neutral_pion],
                                                    antiparts=[False, False, False],
                                                    ME=matrix_elements['KaonSTwoPi0'],
                                                    kind=kind)
            inters += cls.three_particle_meson_decay(parts=[kaon_short, charged_pion, charged_pion],
                                                    antiparts=[False, False, True],
                                                    ME=matrix_elements['KaonSPiPlusPiMin'],
                                                    kind=kind)

        def eta_decay(inters=None):
            if not already_there(inters, 'Charged pion'):
                charged_pion_decay(inters)

            if not already_there(inters, 'Neutral pion'):
                neutral_pion_decay(inters)

            inters += cls.three_particle_meson_decay(parts=[eta, photon, photon],
                                                    antiparts=[False, False, False],
                                                    ME=matrix_elements['EtaTwoGamma'],
                                                    kind=kind)
            inters += cls.four_particle_meson_decay(parts=[eta, neutral_pion, neutral_pion, neutral_pion],
                                                    antiparts=[False, False, False, False],
                                                    ME=matrix_elements['EtaThreePi0'],
                                                    kind=kind)
            inters += cls.four_particle_meson_decay(parts=[eta, neutral_pion, charged_pion, charged_pion],
                                                    antiparts=[False, False, False, True],
                                                    ME=matrix_elements['EtaPiPlusPiMinPi0'],
                                                    kind=kind)
            inters += cls.four_particle_meson_decay(parts=[eta, photon, charged_pion, charged_pion],
                                                    antiparts=[False, False, False, True],
                                                    ME=matrix_elements['EtaPiPlusPiMinGamma'],
                                                    kind=kind)

        def neutral_rho_decay(inters=None):
            if not already_there(inters, 'Charged pion'):
                charged_pion_decay(inters)

            inters += cls.three_particle_meson_decay(parts=[neutral_rho, charged_pion, charged_pion],
                                                    antiparts=[False, True, False],
                                                    ME=matrix_elements['Rho0PiPlusPiMin'],
                                                    kind=kind)

        def charged_rho_decay(inters=None):
            if not already_there(inters, 'Charged pion'):
                charged_pion_decay(inters)

            if not already_there(inters, 'Neutral pion'):
                neutral_pion_decay(inters)

            inters += cls.three_particle_meson_decay(parts=[charged_rho, charged_pion, neutral_pion],
                                                    antiparts=[False, False, False],
                                                    ME=matrix_elements['RhoPlusPiPlusPi0'],
                                                    kind=kind)

        def omega_decay(inters=None):
            if not already_there(inters, 'Charged pion'):
                charged_pion_decay(inters)

            if not already_there(inters, 'Neutral pion'):
                neutral_pion_decay(inters)

            inters += cls.four_particle_meson_decay(parts=[omega, charged_pion, charged_pion, neutral_pion],
                                                    antiparts=[False, True, False, False],
                                                    ME=matrix_elements['OmegaPiPlusPiMinPi0'],
                                                    kind=kind)
            inters += cls.three_particle_meson_decay(parts=[omega, neutral_pion, photon],
                                                    antiparts=[False, False, False],
                                                    ME=matrix_elements['OmegaPi0Gamma'],
                                                    kind=kind)

        def eta_prime_decay(inters=None):
            if not already_there(inters, 'Neutral rho'):
                neutral_rho_decay(inters)

            if not already_there(inters, 'Eta'):
                eta_decay(inters)

            if not already_there(inters, 'Charged pion'):
                charged_pion_decay(inters)

            if not already_there(inters, 'Neutral pion'):
                neutral_pion_decay(inters)

            inters += cls.four_particle_meson_decay(parts=[eta_prime, eta, charged_pion, charged_pion],
                                                    antiparts=[False, False, False, True],
                                                    ME=matrix_elements['EtaPrimePiPlusPiMinEta'],
                                                    kind=kind)
            inters += cls.four_particle_meson_decay(parts=[eta_prime, neutral_pion, neutral_pion, eta],
                                                    antiparts=[False, False, False, False],
                                                    ME=matrix_elements['EtaPrimeTwoPi0Eta'],
                                                    kind=kind)
            inters += cls.three_particle_meson_decay(parts=[eta_prime, neutral_rho, photon],
                                                    antiparts=[False, False, False],
                                                    ME=matrix_elements['EtaPrimeRho0Gamma'],
                                                    kind=kind)

        def phi_decay(inters=None):
            if not already_there(inters, 'Neutral rho'):
                neutral_rho_decay(inters)

            if not already_there(inters, 'Charged kaon'):
                charged_kaon_decay(inters)

            if not already_there(inters, 'Kaon long'):
                kaon_long_decay(inters)

            if not already_there(inters, 'Kaon short'):
                kaon_short_decay(inters)

            if not already_there(inters, 'Neutral pion'):
                neutral_pion_decay(inters)

            inters += cls.three_particle_meson_decay(parts=[phi, charged_kaon, charged_kaon],
                                                    antiparts=[False, True, False],
                                                    ME=matrix_elements['PhiKPlusKMin'],
                                                    kind=kind)
            inters += cls.three_particle_meson_decay(parts=[phi, kaon_long, kaon_short],
                                                    antiparts=[False, False, False],
                                                    ME=matrix_elements['PhiK0LK0S'],
                                                    kind=kind)
            inters += cls.three_particle_meson_decay(parts=[phi, neutral_rho, neutral_pion],
                                                    antiparts=[False, False, False],
                                                    ME=matrix_elements['PhiRho0Pi0'],
                                                    kind=kind)

        for meson in primary_mesons:
            if meson.name == 'Neutral pion':
                neutral_pion_decay(inters)

            if meson.name == 'Charged pion':
                charged_pion_decay(inters)

            if meson.name == 'Charged kaon':
                charged_kaon_decay(inters)

            if meson.name == 'Kaon long':
                kaon_long_decay(inters)

            if meson.name == 'Kaon short':
                kaon_short_decay(inters)

            if meson.name == 'Eta':
                eta_decay(inters)

            if meson.name == 'Neutral rho':
                neutral_rho_decay(inters)

            if meson.name == 'Charged rho':
                charged_rho_decay(inters)

            if meson.name == 'Omega':
                omega_decay(inters)

            if meson.name == 'Eta prime':
                eta_prime_decay(inters)

            if meson.name == 'Phi':
                phi_decay(inters)

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
