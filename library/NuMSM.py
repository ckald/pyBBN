# -*- coding: utf-8 -*-
"""
# Neutrino Minimal Standard Model particles and interactions
"""

from __future__ import division

import itertools

from common import UNITS, CONST, statistics as STATISTICS
from interactions import CrossGeneratingInteraction
from interactions.three_particle import ThreeParticleM, ThreeParticleIntegral
from interactions.four_particle import FourParticleIntegral
from library.SM import WeakM, particles as SM_particles


class SterileM(WeakM):

    """ ## Sterile interactions matrix element
        Sterile processes usually include a common factor of $32 θ^2 G_F^2$ """

    def __init__(self, theta=1., *args, **kwargs):
        super(SterileM, self).__init__(*args, **kwargs)

        self.theta = theta
        self.K1 *= self.theta**2
        self.K2 *= self.theta**2

    def __str__(self):
        """ String-like representation of the matrix element """
        ret = ""
        if self.K1:
            ret += "K1=32 θ^2 G_F^2 {: .3e} ".format(self.K1 / self.const / self.theta**2)
        if self.K2:
            ret += "K2=32 θ^2 G_F^2 {: .3e} ".format(self.K2 / self.const / self.theta**2)
        return ret + self.order_format() + " θ ={: .3e}".format(self.theta)


class particles(object):
    @staticmethod
    def sterile_neutrino(mass=33.9 * UNITS.MeV):
        return {
            'name': 'Sterile neutrino (Majorana)',
            'symbol': 'N',
            'statistics': STATISTICS.FERMION,
            'mass': mass,
            'dof': 2,
            'decoupling_temperature': 1. * UNITS.GeV,
            'majorana': True
        }

    @staticmethod
    def dirac_sterile_neutrino(mass=33.9 * UNITS.MeV):
        return {
            'name': 'Sterile neutrino (Dirac)',
            'symbol': 'N',
            'statistics': STATISTICS.FERMION,
            'mass': mass,
            'dof': 4,
            'decoupling_temperature': 1. * UNITS.GeV,
            'majorana': False
        }


class interactions(object):

    @staticmethod
    def sterile_active_scattering(theta=1., sterile=None, active_a=None,
                                 active_b=None, kind=None):
        """ \begin{align}
                N_S + \nu_{\beta} &\to \nu_{\alpha} + \nu_{\beta}
            \end{align}
        """

        # SOMETHING IS CLEARLY NOT RIGHT
        K1 = 2. if active_a == active_b else 1

        return CrossGeneratingInteraction(
            name="Sterile-active neutrino scattering",
            particles=((sterile, active_b), (active_a, active_b)),
            antiparticles=((False, False), (False, False)),
            Ms=(SterileM(theta=theta, K1=K1, order=(0, 1, 2, 3)), ),
            integral_type=FourParticleIntegral,
            kind=kind
        )

    @staticmethod
    def sterile_active_to_leptons_NC(theta=1., g_L=CONST.g_R+0.5, sterile=None,
                                     active=None, lepton=None, kind=None):
        """ \begin{align}
                N_S + \overline{\nu_{\alpha}} \to l^+ + l^-
            \end{align}
        """

        return CrossGeneratingInteraction(
            name="Sterile decay into neutrino and leptons pair",
            particles=((sterile, active), (lepton, lepton)),
            antiparticles=((False, True), (True, False)),
            Ms=(
                SterileM(theta=theta, K1=4. * g_L**2, order=(0, 2, 1, 3)),
                SterileM(theta=theta, K1=4. * CONST.g_R**2, order=(0, 3, 1, 2)),
                SterileM(theta=theta, K2=4. * CONST.g_R * g_L, order=(2, 3, 0, 1))
            ),
            integral_type=FourParticleIntegral,
            kind=kind
        )

    @staticmethod
    def sterile_active_to_leptons_CC(theta=1., sterile=None,
                                     active=None, lepton1=None, lepton2=None, kind=None):
        """ \begin{align}
                N_S + \overline{\nu_{\alpha}} \to l1^+ + l2^-
            \end{align}
        """

        return CrossGeneratingInteraction(
            name="Sterile decay into neutrino and leptons",
            particles=((sterile, active), (lepton2, lepton1)),
            antiparticles=((False, True), (True, False)),
            Ms=(
                SterileM(theta=theta, K1=4., order=(0, 2, 1, 3)),
            ),
            integral_type=FourParticleIntegral,
            kind=kind
        )

    @classmethod
    def sterile_leptons_interactions(cls, thetas=None, sterile=None,
                                     neutrinos=None, leptons=None, kind=None):

        g_R = CONST.g_R
        inters = []

        # Neutrinos scatterings
        for neutrino_a, neutrino_b in itertools.product(neutrinos, neutrinos):
            if thetas[neutrino_a.flavour]:
                inters.append(cls.sterile_active_scattering(
                    theta=thetas[neutrino_a.flavour],
                    sterile=sterile,
                    active_a=neutrino_a,
                    active_b=neutrino_b,
                    kind=kind
                ))

        # Interactions of neutrinos and leptons
        for lepton in leptons:
            for neutrino in neutrinos:
                if thetas[neutrino.flavour]:
                    g_L = g_R + 0.5 if lepton.flavour == neutrino.flavour else g_R - 0.5
                    # Neutral current
                    inters.append(cls.sterile_active_to_leptons_NC(
                        theta=thetas[neutrino.flavour],
                        g_L=g_L,
                        sterile=sterile,
                        active=neutrino,
                        lepton=lepton,
                        kind=kind
                    ))

                if len(leptons) > 1:
                    other_lepton = [par for par in leptons if par != lepton][0]
                    if thetas[lepton.flavour] and neutrino.flavour == other_lepton.flavour:
                        # Charged current
                        inters.append(cls.sterile_active_to_leptons_CC(
                            theta=thetas[lepton.flavour],
                            sterile=sterile,
                            active=neutrino,
                            lepton1=lepton,
                            lepton2=other_lepton,
                            kind=kind
                        ))
        return inters

    # ## Quark interactions

    @staticmethod
    def sterile_quark_interactions(thetas=None, sterile=None,
                                   neutrinos=None, leptons=None, quarks=None):

        inters = []

        for quark in quarks:
            for active in neutrinos:
                if thetas[active.flavour]:
                    inters.append(interactions.sterile_quark_neutral(
                        theta=thetas[active.flavour], sterile=sterile,
                        active=active, quark=quark
                    ))

        for up, down in itertools.product([q for q in quarks if q.Q == 2./3.],
                                          [q for q in quarks if q.Q == -1./3.]):
            for lepton in leptons:
                if thetas[lepton.flavour]:
                    inters.append(interactions.sterile_quark_charged(
                        theta=thetas[lepton.flavour], sterile=sterile,
                        lepton=lepton, up=up, down=down
                    ))

        return inters

    @staticmethod
    def sterile_quark_neutral(theta=1., sterile=None, active=None, quark=None):
        """ \begin{align}
                N_S + \overline{\nu} &\to q + \overline{q}
            \end{align}
        """

        g_R = CONST.g_R

        if quark.Q == 2./3.:
            x = (3 - 4 * g_R)
        if quark.Q == -1./3.:
            x = (3 - 2 * g_R)

        return CrossGeneratingInteraction(
            name="Sterile-quarks neutral channel interaction",
            particles=((sterile, active), (quark, quark)),
            antiparticles=((False, True), (False, True)),
            Ms=(
                SterileM(theta=theta, K1=32. / 9. * g_R, order=(0, 3, 1, 2)),
                SterileM(theta=theta, K1=2. / 9. * x**2, order=(0, 2, 1, 3)),
                SterileM(theta=theta, K2=16. / 9. * g_R * x, order=(0, 1, 2, 3)),
            ),
            integral_type=FourParticleIntegral
        )

    @staticmethod
    def sterile_quark_charged(theta=1., sterile=None, lepton=None, up=None, down=None):
        """ \begin{align}
                N_S + l^+ &\to u + \overline{d}
            \end{align}
        """

        CKM = SM_particles.quarks.CKM[(up.family, down.family)]

        return CrossGeneratingInteraction(
            name="Sterile-quarks charged channel interaction",
            particles=((sterile, lepton), (up, down)),
            antiparticles=((False, True), (False, True)),
            Ms=(
                SterileM(theta=theta, K1=.5 * CKM**2, order=(0, 2, 1, 3)),
                SterileM(theta=theta, K1=.5 * CKM**2, order=(0, 3, 1, 2)),
            ),
            integral_type=FourParticleIntegral
        )

    @classmethod
    def sterile_hadrons_interactions(cls, thetas=None, sterile=None, neutrinos=None,
                                     leptons=None, mesons=None, kind=None):

        inters = []

        for neutrino in neutrinos:
            if thetas[neutrino.flavour]:
                for particle in mesons:
                    if particle.Q == 0:
                        if particle.type == "scalar":
                            inters += cls.neutral_scalar_meson(theta=thetas[neutrino.flavour],
                                                               sterile=sterile,
                                                               active=neutrino,
                                                               meson=particle,
                                                               kind=kind)

                        if particle.type == "vector":
                            inters += cls.neutral_vector_meson(theta=thetas[neutrino.flavour],
                                                               sterile=sterile,
                                                               active=neutrino,
                                                               meson=particle,
                                                               kind=kind)

        for lepton in leptons:
            if thetas[lepton.flavour]:
                for particle in mesons:
                    if particle.Q == -1:
                        if particle.type == "scalar":
                            inters += cls.charged_scalar_meson(theta=thetas[lepton.flavour],
                                                               sterile=sterile,
                                                               lepton=lepton,
                                                               meson=particle,
                                                               kind=kind)

                        if particle.type == "vector":
                            inters += cls.charged_vector_meson(theta=thetas[neutrino.flavour],
                                                               sterile=sterile,
                                                               lepton=lepton,
                                                               meson=particle,
                                                               kind=kind)

        return inters

    @staticmethod
    def neutral_scalar_meson(theta=1., sterile=None, active=None, meson=None, kind=None):
        """ \begin{align}
                N \to \nu + meson
            \end{align}
        """
        return [CrossGeneratingInteraction(
            name="Sterile neutrino decay to neutral scalar meson and neutrino",
            particles=((sterile, ), (active, meson)),
            antiparticles=antiparticles,
            Ms=(
                ThreeParticleM(K=(CONST.G_F * theta * meson.decay_constant)**2
                                  * sterile.mass**4 * (1 - (meson.mass / sterile.mass)**2)), ),
            integral_type=ThreeParticleIntegral,
            kind=kind
        ) for antiparticles in [
            ((False, ), (False, False)),
            ((True, ), (True, False))
        ]]

    @staticmethod
    def charged_scalar_meson(theta=1., sterile=None, lepton=None, meson=None, kind=None):
        """ \begin{align}
                N \to l + meson
            \end{align}
        """

        CKM = SM_particles.quarks.CKM[(1, 1)]

        return [CrossGeneratingInteraction(
            name="Sterile neutrino decay to charged scalar meson and lepton",
            particles=((sterile, ), (lepton, meson)),
            antiparticles=antiparticles,
            Ms=(ThreeParticleM(
                K=(CONST.G_F * theta * meson.decay_constant * CKM)**2 * sterile.mass**4 * (
                    (1 - (lepton.mass / sterile.mass)**2)**2
                    - (meson.mass / sterile.mass)**2 * (1 + (lepton.mass / sterile.mass)**2)
                )
            ),),
            integral_type=ThreeParticleIntegral,
            kind=kind
        ) for antiparticles in [
            ((False, ), (False, True)),
            ((True, ), (True, False))
        ]]

    @staticmethod
    def neutral_vector_meson(theta=1., sterile=None, active=None, meson=None, kind=None):
        """ \begin{align}
                N \to \nu + meson
            \end{align}
        """

        if meson.name in ["Neutral rho"]:
            return [CrossGeneratingInteraction(
                name="Sterile neutrino decay to neutral vector meson and neutrino",
                particles=((sterile, ), (active, meson)),
                antiparticles=antiparticles,
                Ms=(ThreeParticleM(K=0.5 * (CONST.G_F * theta * meson.decay_constant)**2
                                   * sterile.mass**4 * (1 - 2 * CONST.sin_theta_w_2)**2
                                   * (1 + 2 * (meson.mass / sterile.mass)**2)
                                   * (1 - (meson.mass / sterile.mass)**2)), ),
                integral_type=ThreeParticleIntegral,
                kind=kind
            ) for antiparticles in [
                ((False, ), (False, False)),
                ((True, ), (True, False))
            ]]

        if meson.name in ["Omega", "Phi"]:
            return [CrossGeneratingInteraction(
                name="Sterile neutrino decay to neutral vector meson and neutrino",
                particles=((sterile, ), (active, meson)),
                antiparticles=antiparticles,
                Ms=(ThreeParticleM(K=8 * (CONST.G_F * theta * meson.decay_constant)**2
                                   * sterile.mass**4 * (CONST.sin_theta_w_2 / 3)**2
                                   * (1 + 2 * (meson.mass / sterile.mass)**2)
                                   * (1 - (meson.mass / sterile.mass)**2)), ),
                integral_type=ThreeParticleIntegral,
                kind=kind
            ) for antiparticles in [
                ((False, ), (False, False)),
                ((True, ), (True, False))
            ]]


    @staticmethod
    def charged_vector_meson(theta=1., sterile=None, lepton=None, meson=None, kind=None):
        """ \begin{align}
                N \to l + meson
            \end{align}
        """

        CKM = SM_particles.quarks.CKM[(1, 1)]

        return [CrossGeneratingInteraction(
            name="Sterile neutrino decay to charged vector meson and lepton",
            particles=((sterile, ), (lepton, meson)),
            antiparticles=antiparticles,
            Ms=(ThreeParticleM(
                K=(CONST.G_F * theta * meson.decay_constant * CKM)**2 * sterile.mass**4 * (
                    1 + (lepton.mass / sterile.mass)**4 + (lepton.mass / sterile.mass)**2
                    * ((meson.mass / sterile.mass)**2 - 2) + (meson.mass / sterile.mass)**2
                    - 2 * (meson.mass / sterile.mass)**4
                )
            ),),
            integral_type=ThreeParticleIntegral,
            kind=kind
        ) for antiparticles in [
            ((False, ), (False, True)),
            ((True, ), (True, False))
        ]]
