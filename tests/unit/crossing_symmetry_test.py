from collections import defaultdict, Counter
from . import non_equilibium_setup, with_setup_args, setup
from common import CONST, UNITS
from evolution import Universe
from particles import Particle
from library.SM import particles as SMP
from library.NuMSM import particles as NuP, interactions as NuI


@with_setup_args(non_equilibium_setup)
def neutrino_scattering_amplitude_test(params, universe):

    params.update(universe.total_energy_density(), universe.total_entropy())

    photon, neutrino_e, neutrino_mu = universe.particles

    assert len(universe.interactions) == 1
    assert len(universe.interactions[0].integrals) == 1
    integral = universe.interactions[0].integrals[0]
    assert len(integral.Ms) == 2

    assert any(M.order == (0, 3, 1, 2) and M.K2 == 0 and M.K1 == 4 * 32 * CONST.G_F**2
               for M in integral.Ms)
    assert any(M.order == (0, 1, 2, 3) and M.K2 == 0 and M.K1 == 2 * 32 * CONST.G_F**2
               for M in integral.Ms)


@with_setup_args(setup)
def three_particle_integral_heavy_test(params):
    """ If M_N > M_pi, there should be integrals for the reactions:

        N <--> nu_e + pi^0
        nu_e + pi^0 <--> N
        pi^0 +nu_e <--> N
        pi^0 + anti-nu_e <--> anti-N
    """

    photon = Particle(**SMP.photon)
    neutrino_e = Particle(**SMP.leptons.neutrino_e)
    sterile = Particle(**NuP.dirac_sterile_neutrino(mass=200 * UNITS.MeV))
    neutral_pion = Particle(**SMP.hadrons.neutral_pion)

    theta = 1e-3
    thetas = defaultdict(float, {
        'electron': theta,
    })

    interaction = NuI.sterile_hadrons_interactions(
        thetas=thetas, sterile=sterile,
        neutrinos=[neutrino_e],
        leptons=[],
        mesons=[neutral_pion]
    )

    universe = Universe(params=params)
    universe.add_particles([photon, neutrino_e, sterile, neutral_pion])
    universe.interactions += interaction

    assert len(universe.interactions) == 2
    assert len(universe.interactions[0].integrals) == 2
    assert len(universe.interactions[1].integrals) == 2

    integral = universe.interactions[0].integrals[0]
    assert len(integral.Ms) == 1
    assert isinstance(integral.Ms[0].K, (int, float))


@with_setup_args(setup)
def three_particle_integral_light_test(params):
    """ If M_N < M_pi, there should be integrals for the reactions:

        N + anti-nu_e <--> pion
        nu_e + anti-N <--> pion
        pi^0 <--> anti-nu_e + N
        pi^0 <--> nu_e + anti-N
    """

    photon = Particle(**SMP.photon)
    neutrino_e = Particle(**SMP.leptons.neutrino_e)
    sterile = Particle(**NuP.dirac_sterile_neutrino(mass=100 * UNITS.MeV))
    neutral_pion = Particle(**SMP.hadrons.neutral_pion)

    theta = 1e-3
    thetas = defaultdict(float, {
        'electron': theta,
    })

    interaction = NuI.sterile_hadrons_interactions(
        thetas=thetas, sterile=sterile,
        neutrinos=[neutrino_e],
        leptons=[],
        mesons=[neutral_pion]
    )

    universe = Universe(params=params)
    universe.add_particles([photon, neutrino_e, sterile, neutral_pion])
    universe.interactions += interaction

    assert len(universe.interactions) == 2
    assert len(universe.interactions[0].integrals) == 2
    assert len(universe.interactions[1].integrals) == 2

    integral = universe.interactions[0].integrals[0]
    assert len(integral.Ms) == 1
    assert isinstance(integral.Ms[0].K, (int, float))


@with_setup_args(setup)
def four_particle_integrals_check_test(params):
    """ Checks whether correct amount of collision integrals is created by the
        CrossingGenerating function.

        Expect: 10 integrals for NC (3x HNL, 3x neutrino, 4x electron)
                16 integrals for CC (4x HNL, 4x neutrino, 4x muon, 4x electron)
    """
    neutrino_e = Particle(**SMP.leptons.neutrino_e)
    neutrino_mu = Particle(**SMP.leptons.neutrino_mu)
    neutrino_tau = Particle(**SMP.leptons.neutrino_tau)
    electron = Particle(**SMP.leptons.electron)
    muon = Particle(**SMP.leptons.muon)
    sterile = Particle(**NuP.dirac_sterile_neutrino(mass=150 * UNITS.MeV))

    theta = 1e-3
    thetas = defaultdict(float, {
        'electron': theta,
    })

    interaction = NuI.sterile_leptons_interactions(
        thetas=thetas, sterile=sterile,
        neutrinos=[neutrino_e, neutrino_mu, neutrino_tau],
        leptons=[electron, muon]
    )

    NC_HNL = 0
    NC_nu = 0
    NC_e = 0
    CC_HNL = 0
    CC_nu = 0
    CC_mu = 0
    CC_e = 0

    for inter in interaction:
        for integral in inter.integrals:
            count = Counter(item.specie for item in integral.reaction)
            if count[neutrino_e] and count[electron] and count[sterile]:
                if integral.particle == sterile:
                    NC_HNL += 1
                if integral.particle == neutrino_e:
                    NC_nu += 1
                if integral.particle == electron:
                    NC_e += 1

            if count[neutrino_mu] and count[electron] and count[muon] and count[sterile]:
                if integral.particle == sterile:
                    CC_HNL += 1
                if integral.particle == neutrino_mu:
                    CC_nu += 1
                if integral.particle == muon:
                    CC_mu += 1
                if integral.particle == electron:
                    CC_e += 1

    tpl = (
        "Incorrect amount of collision integrals for {part} in {int}" +
        " interactions created. Expected: {exp}. Obtained: {val}"
    )

    assert NC_HNL == 3, tpl.format(part="HNL", int="NC", exp=3, val=NC_HNL)
    assert NC_nu == 3, tpl.format(part="nu_e", int="NC", exp=3, val=NC_nu)
    assert NC_e == 4, tpl.format(part="e", int="NC", exp=4, val=NC_e)

    assert CC_HNL == 4, tpl.format(part="HNL", int="CC", exp=4, val=CC_HNL)
    assert CC_nu == 4, tpl.format(part="nu", int="CC", exp=4, val=CC_nu)
    assert CC_mu == 4, tpl.format(part="mu", int="CC", exp=4, val=CC_mu)
    assert CC_e == 4, tpl.format(part="e", int="CC", exp=4, val=CC_e)


@with_setup_args(setup)
def neutral_channel_test(params):
    neutrino = Particle(**SMP.leptons.neutrino_e)
    lepton = Particle(**SMP.leptons.electron)
    sterile = Particle(**NuP.dirac_sterile_neutrino(mass=150 * UNITS.MeV))

    theta = 1e-3
    thetas = defaultdict(float, {
        'electron': theta,
    })

    interaction = NuI.sterile_active_to_leptons_NC(
        theta=thetas[neutrino.flavour],
        g_L=CONST.g_L,
        sterile=sterile,
        active=neutrino,
        lepton=lepton
    )

    counter = Counter(integral.particle.name
                      for integral in interaction.integrals)
    assert counter[sterile.name] == 3
    assert counter[neutrino.name] == 3
    assert counter[lepton.name] == 4

    assert len(interaction.integrals) == 3
