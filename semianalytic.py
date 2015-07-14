from common import UNITS


T_BBN = 0.065 * UNITS.MeV


def step_monitor(universe):
    from copy import deepcopy

    if not hasattr(universe, 'neutron_decoupling_parameters') and len(universe.kawano_data):
        datarow = universe.kawano_data.tail(1)
        rates = datarow.values[0][-6:]
        neutron_equilibration = (sum(rates[0::2]) / sum(rates[1::2]))
        if neutron_equilibration > 10 or neutron_equilibration < 0.1:
            universe.neutron_decoupling_parameters = deepcopy(universe.params)
            print "Neutron decoupled at T = {:e} MeV".format(universe.params.T / UNITS.MeV)

    if not hasattr(universe, 'deuterium_generation_parameters') and universe.params.T <= T_BBN:
            print "Deuterium generation"
            universe.deuterium_generation_parameters = deepcopy(universe.paramsx)

    if (hasattr(universe, 'neutron_decoupling_parameters')
            and hasattr(universe, 'deuterium_generation_parameters')):

        print '_' * 80
        print 'BBN estimates'

        nparams = universe.neutron_decoupling_parameters
        dparams = universe.deuterium_generation_parameters

        neutron_decay_time = nparams.t - dparams.t
        neutron_lifetime = 885.7 * UNITS.s

        import math
        neutron_to_proton = (
            math.exp(-universe.kawano.q / nparams.T)
            * math.exp(-neutron_decay_time / neutron_lifetime)
        )

        helium_fraction = 2 * neutron_to_proton / (1 + neutron_to_proton)

        print "Neutron decoupling: T = {:e} MeV, t = {:e} s, N_eff = {:e}".format(
            nparams.T / UNITS.MeV, nparams.t / UNITS.s, nparams.N_eff
        )
        print "Deuterium generation: T = {:e} MeV, t = {:e} s, N_eff = {:e}".format(
            dparams.T / UNITS.MeV, dparams.t / UNITS.s, dparams.N_eff
        )
        print "Neutron decay time: {:e} s".format(neutron_decay_time / UNITS.s)
        print "Neutron-to-proton ratio: {:e}".format(neutron_to_proton)
        print "Helium fraction: {:e}".format(helium_fraction)
        print '_' * 80
