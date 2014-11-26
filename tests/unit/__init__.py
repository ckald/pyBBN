from common import PARAMS, UNITS
from library import StandardModelParticles as SMP


eps = 1e-5


def setup():
    PARAMS.T_initial = SMP.neutrino_e['decoupling_temperature']
    PARAMS.T_final = 0.075 * UNITS.MeV
    PARAMS.dx = 1e-4 * UNITS.MeV
    PARAMS.infer()
