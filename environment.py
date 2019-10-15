import os


defaults = {
    'SPLIT_COLLISION_INTEGRAL': True,

    # Whether the code should use a fixed order Gauss-Legendre quadrature
    # or adaptive NumPy method `quad` (based on QUADPACK library)
    # NOTE: QUADPACK integrator gives incorrect results when using high-resolution grids
    'FIXED_ORDER_1D_QUADRATURE': True,

    # The order of the Gauss-Legendre quadrature used by the code for momentum-space integrations
    'GAUSS_LEGENDRE_ORDER': 100,
    # The order of the Gauss-Laguerre quadrature used by the code for momentum-space integrations
    'GAUSS_LAGUERRE_ORDER': 100,

    # Whether the code should formulate the temperature equation in terms of
    # a logarithm of scale factor or dimensionful scale factor ($x = a * 1 MeV$)
    'LOGARITHMIC_TIMESTEP': True,

    # Whether the code should use Adams-Bashforth or explicit Euler numerical scheme
    # while solving for the temperature evolution
    'ADAMS_BASHFORTH_TEMPERATURE_CORRECTION': True,

    # Whether the code should use Adams-Moulton or implicit Euler numerical scheme
    # while solving for the distribution function evolution
    'ADAMS_MOULTON_DISTRIBUTION_CORRECTION': False,

    # The default number of points on the momentum space grid
    'MOMENTUM_SAMPLES': 401,
    # The maximal value on the momentum space grid in MeV
    'MAX_MOMENTUM_MEV': 100,

    # The maximal value of the scale factor till which collision integrals will be computed
    # based on a photon-electron only universe
    'MAX_SCALE_FACTOR': 10,

    # The resolution on the four particle momentum space grid in MeV / step
    'FOUR_PARTICLE_GRID_RESOLUTION': 0.24999,

    # The hierarchy of neutrinos (normal or inverted)
    'NORMAL_HIERARCHY_NEUTRINOS' : True,

    # The fraction $\gamma = m/T$ at which the equilibrium particle switch
    'REGIME_SWITCHING_FACTOR': 1e3,

    # Use Simpson's rule for computation of collision integrals instead of Gaussian quadrature
    'SIMPSONS_INTEGRATION': False,

    # Use Simpson's rule for computation of thermodynamical quantities of non-equilibrium particles
    # instead of Gaussian quadrature
    # NOTE: Gaussian quadrature (wrongly) gives increasing comoving number density and entropy for decoupled species.
    'SIMPSONS_NONEQ_PARTICLES': True,

    'LAGUERRE_GAUSS_FOR_MASSIVE_EQUILIBRIUM_PARTICLES': True,

}


def get(name):
    default = defaults.get(name)
    val = os.environ.get(name, default)

    if default:
        typ = type(default)
        val = typ(val)

    return val
