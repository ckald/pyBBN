import os


defaults = {
    # Whether the code should use a fixed order Gauss-Legendre quadrature
    # or adaptive NumPy method `quad` (based on QUADPACK library)
    'FIXED_ORDER_1D_QUADRATURE': True,

    # The order of the Gauss-Legendre quadrature used by the code for momentum-space integrations
    'GAUSS_LEGENDRE_ORDER': 30,
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
    'ADAMS_MOULTON_DISTRIBUTION_CORRECTION': True,

    # The default number of points on the momentum space grid
    'MOMENTUM_SAMPLES': 51,
    # The maximal value on the momentum space grid in MeV
    'MAX_MOMENTUM_MEV': 20,

    # The fraction $\gamma = m/T$ at which the equilibrium particle switch
    'REGIME_SWITCHING_FACTOR': 1e3,

    # Maximum momentum used in integration of massive equilibrium particle
    # contribution to temperature
    'TEMPERATURE_INTEGRAL_MAX_MOMENTUM_MEV': 200,

    # Use Simpson's rule for computation of collision integrals instead of Gaussian quadrature
    'SIMPSONS_INTEGRATION': False,

    'LAGUERRE_GAUSS_FOR_MASSIVE_EQUILIBRIUM_PARTICLES': False,

}


def get(name):
    default = defaults.get(name)
    val = os.environ.get(name, default)

    if default:
        typ = type(default)
        val = typ(val)

    return val
