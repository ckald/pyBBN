import numpy


def Fermi(e, mu=0):
    """ Fermi-Dirac:
        \begin{equation}
            \frac{1}{e^E + 1}
        \end{equation}
    """
    return 1. / (numpy.exp(e - mu) + 1.)


def Bose(e, mu=0):
    """ Bose-Einstein:
        \begin{equation}
            \frac{1}{e^E - 1}
        \end{equation}
    """
    return 1. / (numpy.exp(e - mu) - 1.)

BOSON = Bose
FERMION = Fermi
