# -*- coding: utf-8 -*-
import numpy
from common import PicklableObject


class DistributionFunctional(object):
    """ == $\mathcal{F}(f_\alpha)$ functional == """

    """ === Naive form ===

        \begin{align}
            \mathcal{F} &= (1 \pm f_1)(1 \pm f_2) f_3 f_4 - f_1 f_2 (1 \pm f_3) (1 \pm f_4)
            \\\\ &= \mathcal{F}_B + \mathcal{F}_A
        \end{align}
    """

    def F_A(self, p, skip_index=None):
        """
        Forward reaction distribution functional term

        \begin{equation}
            \mathcal{F}_A = - f_1 f_2 (1 \pm f_3) (1 \pm f_4)
        \end{equation}

        :param skip_index: Particle to skip in the expression
        """
        temp = -1.

        for i, particle in enumerate(self.particles):
            if skip_index is None or i != skip_index:
                if self.signs[i] == 1:
                    temp *= particle.distribution(p[i])
                else:
                    temp *= 1. - particle.eta * particle.distribution(p[i])

        return temp

    def F_B(self, p, skip_index=None):
        """
        Backward reaction distribution functional term

        \begin{equation}
            \mathcal{F}_B = f_3 f_4 (1 \pm f_1) (1 \pm f_2)
        \end{equation}

        :param skip_index: Particle to skip in the expression
        """
        temp = 1.

        for i, particle in enumerate(self.particles):
            if skip_index is None or i != skip_index:
                if self.signs[i] == -1:
                    temp *= particle.distribution(p[i])
                else:
                    temp *= 1. - particle.eta * particle.distribution(p[i])

        return temp

    """
    === Linearized in $\, f_1$ form ===

    \begin{equation}
        \mathcal{F}(f) = f_3 f_4 (1 \pm f_1) (1 \pm f_2) - f_1 f_2 (1 \pm f_3) (1 \pm f_4)
    \end{equation}

    \begin{equation}
        \mathcal{F}(f) = f_1 (\mp f_3 f_4 (1 \pm f_2) - f_2 (1 \pm f_3) (1 \pm f_4)) \
        + f_3 f_4 (1 \pm f_2)
    \end{equation}

    \begin{equation}
        \mathcal{F}(f) = \mathcal{F}_B^{(1)} + f_1 (\mathcal{F}_A^{(1)} \pm_1 \mathcal{F}_B^{(1)})
    \end{equation}

    $^{(i)}$ in $\mathcal{F}^{(i)}$ means that the distribution function $f_i$ was omitted in the\
    corresponding expression. $\pm_j$ represents the $\eta$ value of the particle $j$.
    """
    def F_f(self, p):
        """ Variable part of the distribution functional """
        return (
            self.F_A(p=p, skip_index=0) - self.in_particles[0].eta * self.F_B(p=p, skip_index=0)
        )

    def F_1(self, p):
        """ Constant part of the distribution functional """
        return self.F_B(p=p, skip_index=0)


class BoltzmannIntegral(PicklableObject, DistributionFunctional):

    """ == Integral ==
        Representation of the concrete collision integral for a specific particle \
        `Integral.particles[0]` """

    _saveable_fields = [
        'particles', 'in_particles', 'out_particles',
        'signs', 'decoupling_temperature', 'constant', 'Ms',
    ]

    in_particles = []  # Incoming particles
    out_particles = []  # Outgoing particles
    particles = []  # All particles involved

    """ === Energy conservation law of the integral ===

        \begin{equation}
            0 = \vec{s} \cdot \vec{E} = \sum_i s_i E_i \sim E_0 + E_1 - E_2 - E_3
        \end{equation}
    """
    signs = [1, 1, -1, -1]

    # Temperature when the typical interaction time exceeds the Hubble expansion time
    decoupling_temperature = 0.
    constant = 0.

    """ Four-particle interactions of the interest can all be rewritten in a form

        \begin{equation}
            |\mathcal{M}|^2 = \sum_{\{i \neq j \neq k \neq l\}} K_1 (p_i \cdot p_j) (p_k \cdot p_l)\
                 + K_2 m_i m_j (p_k \cdot p_l)
        \end{equation} """

    Ms = []

    def __init__(self, **kwargs):
        """ Update self with configuration `kwargs`, construct particles list and \
            energy conservation law of the integral. """
        for key in kwargs:
            setattr(self, key, kwargs[key])

        self.particles = self.in_particles + self.out_particles
        self.signs = [1] * len(self.in_particles) + [-1] * len(self.out_particles)

        self.integrand = numpy.vectorize(self.integrand, otypes=[numpy.float_])

    def __str__(self):
        """ String-like representation of the integral. Corresponds to the first particle """
        return " + ".join([p.symbol for p in self.in_particles]) \
            + " ⟶  " + " + ".join([p.symbol for p in self.out_particles]) \
            + "\t({})".format(len(self.Ms))

    def __repr__(self):
        return self.__str__()

    def initialize(self):
        """
        Initialize collision integral constants and save them to the first involved particle
        """
        raise NotImplementedError()

    def calculate_kinematics(self, p):
        """ Helper procedure that caches conformal energies and masses of the particles """
        particle_count = len(self.particles)
        p = (p + [0., 0., 0., 0.])[:particle_count]
        E = []
        m = []
        for i, particle in enumerate(self.particles):
            E.append(particle.conformal_energy(p[i]))
            m.append(particle.conformal_mass)

        """ Parameters of one particle can be inferred from the energy conservation law
            \begin{equation}E_3 = -s_3 \sum_{i \neq 3} s_i E_i \end{equation} """
        E[particle_count-1] = -self.signs[particle_count-1] \
            * sum([self.signs[i] * E[i] for i in range(particle_count-1)])
        p[particle_count-1] = numpy.sqrt(numpy.abs(E[particle_count-1]**2 - m[particle_count-1]**2))
        return p, E, m

    def correction(self, p0):
        """ === Particle collisions integration === """

        particle = self.particles[0]

        integral_1, _ = self.integral_1(p0)
        integral_f, _ = self.integral_f(p0)

        order = min(len(particle.data['collision_integral']) + 1, 5)

        index = numpy.argwhere(GRID.TEMPLATE == p0)[0][0]
        fs = [i[index] for i in particle.data['collision_integral'][-order:]]

        prediction = integrators.adams_moulton_solver(
            y=particle.distribution(p0), fs=fs,
            A=integral_1, B=integral_f,
            h=particle.params.dy, order=order
        )

        total_integral = (prediction - particle.distribution(p0)) / particle.params.dy

        return total_integral

    def integrate(self, p0, integrand, bounds=None, kwargs=None):
        raise NotImplementedError()

    def integral_1(self, p0, bounds=None):
        return self.integrate(
            p0, self.integrand,
            kwargs={'fau': self.F_1},
            name=str(self), bounds=bounds
        )

    def integral_f(self, p0, bounds=None):
        return self.integrate(
            p0, self.integrand,
            kwargs={'fau': self.F_f},
            name=str(self), bounds=bounds
        )

    def integrand(self, *args, **kwargs):
        """ Collision integral interior. """
        raise NotImplementedError()

    """ === Integration region bounds methods === """

    def in_bounds(self, p, E=None, m=None):
        raise NotImplementedError()

    def bounds(self, p0):
        """ Coarse integration region based on the `GRID` points. Assumes that integration region\
            is connected. """
        points = []
        for p1 in GRID.TEMPLATE:
            points.append((p1, self.lower_bound(p0, p1),))
            points.append((p1, self.upper_bound(p0, p1),))

        return points

    def lower_bound(self, p0, p1):
        """ Find the first `GRID` point in the integration region """

        index = 0
        while index < GRID.MOMENTUM_SAMPLES and not self.in_bounds([p0, p1, GRID.TEMPLATE[index]]):
            index += 1

        if index == GRID.MOMENTUM_SAMPLES:
            return GRID.MIN_MOMENTUM

        return GRID.TEMPLATE[index]

    def upper_bound(self, p0, p1):
        """ Find the last `GRID` point in the integration region """

        index = int((min(p0 + p1, GRID.MAX_MOMENTUM) - GRID.MIN_MOMENTUM) / GRID.MOMENTUM_STEP)

        while index >= 0 and not self.in_bounds([p0, p1, GRID.TEMPLATE[index]]):
            index -= 1

        if index == -1:
            return GRID.MIN_MOMENTUM

        return GRID.TEMPLATE[index]
