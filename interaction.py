# -*- coding: utf-8 -*-
import copy
import numpy
from common import PARAMS, GRID, CONST, UNITS, integrators
from common.utils import PicklableObject, benchmark
from ds import D, Db1, Db2


class INTERACTIONS:
    DECAY = 'decay'


"""
== Boltzmann collision integral ==

\begin{multline}
    I_{coll}(t,p_i) = \frac{1}{2 E_i} \sum_{reactions} \int \cdots \int
        \prod_{j} \frac{d^3 p_{in,j}}{(2 \pi)^3 2 E_{in,j}}
        \prod_{k} \frac{d^3 p_{out,k}}{(2 \pi)^3 2 E_{out,k}} \\
        S |\mathcal{M}|^2 \mathcal{F}(\{f_\alpha\}) (2 \pi)^4 \\
        \delta^4(\sum_{m} p_{in,m}^\mu - \sum_{n} p_{out,n}^\mu)
\end{multline}
"""


class M(object):

    """ == Matrix element ==
        All four-particle interactions of the interest can be rewritten in a form

        \begin{equation}
            |\mathcal{M}|^2 = \sum_{\{i \neq j \neq k \neq l\}} K_1 (p_i \cdot p_j) (p_k \cdot p_l)\
                 + K_2 m_i m_j (p_k \cdot p_l)
        \end{equation} """
    K1 = 0.
    K2 = 0.
    # Order defines the values of the $(i, j, k, l)$ indices
    order = (0, 1, 2, 3)

    def __init__(self, *args, **kwargs):
        """ Configure matrix element, check that it makes sense """
        for key in kwargs:
            setattr(self, key, kwargs[key])

        if len(set(self.order)) != len(self.order):
            raise Exception("Meaningless order of momenta of the matrix element: {}"
                            .format(self.order))

    def __str__(self):
        """ String-like representation of the matrix element """
        return "K1={: .2e}, K2={: .2e}, {}".format(self.K1, self.K2, self.order)


class WeakM(M):

    """ == Weak interactions matrix element ==
        Weak processes usually include a common factor of $32 G_F^2$ """

    def __init__(self, *args, **kwargs):
        super(WeakM, self).__init__(*args, **kwargs)

        self.K1 *= 32 * CONST.G_F**2
        self.K2 *= 32 * CONST.G_F**2


class Interaction(object):

    """
    == Interaction ==
    Helper class that takes care of creating all necessary `Integral`s for the actual interaction.
    """

    # Human-friendly interaction identifier
    name = "Particle interaction"

    integrals = []

    in_particles = []  # Incoming particles
    out_particles = []  # Outgoing particles
    particles = []  # All particles involved

    # Temperature when the typical interaction time exceeds the Hubble expansion time
    decoupling_temperature = 0.

    # Matrix elements of the interaction
    Ms = []

    def __init__(self, *args, **kwargs):
        """ Create an `Integral` object for all particle species involved in the interaction.

            Precise expressions for all integrals can be derived by permuting all particle-related\
            functions in the distribution functions, matrix elements, momenta, energy and mass\
            arrays.

            To avoid double-counting, one should create an integral for each particle specie only \
            once.
        """

        for key in kwargs:
            setattr(self, key, kwargs[key])

        self.integrals = []
        self.particles = self.in_particles + self.out_particles

        Ms = copy.deepcopy(self.Ms)

        # Remember all particle species we've already considered to avoid double-counting
        accounted_particles = set()

        # === Interaction integrals initialization strategy ===

        # ==== 1. Permute all in_particles ====
        in_particles = self.in_particles
        out_particles = self.out_particles
        accounted_particles = self.init_integrals(in_particles, out_particles,
                                                  Ms, accounted_particles)

        # ==== 2. Turn to the backward process ====
        # `(0, 1, 2, 3) -> (2, 3, 0, 1)`
        for M in Ms:
            M.order = M.order[len(self.in_particles):] + M.order[:len(self.in_particles)]

        # ==== 3. Permute all new `in_particles` (former `out_particles`) ====
        accounted_particles = self.init_integrals(out_particles, in_particles,
                                                  Ms, accounted_particles)

    def __str__(self):
        """ String-like representation of the interaction and all its integral """
        return self.name + "\n\t" + "\n\t".join([str(integral) for integral in self.integrals])

    def init_integrals(self, in_particles, out_particles, Ms, accounted_particles):
        for i, particle in enumerate(in_particles):

            # Skip already accounted species
            if particle in accounted_particles:
                continue

            particle_Ms = copy.deepcopy(Ms)

            for M in particle_Ms:
                M.order = M.order[i:i+1] + M.order[:i] + M.order[i+1:]

            # Add interaction integrals by putting each incoming particle as the first one
            self.integrals.append(Integral(
                in_particles=in_particles[i:i+1] + in_particles[:i] + in_particles[i+1:],
                out_particles=out_particles,
                decoupling_temperature=self.decoupling_temperature,
                Ms=particle_Ms
            ))
            accounted_particles.add(particle)
        return accounted_particles

    def initialize(self):
        """ Proxy method """

        for integral in self.integrals:
            integral.initialize()


class Integral:

    """ == Integral ==
        Representation of the concrete collision integral for a specific particle \
        `Integral.particles[0]` """

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

    DETAILED_OUTPUT = False
    COLLECTIVE_INTEGRATION = False

    def __init__(self, *args, **kwargs):
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
        if PARAMS.T > self.decoupling_temperature and not self.in_particles[0].in_equilibrium:
            self.constant = 1./64. / numpy.pi**3 * PARAMS.m**5 / PARAMS.x**5 / PARAMS.H
            self.particles[0].collision_integrals.append(self)

    def calculate_kinematics(self, p):
        """ Helper procedure that caches conformal energies and masses of the particles """
        p = (p + [0., 0., 0., 0.])[:4]
        E = []
        m = []
        for i, particle in enumerate(self.particles):
            E.append(particle.conformal_energy(p[i]))
            m.append(particle.conformal_mass)

        """ Parameters of one particle can be inferred from the energy conservation law
            \begin{equation}E_3 = -\frac{1}{s_3} \sum_{i \neq 3} s_i E_i \end{equation} """
        E[3] = -self.signs[3] * sum([self.signs[i] * E[i] for i in range(3)])
        p[3] = numpy.sqrt(numpy.abs(E[3]**2 - m[3]**2))
        return p, E, m

    def benchmarked_integration(self, p0, integrand, name, bounds=None, kwargs=None):
        kwargs = kwargs if kwargs else {}

        if bounds is None:
            bounds = (
                (GRID.MIN_MOMENTUM,
                 GRID.MAX_MOMENTUM),
                (lambda p1: GRID.MIN_MOMENTUM,
                 lambda p1: min(p0 + p1, GRID.MAX_MOMENTUM)),
            )

        if self.DETAILED_OUTPUT:
            with benchmark("\t"):
                integral, error = integrators.integrate_2D(
                    lambda p1, p2: integrand(p0, p1, p2, **kwargs),
                    bounds=bounds
                )

                print '{name:}\t\tI( {p0:5.2f} ) = {integral: .5e}\t'\
                    .format(name=name, integral=integral, p0=p0 / UNITS.MeV),

        else:
            integral, error = integrators.integrate_2D(
                lambda p1, p2: integrand(p0, p1, p2, **kwargs),
                bounds=bounds
            )

        return integral, error

    def integrate(self, p0):
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
            h=PARAMS.dy, order=order
        )

        total_integral = (prediction - particle.distribution(p0)) / PARAMS.dy

        return total_integral

    def integral_1(self, p0, bounds=None):
        return self.benchmarked_integration(
            p0, self.integrand,
            kwargs={'F_A': False, 'F_B': False, 'F_1': True, 'F_f': False},
            name=str(self), bounds=bounds
        )

    def integral_f(self, p0, bounds=None):
        return self.benchmarked_integration(
            p0, self.integrand,
            kwargs={'F_A': False, 'F_B': False, 'F_1': False, 'F_f': True},
            name=str(self), bounds=bounds
        )

    def integrand(self, p0, p1, p2, F_A=True, F_B=True, F_1=False, F_f=False):

        """
        Collision integral interior.

        :param F_A, F_B: `bool`. Include/skip naive form terms of the distribution functional
        :param F_1, F_f: `bool`. Include/skip terms of the linearized distribution functional

        WARNING: $F_1 + F_f \neq F_A + F_B = F_1 + f(p_0) F_f$
        """

        p = [p0, p1, p2, 0]
        p, E, m = self.calculate_kinematics(p)
        if not self.in_bounds(p, E, m):
            return 0.

        integrand = self.constant

        ds = 0.
        if p[0] != 0:
            for M in self.Ms:
                ds += D(p=p, E=E, m=m, K1=M.K1, K2=M.K2, order=M.order) / p[0] / E[0]
        else:
            for M in self.Ms:
                ds += Db1(*p[1:]) + m[1] * (E[2] * E[3] + Db2(*p[1:]))

        integrand *= ds

        # Avoid rounding errors and division by zero
        for i in [1, 2, 3]:
            if m[i] != 0:
                integrand *= p[i] / E[i]

        if integrand == 0:
            return 0

        fau = 0
        if F_B:
            fau += self.F_B(p)
        if F_A:
            fau += self.F_A(p)

        # Be on the safe side
        if F_1 and not F_f:
            fau += self.F_1(p)
        elif F_f and not F_1:
            fau += self.F_f(p)
        elif F_1 and F_f:
            raise Exception("F_1 and F_f can't be used at the same time")

        integrand *= fau

        return integrand

    """ === Integration region bounds methods === """

    def in_bounds(self, p, E=None, m=None):
        """ $D$-functions involved in the interactions imply a cut-off region for the collision\
            integrand. In the general case of arbitrary particle masses, this is a set of \
            irrational inequalities that can hardly be solved (at least, Wolfram Mathematica does\
            not succeed in this). To avoid excessive computations, it is convenient to do an early\
            `return 0` when the particles kinematics lay out of the cut-off region """
        if not E or not m:
            p, E, m = self.calculate_kinematics(p)

        q1, q2 = (p[0], p[1]) if p[0] > p[1] else (p[1], p[0])
        q3, q4 = (p[2], p[3]) if p[2] > p[3] else (p[3], p[2])

        is_in = E[3] >= m[3] \
            and q1 <= q2 + q3 + q4 \
            and q3 <= q1 + q2 + q4

        return is_in

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
