#include "integral.h"


dbl energy(dbl y, dbl mass=0) {
    return sqrt(y*y + mass*mass);
}


std::pair<int, int> binary_find(const std::vector<dbl> &grid, dbl x) {
    int head(0), tail(grid.size() - 1);
    int middle;

    if (grid[tail] < x) {
        return std::make_pair(tail, -1);
    }
    if (grid[head] > x) {
        return std::make_pair(-1, head);
    }

    while (tail - head > 1) {
        middle = (tail + head) / 2;
        if (grid[middle] > x) {
            tail = middle;
        } else {
            head = middle;
        }
    }

    if (grid[tail] == x) {
        return std::make_pair(tail, tail);
    }

    if (grid[head] == x) {
        return std::make_pair(head, head);
    }

    return std::make_pair(head, tail);
}


dbl distribution_interpolation(const std::vector<dbl> &grid,
                               const std::vector<dbl> &distribution,
                               dbl p, dbl m=0., int eta=1, dbl T=1.,
                               bool in_equilibrium=false) {

    if (in_equilibrium) {
        return 1. / (
           exp(energy(p, m) / T)
           + eta
       );
    }

    int i_lo, i_hi;
    std::tie(i_lo, i_hi) = binary_find(grid, p);
    if(i_lo == -1) {
        throw std::runtime_error("Input momentum is too small for the given grid");
    }

    if(i_hi == -1) {
        return distribution[grid.size()-1]
            * exp((energy(grid[grid.size() -1], m) - energy(p, m)) / T);
    }

    if(i_lo == i_hi) {
        return distribution[i_lo];
    }

    dbl p_lo = grid[i_lo],
        p_hi = grid[i_hi];

    // === Exponential interpolation ===
    dbl E_p  = energy(p, m),
        E_lo = energy(p_lo, m),
        E_hi = energy(p_hi, m);

    /*
    \begin{equation}
        g = \frac{ (E_p - E_{low}) g_{high} + (E_{high} - E_p) g_{low} }\
        { (E_{high} - E_{low}) }
    \end{equation}
    */

    dbl g_hi = distribution[i_hi],
        g_lo = distribution[i_lo];

    g_hi = (1. / g_hi - eta);
    if (g_hi > 0) {
        g_hi = log(g_hi);
    } else {
        return 0.;
    }

    g_lo = (1. / g_lo - eta);
    if (g_lo > 0) {
        g_lo = log(g_lo);
    } else {
        return 0.;
    }

    dbl g = ((E_p - E_lo) * g_hi + (E_hi - E_p) * g_lo) / (E_hi - E_lo);

    g = 1. / (exp(g) + eta);
    if (isnan(g)) {
        return 0.;
    }
    return g;

}


/* ## $\mathcal{F}(f_\alpha)$ functional */

/* ### Naive form

    \begin{align}
        \mathcal{F} &= (1 \pm f_1)(1 \pm f_2) f_3 f_4 - f_1 f_2 (1 \pm f_3) (1 \pm f_4)
        \\\\ &= \mathcal{F}_B + \mathcal{F}_A
    \end{align}
*/

dbl F_A(const std::vector<reaction_t> &reaction, const std::array<dbl, 4> &f, int skip_index=-1) {
    /*
    Forward reaction distribution functional term

    \begin{equation}
        \mathcal{F}_A = - f_1 f_2 (1 \pm f_3) (1 \pm f_4)
    \end{equation}

    :param skip_index: Particle to skip in the expression
    */

    dbl temp(-1.);

    for (int i = 0; i < 4; ++i) {
        if (i != skip_index) {
            if (f[i] < 0) {
                throw std::runtime_error("Negative value of distribution function");
            }
            if (reaction[i].side == -1) {
                temp *= f[i];
            } else {
                temp *= 1. - reaction[i].specie.eta * f[i];
            }
        }
    }

    return temp;
}

dbl F_B(const std::vector<reaction_t> &reaction, const std::array<dbl, 4> &f, int skip_index=-1) {
    /*
    Backward reaction distribution functional term

    \begin{equation}
        \mathcal{F}_B = f_3 f_4 (1 \pm f_1) (1 \pm f_2)
    \end{equation}

    :param skip_index: Particle to skip in the expression
    */

    dbl temp(1.);

    for (int i = 0; i < 4; ++i) {
        if (i != skip_index) {
            if (f[i] < 0) {
                throw std::runtime_error("Negative value of distribution function");
            }
            if (reaction[i].side == 1) {
                temp *= f[i];
            } else {
                temp *= 1. - reaction[i].specie.eta * f[i];
            }
        }
    }

    return temp;
}

/*
### Linearized in $\, f_1$ form

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
*/
dbl F_f(const std::vector<reaction_t> &reaction, const std::array<dbl, 4> &f) {
    /* Variable part of the distribution functional */
    return F_A(reaction, f, 0) - reaction[0].specie.eta * F_B(reaction, f, 0);
}

dbl F_1(const std::vector<reaction_t> &reaction, const std::array<dbl, 4> &f) {
    /* Constant part of the distribution functional */
    return F_B(reaction, f, 0);
}


dbl in_bounds(const std::array<dbl, 4> p, const std::array<dbl, 4> E, const std::array<dbl, 4> m) {
    /* $D$-functions involved in the interactions imply a cut-off region for the collision\
        integrand. In the general case of arbitrary particle masses, this is a set of \
        irrational inequalities that can hardly be solved (at least, Wolfram Mathematica does\
        not succeed in this). To avoid excessive computations, it is convenient to do an early\
        `return 0` when the particles kinematics lay out of the cut-off region */
    dbl q1, q2, q3, q4;
    q1 = p[0];
    q2 = p[1];
    q3 = p[2];
    q4 = p[3];

    if (q1 < q2) { std::swap(q1, q2); }
    if (q3 < q4) { std::swap(q3, q4); }

    return (E[3] >= m[3] && q1 <= q2 + q3 + q4 && q3 <= q1 + q2 + q4);
}


std::pair<npdbl, npdbl> integrand(
    dbl p0, npdbl &p1_buffer, npdbl &p2_buffer,
    const std::vector<reaction_t> &reaction, const std::vector<M_t> &Ms
) {
    /*
    Collision integral interior.
    */

    auto p1s = p1_buffer.unchecked<1>(),
         p2s = p2_buffer.unchecked<1>();

    if (p1s.ndim() != 1 || p2s.ndim() != 1) {
        throw std::runtime_error("p1s and p2s must be 1-dimensional");
    }
    if (p1s.shape(0) != p2s.shape(0)) {
        throw std::runtime_error("p1s and p2s must be of the same size!");
    }
    size_t length = p1s.size();

    auto integrands_1_buffer = npdbl(length),
         integrands_f_buffer = npdbl(length);

    auto integrands_1 = integrands_1_buffer.mutable_unchecked<1>(),
         integrands_f = integrands_f_buffer.mutable_unchecked<1>();

    std::array<dbl, 4> m;
    std::array<int, 4> sides;

    for (int i = 0; i < 4; ++i) {
        sides[i] = reaction[i].side;
        m[i] = reaction[i].specie.m;
    }

    #pragma omp parallel for default(none) shared(length, p0, p1s, p2s, m, sides, Ms, reaction, integrands_1, integrands_f)
    for (size_t i = 0; i < length; ++i) {
        dbl p1 = p1s(i),
            p2 = p2s(i);

        integrands_1(i) = 0.;
        integrands_f(i) = 0.;

        if (p2 > p0 + p1) { continue; }

        std::array<dbl, 4> p, E;
        p[0] = p0;
        p[1] = p1;
        p[2] = p2;
        p[3] = 0.;
        E[3] = 0.;
        for (int j = 0; j < 3; ++j) {
            E[j] = energy(p[j], m[j]);
            E[3] += sides[j] * E[j];
        }

        E[3] *= -sides[3];

        if (E[3] < m[3]) { continue; }

        p[3] = sqrt(pow(E[3], 2) - pow(m[3], 2));

        if (!in_bounds(p, E, m)) { continue; }

        dbl temp = 1.;

        // Avoid rounding errors and division by zero
        for (int k = 1; k < 3; ++k) {
            if (m[k] != 0.) {
                temp *= p[k] / E[k];
            }
        }

        if (temp == 0.) { continue; }

        dbl ds = 0.;
        if (p[0] != 0.) {
            for (const M_t &M : Ms) {
                ds += D(p, E, m, M.K1, M.K2, M.order, sides);
            }
            ds /= p[0] * E[0];
        } else {
            for (const M_t &M : Ms) {
                ds += Db(p, E, m, M.K1, M.K2, M.order, sides);
            }
        }
        temp *= ds;

        if (temp == 0.) { continue; }

        std::array<dbl, 4> f;
        // The distribution function of the main particle is not required here
        f[0] = -1;
        for (int k = 1; k < 4; ++k) {
            const particle_t &specie = reaction[k].specie;
            f[k] = distribution_interpolation(
                specie.grid.grid, specie.grid.distribution,
                p[k],
                specie.m, specie.eta,
                specie.T, specie.in_equilibrium
            );
        }

        integrands_1(i) = temp * F_1(reaction, f);
        integrands_f(i) = temp * F_f(reaction, f);

    }

    return std::make_pair(integrands_1_buffer, integrands_f_buffer);
}


PYBIND11_MODULE(integral, m) {
    m.def("distribution_interpolation", &distribution_interpolation,
          "Exponential interpolation of distribution function",
          "grid"_a, "distribution"_a,
          "p"_a, "m"_a=0, "eta"_a=1, "T"_a=1., "in_equilibrium"_a=false);
    m.def("binary_find", &binary_find,
          "grid"_a, "x"_a);

    m.def("D1", &D1);
    m.def("D2", &D2);
    m.def("D3", &D3);
    m.def("D", &D,
          "p"_a, "E"_a, "m"_a,
          "K1"_a, "K2"_a, "order"_a, "sides"_a);
    m.def("Db", &Db,
          "p"_a, "E"_a, "m"_a,
          "K1"_a, "K2"_a, "order"_a, "sides"_a);
    m.def("integrand", &integrand,
          "p0"_a, "p1s"_a, "p2s"_a,
          "reaction"_a, "Ms"_a);

    py::class_<M_t>(m, "M_t")
        .def(py::init<std::array<int, 4>, dbl, dbl>(),
             "order"_a, "K1"_a=0., "K2"_a=0.);

    py::class_<grid_t>(m, "grid_t")
        .def(py::init<std::vector<dbl>, std::vector<dbl>>(),
             "grid"_a, "distribution"_a);

    py::class_<particle_t>(m, "particle_t")
        .def(py::init<int, dbl, grid_t, int, dbl>(),
             "eta"_a, "m"_a, "grid"_a, "in_equilibrium"_a, "T"_a);

    py::class_<reaction_t>(m, "reaction_t")
        .def(py::init<particle_t, int>(),
             "specie"_a, "side"_a);
}