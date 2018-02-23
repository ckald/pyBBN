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
        \mathcal{F} &= (1 \pm f_1)(1 \pm f_2) f_3 - f_1 f_2 (1 \pm f_3)
        \\\\ &= \mathcal{F}_B + \mathcal{F}_A
    \end{align}
*/

dbl F_A(const std::vector<reaction_t3> &reaction, const std::array<dbl, 3> &f, int skip_index=-1) {
    /*
    Forward reaction distribution functional term

    \begin{equation}
        \mathcal{F}_A = - f_1 f_2 (1 \pm f_3)
    \end{equation}

    :param skip_index: Particle to skip in the expression
    */

    dbl temp(-1.);

    for (int i = 0; i < 3; ++i) {
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

dbl F_B(const std::vector<reaction_t3> &reaction, const std::array<dbl, 3> &f, int skip_index=-1) {
    /*
    Backward reaction distribution functional term

    \begin{equation}
        \mathcal{F}_B = f_3 (1 \pm f_1) (1 \pm f_2)
    \end{equation}

    :param skip_index: Particle to skip in the expression
    */

    dbl temp(1.);

    for (int i = 0; i < 3; ++i) {
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
    \mathcal{F}(f) = f_3 (1 \pm f_1) (1 \pm f_2) - f_1 f_2 (1 \pm f_3)
\end{equation}

\begin{equation}
    \mathcal{F}(f) = f_1 (\mp f_3 (1 \pm f_2) - f_2 (1 \pm f_3) \
    + f_3 (1 \pm f_2)
\end{equation}

\begin{equation}
    \mathcal{F}(f) = \mathcal{F}_B^{(1)} + f_1 (\mathcal{F}_A^{(1)} \pm_1 \mathcal{F}_B^{(1)})
\end{equation}

$^{(i)}$ in $\mathcal{F}^{(i)}$ means that the distribution function $f_i$ was omitted in the\
corresponding expression. $\pm_j$ represents the $\eta$ value of the particle $j$.
*/
dbl F_f(const std::vector<reaction_t3> &reaction, const std::array<dbl, 3> &f) {
    /* Variable part of the distribution functional */
    return -1;//F_A(reaction, f, 0); //- reaction[0].specie.eta * F_B(reaction, f, 0); // For the decay test
}

dbl F_1(const std::vector<reaction_t3> &reaction, const std::array<dbl, 3> &f) {
    /* Constant part of the distribution functional */
    return 0; //F_B(reaction, f, 0); // For the decay test
}

int Sgn(double lam) {
  return (lam > 0) - (lam < 0);
}

int in_bounds(const std::array<dbl, 3> p, const std::array<dbl, 3> E, const std::array<dbl, 3> m) {
    /* $D$-functions involved in the interactions imply a cut-off region for the collision\
        integrand. In the general case of arbitrary particle masses, this is a set of \
        irrational inequalities that can hardly be solved (at least, Wolfram Mathematica does\
        not succeed in this). To avoid excessive computations, it is convenient to do an early\
        `return 0` when the particles kinematics lay out of the cut-off region */
    dbl q1, q2, q3;

    std::array<dbl, 3> mom = p;

    for (int i=0; i<3; i++) {
        for (int j=i; j<3; j++) {
            if (mom[j] > mom[i]) {
                dbl tmp = mom[i];
                mom[i] = mom[j];
                mom[j] = tmp;
            }
        }
    }

    q1 = mom[0];
    q2 = mom[1];
    q3 = mom[2];

    return (Sgn(q1+q2-q3) + Sgn(q1-q2+q3) - Sgn(q1-q2-q3) - Sgn(q1+q2+q3));
}


std::pair<npdbl, npdbl> integrand_3(
    dbl p0, npdbl &p1_buffer, dbl Ms,
    const std::vector<reaction_t3> &reaction
) {
    /*
    Collision integral interior.
    */
    auto p1s = p1_buffer.unchecked<1>();

    if (p1s.ndim() != 1) {
        throw std::runtime_error("p1s must be 1-dimensional");
    }

    size_t length = p1s.size();

    auto integrands_1_buffer = npdbl(length),
         integrands_f_buffer = npdbl(length);

    auto integrands_1 = integrands_1_buffer.mutable_unchecked<1>(),
         integrands_f = integrands_f_buffer.mutable_unchecked<1>();

    std::array<dbl, 3> m;
    std::array<int, 3> sides;

    for (int i = 0; i < 3; ++i) {
        sides[i] = reaction[i].side;
        m[i] = reaction[i].specie.m;
    }

    if (p0 == 0) {

        std::array<dbl, 3> p, E;

        p[0] = p0;
        p[1] = sqrt(pow(pow(m[0],2) - pow(m[1],2) - pow(m[2],2), 2) - 4. * pow(m[1],2) * pow(m[2],2)) / (2. * m[0]);
        p[2] = 0.;
        E[2] = 0.;
        for (int j = 0; j < 2; ++j) {
            E[j] = energy(p[j], m[j]);
            E[2] += sides[j] * E[j];
        }

        E[2] *= -sides[2];

        if (E[2] < m[2]) {return std::make_pair(p1_buffer, p1_buffer); } // Don't really need this

        p[2] = sqrt(pow(E[2], 2) - pow(m[2], 2));

        dbl temp = 1.;

        temp *= p[1] * Ms;

        if (temp == 0.) {return std::make_pair(p1_buffer, p1_buffer); }

        std::array<dbl, 3> f;
        // The distribution function of the main particle is not required here
        f[0] = -1;
        for (int k = 1; k < 3; ++k) {
            const particle_t3 &specie = reaction[k].specie;
            f[k] = distribution_interpolation(
                specie.grid.grid, specie.grid.distribution,
                p[k],
                specie.m, specie.eta,
                specie.T, specie.in_equilibrium
            );
        }

        integrands_1(0) = temp * F_1(reaction, f);
        integrands_f(0) = temp * F_f(reaction, f);

        return std::make_pair(integrands_1_buffer, integrands_f_buffer);
    }

    else {
        #pragma omp parallel for default(none) shared(std::cout, length, p0, p1s, m, sides, Ms, reaction, integrands_1, integrands_f)
        for (size_t i = 0; i < length; ++i) {
            dbl p1 = p1s(i);

            integrands_1(i) = 0.;
            integrands_f(i) = 0.;

            std::array<dbl, 3> p, E;
            p[0] = p0;
            p[1] = p1;
            p[2] = 0.;
            E[2] = 0.;
            for (int j = 0; j < 2; ++j) {
                E[j] = energy(p[j], m[j]);
                E[2] += sides[j] * E[j];
            }

            E[2] *= -sides[2];

            if (E[2] < m[2]) {continue; }

            p[2] = sqrt(pow(E[2], 2) - pow(m[2], 2));

            dbl temp = 1.;

            temp *= in_bounds(p, E, m);

            if (temp == 0. ) {continue; }

//            std::cout << p[0] << "\t" << p[1]  << "\t"  << p[2]  << "\t"  << E[0] << "\t" << E[1] << "\t" << E[2] << "\t" << m[2] << "\t" << in_bounds(p, E, m) << "\n";

            // Avoid rounding errors and division by zero
            if (m[1] != 0.) { temp *= p[1] / E[1]; }

            if (temp == 0.) {continue; }

            temp *= Ms / p[0] / E[0];

            if (temp == 0.) {continue; } // Can get rid of this later

            std::array<dbl, 3> f;
            // The distribution function of the main particle is not required here
            f[0] = -1;
            for (int k = 1; k < 3; ++k) {
                const particle_t3 &specie = reaction[k].specie;
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
}


PYBIND11_MODULE(integral, m) {
    m.def("distribution_interpolation", &distribution_interpolation,
          "Exponential interpolation of distribution function",
          "grid"_a, "distribution"_a,
          "p"_a, "m"_a=0, "eta"_a=1, "T"_a=1., "in_equilibrium"_a=false);
    m.def("binary_find", &binary_find,
          "grid"_a, "x"_a);

    m.def("integrand_3", &integrand_3,
          "p0"_a, "p1s"_a, "Ms"_a,
          "reaction"_a);

    py::class_<grid_t3>(m, "grid_t3")
        .def(py::init<std::vector<dbl>, std::vector<dbl>>(),
             "grid"_a, "distribution"_a);

    py::class_<particle_t3>(m, "particle_t3")
        .def(py::init<int, dbl, grid_t3, int, dbl>(),
             "eta"_a, "m"_a, "grid"_a, "in_equilibrium"_a, "T"_a);

    py::class_<reaction_t3>(m, "reaction_t3")
        .def(py::init<particle_t3, int>(),
             "specie"_a, "side"_a);
}
