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


dbl distribution_interpolation(const particle_t3 &specie, dbl p) {
    return distribution_interpolation(
        specie.grid.grid, specie.grid.distribution,
        p,
        specie.m, specie.eta,
        specie.T, specie.in_equilibrium
    );
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
    return F_A(reaction, f, 0) - reaction[0].specie.eta * F_B(reaction, f, 0);
}

dbl F_1(const std::vector<reaction_t3> &reaction, const std::array<dbl, 3> &f) {
    /* Constant part of the distribution functional */
    return F_B(reaction, f, 0);
}

dbl F_f_vacuum_decay(const std::vector<reaction_t3> &reaction, const std::array<dbl, 3> &f) {
    /* Variable part of the distribution functional */
    return -1;
}

dbl F_1_vacuum_decay(const std::vector<reaction_t3> &reaction, const std::array<dbl, 3> &f) {
    /* Constant part of the distribution functional */
    return f[2];
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

    std::sort(mom.begin(), mom.end(), std::greater<dbl>());

    q1 = mom[0];
    q2 = mom[1];
    q3 = mom[2];

    return (Sgn(q1+q2-q3) + Sgn(q1-q2+q3) - Sgn(q1-q2-q3) - Sgn(q1+q2+q3));
}


dbl integrand_full(
    dbl p0, dbl p1, int kind,
    const std::vector<reaction_t3> &reaction
) {
    /*
    Collision integral interior.
    */

    dbl integrand = 0.;

    std::array<dbl, 3> m;
    std::array<int, 3> sides;

    for (int i = 0; i < 3; ++i) {
        sides[i] = reaction[i].side;
        m[i] = reaction[i].specie.m;
    }

    std::array<dbl, 3> p, E;

    dbl temp = 1.;

    if (p0 == 0) {

        p[0] = p0;
        p[1] = p1;
        p[2] = 0.;
        E[2] = 0.;
        for (int j = 0; j < 2; ++j) {
            E[j] = energy(p[j], m[j]);
            E[2] += sides[j] * E[j];
        }

        E[2] *= -sides[2];

        if (E[2] < m[2]) {return integrand; }

        p[2] = sqrt(pow(E[2], 2) - pow(m[2], 2));

        temp *= p[1];

        if (temp == 0.) {return integrand; }
    }

    else {

        p[0] = p0;
        p[1] = p1;
        p[2] = 0.;
        E[2] = 0.;
        for (int j = 0; j < 2; ++j) {
            E[j] = energy(p[j], m[j]);
            E[2] += sides[j] * E[j];
        }

        E[2] *= -sides[2];

        if (E[2] < m[2]) {return integrand; }

        p[2] = sqrt(pow(E[2], 2) - pow(m[2], 2));

        temp *= in_bounds(p, E, m);

        if (temp == 0. ) {return integrand; }

        // Avoid rounding errors and division by zero
        if (m[1] != 0.) {
            temp *= p[1] / E[1];
        }

        if (temp == 0.) {return integrand; }

        temp /= p[0] * E[0];

        if (temp == 0.) {return integrand; }

    }

    std::array<dbl, 3> f;
    // The distribution function of the main particle is not required here
    for (int k = 0; k < 3; ++k) {
        const particle_t3 &specie = reaction[k].specie;
        f[k] = distribution_interpolation(specie, p[k]);
    }

    auto integral_kind = CollisionIntegralKind_3(kind);

    switch (integral_kind) {
        case CollisionIntegralKind_3::F_1_vacuum_decay:
             return temp * F_1_vacuum_decay(reaction, f);
        case CollisionIntegralKind_3::F_f_vacuum_decay:
            return temp * F_f_vacuum_decay(reaction, f);
        case CollisionIntegralKind_3::Full_vacuum_decay:
            return temp * (F_1_vacuum_decay(reaction, f) + f[0] * F_f_vacuum_decay(reaction, f));

        case CollisionIntegralKind_3::F_1:
            return temp * F_1(reaction, f);
        case CollisionIntegralKind_3::F_f:
            return temp * F_f(reaction, f);
        case CollisionIntegralKind_3::Full:
        default:
            return temp * (F_1(reaction, f) + f[0] * F_f(reaction, f));
    }
}


Kinematics_3 get_reaction_type(const std::vector<reaction_t3> &reaction) {
    int reaction_type = 0;
    for (const reaction_t3 &reactant : reaction) {
        reaction_type += reactant.side;
    }
    return static_cast<Kinematics_3>(reaction_type);
}


struct integration_params {
    dbl p0;
    dbl p1;
    const std::vector<reaction_t3> *reaction;
    dbl min_1;
    dbl max_1;
    int kind;
    dbl releps;
    dbl abseps;
    size_t subdivisions;
};


dbl integrand_integration(
    dbl p1, void *p
) {
    struct integration_params &params = *(struct integration_params *) p;
    dbl p0 = params.p0;
    return integrand_full(p0, p1, params.kind, *params.reaction);
}


dbl p1_bounds_1(const std::vector<reaction_t3> &reaction
) {
    dbl m0 = reaction[0].specie.m;
    dbl m1 = reaction[1].specie.m;
    dbl m2 = reaction[2].specie.m;
    if (m0 == 0) {return 0.;}
    return sqrt(
                pow(
                    pow(m0, 2) - pow(m1,2) - pow(m2, 2)
                , 2)
                - 4. * pow(m1, 2) * pow(m2, 2)
            )
            / (2. * m0);
}


dbl p1_bounds_2(const std::vector<reaction_t3> &reaction, dbl p0
) {
    dbl m0 = reaction[0].specie.m;
    dbl m1 = reaction[1].specie.m;
    dbl m2 = reaction[2].specie.m;
    if (m0 == 0 && p0 == 0) {return 0.;}
    return (
            pow(m1, 4) + pow(m2, 4)
            -2 * pow(m1, 2) * (pow(m2, 2) + 2 * pow(p0, 2))
        ) / (
            4 * p0 * (pow(m1, 2) - pow(m2, 2))
        );
}


dbl p1_bounds_3(const std::vector<reaction_t3> &reaction, dbl p0,
    int sign1, int sign2
) {
    dbl m0 = reaction[0].specie.m;
    dbl m1 = reaction[1].specie.m;
    dbl m2 = reaction[2].specie.m;

    dbl temp1 = (pow(p0, 2) + pow(m0, 2))
                * (
                    pow(m1, 4) + pow(pow(m2, 2) - pow(m0, 2), 2)
                    - 2 * pow(m1, 2) * (pow(m2, 2) + pow(m0, 2))
                );

        if (temp1 < 0) {
            temp1 = 0;
        }
    dbl temp2 = p0 * (pow(m2, 2) - pow(m1, 2) - pow(m0, 2));

    return  (sign1 * temp2 + sign2 * sqrt(temp1)) / (2 * pow(m0, 2));

}


dbl p1_bounds_4(const std::vector<reaction_t3> &reaction, dbl p0, dbl max_2
) {
    dbl m0 = reaction[0].specie.m;
    dbl m1 = reaction[1].specie.m;
    dbl m2 = reaction[2].specie.m;

    dbl max = energy(max_2, m2) - energy(p0, m0);
    dbl max2 = pow(max, 2) - pow(m1, 2);
    if (max <= 0 || max2 <= 0) {
        return -1;
    }
    return sqrt(max2);
}


std::vector<dbl> integration_3(
    std::vector<dbl> ps, dbl min_1, dbl max_1, dbl max_2, const std::vector<reaction_t3> &reaction,
    dbl stepsize, int kind=0
) {

    std::vector<dbl> integral(ps.size(), 0.);

    auto reaction_type = get_reaction_type(reaction);

    // Note firstprivate() clause: those variables will be copied for each thread
    #pragma omp parallel for default(none) shared(std::cout,ps, reaction, integral, stepsize, kind, reaction_type) firstprivate(min_1, max_1, max_2)
    for (size_t i = 0; i < ps.size(); ++i) {
        dbl p0 = ps[i];

        if (p0 == 0) {
            dbl p1 = p1_bounds_1(reaction);
            integral[i] = integrand_full(p0, p1, kind, reaction);
        }

        else {
            if (reaction_type == Kinematics_3::CREATION) {
                if (reaction[0].specie.m == 0) {
                    dbl min = p1_bounds_2(reaction, p0);
                    min_1 = std::max(min, -min);
                    dbl temp = p1_bounds_4(reaction, p0, max_2);
                    if (temp == -1) {continue; }
                    max_1 = temp;
                    if (min_1 >= max_1) {continue; }
                }
                else {
                    dbl min = p1_bounds_3(reaction, p0, -1, 1);
                    min_1 = std::max(min, -min);
                    dbl temp = p1_bounds_4(reaction, p0, max_2);
                    if (temp == -1) {continue; }
                    max_1 = std::min(p1_bounds_3(reaction, p0, 1, 1), temp);
                    if (min_1 >= max_1) {continue; }
                }
            }

            else {
                dbl min = p1_bounds_3(reaction, p0, 1, 1);
                min_1 = std::max(min, -min);
                max_1 = p1_bounds_3(reaction, p0, -1, 1);
            }

            dbl result(0.), error(0.);
            size_t status;
            gsl_function F;
            F.function = &integrand_integration;

            dbl releps = 1e-2;
            dbl abseps = releps / stepsize;

            size_t subdivisions = 10000;
            gsl_integration_workspace *w = gsl_integration_workspace_alloc(subdivisions);
            struct integration_params params = {
                p0, 0.,
                &reaction,
                min_1, max_1,
                kind, releps, abseps,
                subdivisions
            };
            F.params = &params;

            gsl_set_error_handler_off();
            status = gsl_integration_qag(&F, min_1, max_1, abseps, releps, subdivisions, GSL_INTEG_GAUSS15, w, &result, &error);

            if (status) {
                printf("integration result: %e ± %e. %i intervals. %s\n", result, error, (int) w->size, gsl_strerror(status));
                throw std::runtime_error("Integrator failed to reach required accuracy");
            }

            gsl_integration_workspace_free(w);
            integral[i] += result;
        }
    }
    return integral;
}


PYBIND11_MODULE(integral, m) {
    m.def("distribution_interpolation", [](
        const std::vector<dbl> &grid,
        const std::vector<dbl> &distribution,
        const py::array_t<double> &ps, dbl m=0., int eta=1, dbl T=1.,
        bool in_equilibrium=false) {
            auto v = [grid, distribution, m, eta, T, in_equilibrium](double p) {
                return distribution_interpolation(grid, distribution, p, m, eta, T, in_equilibrium);
            };
            return py::vectorize(v)(ps);
        },
        "Exponential interpolation of distribution function",
        "grid"_a, "distribution"_a,
        "p"_a, "m"_a=0, "eta"_a=1, "T"_a=1., "in_equilibrium"_a=false
    );

    m.def("binary_find", &binary_find,
          "grid"_a, "x"_a);

    m.def("integration_3", &integration_3,
          "ps"_a, "min_1"_a, "max_1"_a, "max_2"_a,
          "reaction"_a, "stepsize"_a, "kind"_a);

    py::enum_<CollisionIntegralKind_3>(m, "CollisionIntegralKind_3")
        .value("Full", CollisionIntegralKind_3::Full)
        .value("F_1", CollisionIntegralKind_3::F_1)
        .value("F_f", CollisionIntegralKind_3::F_f)
        .value("Full_vacuum_decay", CollisionIntegralKind_3::Full_vacuum_decay)
        .value("F_1_vacuum_decay", CollisionIntegralKind_3::F_1_vacuum_decay)
        .value("F_f_vacuum_decay", CollisionIntegralKind_3::F_f_vacuum_decay)
        .enum_::export_values();

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