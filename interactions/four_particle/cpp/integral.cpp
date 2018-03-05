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

dbl F_f_vacuum_decay() {
    /* Variable part of the distribution functional */
    return -1;
}

dbl F_1_vacuum_decay() {
    /* Constant part of the distribution functional */
    return 0;
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


dbl integrand_full(
    dbl p0, dbl p1, dbl p2,
    const std::vector<reaction_t> reaction, const std::vector<M_t> Ms,
    int kind
) {
    /*
    Collision integral interior.
    */

    dbl integrand = 0.;

    std::array<dbl, 4> m;
    std::array<int, 4> sides;

    for (int i = 0; i < 4; ++i) {
        sides[i] = reaction[i].side;
        m[i] = reaction[i].specie.m;
    }

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

    if (E[3] < m[3]) { return integrand; }

    p[3] = sqrt(pow(E[3], 2) - pow(m[3], 2));

    if (!in_bounds(p, E, m)) { return integrand; }

    dbl temp = 1.;

    // Avoid rounding errors and division by zero
    for (int k = 1; k < 3; ++k) {
        if (m[k] != 0.) {
            temp *= p[k] / E[k];
        }
    }

    if (temp == 0.) { return integrand; }

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
        ds /= m[0];
    }
    temp *= ds;

    if (temp == 0.) { return integrand; }

    std::array<dbl, 4> f;
    for (int k = 0; k < 4; ++k) {
        const particle_t &specie = reaction[k].specie;
        f[k] = distribution_interpolation(
            specie.grid.grid, specie.grid.distribution,
            p[k],
            specie.m, specie.eta,
            specie.T, specie.in_equilibrium
        );
    }

    if (kind == 1) {
        return temp * F_1(reaction, f);
    }
    if (kind == 2) {
        return temp * f[0] * F_f(reaction, f);
    }
    if (kind == 3) {
        return temp * F_1_vacuum_decay();
    }
    if (kind == 4) {
        return temp * f[0] * F_f_vacuum_decay();
    }
    if (kind == 5) {
        return temp * (F_1_vacuum_decay() + f[0] * F_f_vacuum_decay());
    }

    return temp * (F_1(reaction, f) + f[0] * F_f(reaction, f));
}


struct integration_params {
    dbl p0;
    dbl p1;
    dbl p2;
    const std::vector<reaction_t> *reaction;
    const std::vector<M_t> *Ms;
    dbl low_1;
    dbl up_1;
    dbl low_2;
    dbl up_2;
    int kind;
    dbl releps;
    dbl abseps;
    size_t subdivisions;
};


// Adaptive Gauss-Kronrod
dbl integrand_1st_integration(
    dbl p2, void *p
) {
    struct integration_params &params = *(struct integration_params *) p;
    dbl p0 = params.p0;
    dbl p1 = params.p1;
    return integrand_full(p0, p1, p2, *params.reaction, *params.Ms, params.kind);
}

dbl y2_min_dec(const std::vector<reaction_t> reaction,
           int i=0, int j=1, int k=2, int l=3) {
    return sqrt(
        (
            pow(
                pow(reaction[i].specie.m - reaction[j].specie.m, 2)
                - pow(reaction[k].specie.m, 2)
                - pow(reaction[l].specie.m, 2)
            , 2)
            - 4 * pow(reaction[k].specie.m * reaction[l].specie.m, 2)
        ) / (4 * pow(reaction[i].specie.m - reaction[j].specie.m, 2))
    );
}

dbl y2_max_dec(const std::vector<reaction_t> reaction,
            dbl p0, dbl p1) {
    return sqrt(
            pow(
                energy(p0, reaction[0].specie.m)
                - energy(p1, reaction[1].specie.m)
                - reaction[3].specie.m
            , 2)
            - pow(reaction[2].specie.m, 2)
        );
}

dbl y2_min_scat(const std::vector<reaction_t> reaction, dbl p0, dbl p1) {
    dbl m0 = reaction[0].specie.m;
    dbl m1 = reaction[1].specie.m;
    dbl m2 = reaction[2].specie.m;
    dbl m3 = reaction[3].specie.m;

    if (p0 != 0){
        return 0.;
    }

    if (m0 != 0) {
        dbl temp1 = m0 + sqrt(pow(m1,2) + pow(p1,2));
        dbl temp2 = pow(m3,2) + pow(p1,2);
        dbl temp3 = ((temp2 - pow(temp1,2) - pow(m2,2)) * p1
                    + sqrt(
                        pow(temp1,2) * (pow(temp1,4) + pow(temp2,2)
                        + (4 * pow(p1,2) - 2 * temp2 + pow(m2,2)) * pow(m2,2)
                        - 2 * pow(temp1,2) * (temp2 + pow(m2,2)))
                    )) / (2 * (pow(temp1,2) - pow(p1,2)));

        return std::max(temp3, 0.);
    }

    if (m1 != 0) {
        dbl temp1 = (-p1 * (pow(m1,2) + pow(m2,2) - pow(m3,2))
                    + sqrt(
                        (pow(p1,2) + pow(m1,2)) * (pow(m1,4)
                        + pow(pow(m2,2)-pow(m3,2),2) - 2 * pow(m1,2)
                        * (pow(m2,2) + pow(m3,2)))
                        )
                    ) / (2 * pow(m1,2));

        dbl temp2 = -temp1;

        return std::max(temp1, temp2);
    }

    return 0.;
}

dbl y2_max_scat(const std::vector<reaction_t> reaction, dbl p0, dbl p1) {
    dbl m0 = reaction[0].specie.m;
    dbl m1 = reaction[1].specie.m;
    dbl m2 = reaction[2].specie.m;
    dbl m3 = reaction[3].specie.m;

    if (p0 != 0) {
        return sqrt(
                pow(
                    energy(p0, reaction[0].specie.m)
                    + energy(p1, reaction[1].specie.m)
                    - reaction[3].specie.m
                , 2)
                - pow(reaction[2].specie.m, 2)
            );
    }

    if (m0 != 0) {
        dbl temp1 = m0 + sqrt(pow(m1,2) + pow(p1,2));
        dbl temp2 = pow(m3,2) + pow(p1,2);
        dbl temp3 = ((-temp2 + pow(temp1,2) + pow(m2,2)) * p1
                    + sqrt(
                        pow(temp1,2) * (pow(temp1,4) + pow(temp2,2)
                        + (4 * pow(p1,2) - 2 * temp2 + pow(m2,2)) * pow(m2,2)
                        - 2 * pow(temp1,2) * (temp2 + pow(m2,2)))
                    )) / (2 * (pow(temp1,2) - pow(p1,2)));

        return std::max(temp3, 0.);
    }

    if (m1 != 0) {
        return (p1 * (pow(m1,2) + pow(m2,2) - pow(m3,2))
                    + sqrt(
                        (pow(p1,2) + pow(m1,2)) * (pow(m1,4)
                        + pow(pow(m2,2)-pow(m3,2),2) - 2 * pow(m1,2)
                        * (pow(m2,2) + pow(m3,2)))
                        )
                    ) / (2 * pow(m1,2));
    }

    return p1;
}


dbl integrand_2nd_integration(
    dbl p1, void *p
) {
    struct integration_params &params = *(struct integration_params *) p;
    dbl low_2 = params.low_2;
    dbl up_2 = params.up_2;
    dbl p0 = params.p0;
    auto reaction = *params.reaction;

    int reaction_type = 0;
    for (const reaction_t &reactant : reaction) {
        reaction_type += reactant.side;
    }

    if (reaction_type == 2) {
        if (p0 == 0) {
            dbl y2_min_1 = y2_min_dec(reaction, 0, 1, 2, 3);
            dbl y2_min_2 = y2_min_dec(reaction, 0, 2, 1, 3);

            low_2 = std::max(y2_min_1 - p1 * (y2_min_1 / y2_min_2), 0.);
        }

        else {
            low_2 = 0;
        }

        up_2 = y2_max_dec(reaction, p0, p1);
    }

    else {
        low_2 = y2_min_scat(reaction, p0, p1);
        up_2 = y2_max_scat(reaction, p0, p1);

    }

    gsl_function F;
    params.p1 = p1;
    params.low_2 = low_2;
    params.up_2 = up_2;
    F.params = &params;

    // return result;
    gsl_set_error_handler_off();

    dbl result, error;
    size_t status;
    F.function = &integrand_1st_integration;

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(params.subdivisions);
    //status = gsl_integration_qags(&F, low_2, up_2, params.abseps, params.releps, params.subdivisions, w, &result, &error);
    status = gsl_integration_qag(&F, low_2, up_2, params.abseps, params.releps, params.subdivisions, GSL_INTEG_GAUSS21, w, &result, &error);
    //if (status) {
    //    printf("(p0=%e, p1=%e) 1st integration result: %e ± %e. %i intervals. %s\n", params.p0, p1, result, error, (int) w->size, gsl_strerror(status));
    //}
    gsl_integration_workspace_free(w);

    return result;
}


std::vector<dbl> integration(
    std::vector<dbl> ps, dbl min_1, dbl max_1, dbl min_2, dbl max_2,
    const std::vector<reaction_t> &reaction,
    const std::vector<M_t> &Ms,
    dbl stepsize, int kind=0
) {

    std::vector<dbl> integral(ps.size(), 0.);

    // Determine the integration bounds
    int reaction_type = 0;
    for (const reaction_t &reactant : reaction) {
        reaction_type += reactant.side;
    }

    #pragma omp parallel for default(none) shared(ps, Ms, reaction, min_1, max_1, min_2, max_2, integral, stepsize, kind, reaction_type)
    for (size_t i = 0; i < ps.size(); ++i) {
        dbl p0 = ps[i];

        dbl low_1(min_1), up_1(max_1), low_2(min_1), up_2(max_1);

        if (reaction_type == 2) {
            dbl max = sqrt(
                pow(energy(p0, reaction[0].specie.m) - reaction[2].specie.m - reaction[3].specie.m, 2)
                - pow(reaction[1].specie.m, 2)
            );
            up_1 = std::min(up_1, max);
        }
        if (reaction_type == 0) {
            dbl min = pow(reaction[2].specie.m + reaction[3].specie.m - energy(p0, reaction[0].specie.m), 2)
                    - pow(reaction[1].specie.m, 2);
            if (min > 0) {
                low_1 = std::max(low_1, sqrt(min));
            }
        }

        dbl result, error;
        size_t status;
        gsl_function F;
        F.function = &integrand_2nd_integration;

        dbl f = distribution_interpolation(
            reaction[0].specie.grid.grid, reaction[0].specie.grid.distribution,
            p0,
            reaction[0].specie.m, reaction[0].specie.eta,
            reaction[0].specie.T, reaction[0].specie.in_equilibrium
        );

        dbl releps = 1e-2;
        //dbl abseps = 1e-10;
        dbl abseps = std::max(f / stepsize * releps, 1e-15);

        size_t subdivisions = 100;
        struct integration_params params = {
            p0, low_1, low_2,
            &reaction, &Ms,
            low_1, up_1, low_2, up_2,
            kind, releps, abseps,
            subdivisions
        };
        F.params = &params;

        gsl_set_error_handler_off();
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(subdivisions);
        //status = gsl_integration_qags(&F, low_1, up_1, abseps, releps, subdivisions, w, &result, &error);
        status = gsl_integration_qag(&F, low_1, up_1, abseps, releps, subdivisions, GSL_INTEG_GAUSS21, w, &result, &error);
        if (status) {
            printf("2nd integration_1 result: %e ± %e. %i intervals. %s\n", result, error, (int) w->size, gsl_strerror(status));
            throw std::runtime_error("Integrator failed to reach required accuracy");
        }
        gsl_integration_workspace_free(w);
        integral[i] += result;

    }

    return integral;
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

    m.def("integration", &integration,
          "ps"_a, "min_1"_a, "max_1"_a, "min_2"_a, "max_2"_a,
          "reaction"_a, "Ms"_a, "stepsize"_a, "kind"_a);

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