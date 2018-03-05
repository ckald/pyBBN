// Non-adaptive Gauss-Kronrod

dbl integrand_1st_integration(
    dbl p2, void *p
) {
    struct integration_params &params = *(struct integration_params *) p;
    dbl p0 = params.p0;
    dbl p1 = params.p1;
    return integrand_full(p0, p1, p2, *params.reaction, *params.Ms, params.kind);
}

dbl integrand_2nd_integration(
    dbl p1, void *p
) {
    struct integration_params &params = *(struct integration_params *) p;
    dbl g = params.g;
    dbl h = params.h;

    params.p1 = p1;

    gsl_function F;
    F.function = &integrand_1st_integration;
    struct integration_params new_params = {
        params.p0, p1, 0.,
        params.reaction, params.Ms,
        params.a, params.b, params.g, params.h,
        params.kind, params.releps, params.abseps,
        params.subdivisions
    };
    F.params = &new_params;

    dbl result, error;
    size_t neval, status;
    gsl_set_error_handler_off();
    status = gsl_integration_qng(&F, g, h, params.abseps, params.releps, &result, &error, &neval);
    if (status) {
        printf("1st integration_%i result: %e ± %e. %i evaluations. %s\n",
               params.kind,
               result, error, (int) neval, gsl_strerror(status));
    }

    return result;
}


std::vector<dbl> integration(
    std::vector<dbl> ps, dbl a, dbl b, dbl g, dbl h,
    const std::vector<reaction_t> &reaction,
    const std::vector<M_t> &Ms,
    dbl stepsize
) {

    std::vector<dbl> integral(ps.size(), 0.);

    #pragma omp parallel for default(none) shared(ps, Ms, reaction, a, b, g, h, integral, stepsize)
    for (size_t i = 0; i < ps.size(); ++i) {
        dbl p0 = ps[i];

        dbl result, error;
        size_t neval, status;
        gsl_function F;
        F.function = &integrand_2nd_integration;

        dbl f = distribution_interpolation(
            reaction[0].specie.grid.grid, reaction[0].specie.grid.distribution,
            p0,
            reaction[0].specie.m, reaction[0].specie.eta,
            reaction[0].specie.T, reaction[0].specie.in_equilibrium
        );
        dbl releps = 1e-2;
        // dbl abseps = 1e-10;
        dbl abseps = f / stepsize * releps;

        size_t subdivisions = 10;
        struct integration_params params = {
            p0, 0., 0.,
            &reaction, &Ms,
            a, b, g, h,
            0, releps, abseps,
            subdivisions
        };
        F.params = &params;

        gsl_set_error_handler_off();
        status = gsl_integration_qng(&F, a, b, abseps, releps, &result, &error, &neval);

        if (status) {
            printf("2nd integration_1 result: %e ± %e. %i evaluations. %s\n",
                   result, error, (int) neval, gsl_strerror(status));
        }
        integral[i] += result;

        params = {
            p0, 0., 0.,
            &reaction, &Ms,
            a, b, g, h,
            1, releps, abseps,
            subdivisions
        };
        F.params = &params;

        gsl_set_error_handler_off();
        status = gsl_integration_qng(&F, a, b, abseps, releps, &result, &error, &neval);
        if (status) {
            printf("2nd integration_f result: %e ± %e. %i evaluations. %s\n",
                   result, error, (int) neval, gsl_strerror(status));
        }
        integral[i] += result;
    }

    return integral;
}

