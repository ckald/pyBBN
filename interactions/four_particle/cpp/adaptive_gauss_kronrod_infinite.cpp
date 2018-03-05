// Adaptive Gauss-Kronrod on infinite interval

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

    gsl_function F;
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
    F.function = &integrand_1st_integration;

    gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(params.subdivisions);
    status = gsl_integration_qagiu(&F, g, params.abseps, params.releps, params.subdivisions, w, &result, &error);
    if (status) {
        printf("1st integration result: %e ± %e. %i intervals. %s\n", result, error, w->size, gsl_strerror(status));
    }
    gsl_integration_workspace_free(w);


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
        dbl releps = 1e-1;
        dbl abseps = f / stepsize * releps;

        size_t subdivisions = 100;
        struct integration_params params = {
            p0, 0., 0.,
            &reaction, &Ms,
            a, b, g, h,
            0, releps, abseps,
            subdivisions
        };
        F.params = &params;

        gsl_set_error_handler_off();
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(subdivisions);
        status = gsl_integration_qag(&F, a, b, abseps, releps, subdivisions, GSL_INTEG_GAUSS61, w, &result, &error);
        if (status) {
            printf("2nd integration_1 result: %e ± %e. %i intervals\n", result, error, w->size);
        }
        gsl_integration_workspace_free(w);

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
        w = gsl_integration_workspace_alloc(subdivisions);
        status = gsl_integration_qag(&F, a, b, abseps, releps, subdivisions, GSL_INTEG_GAUSS61, w, &result, &error);
        if (status) {
            printf("2nd integration_f result: %e ± %e. %i intervals\n", result, error, w->size);
        }
        gsl_integration_workspace_free(w);

        integral[i] += result;
    }

    return integral;
}