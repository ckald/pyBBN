// Trapezoidal

dbl trapezium(std::function<dbl(dbl)> f, dbl a, dbl b, size_t n,
              std::map<dbl, dbl> *cache=NULL) {
    std::function<dbl(dbl)> F;
    if (cache) {
        F = [&cache, f](dbl x) -> dbl {
            std::map<dbl, dbl>::iterator it = cache->find(x);
            if (it != cache->end()) {
                return it->second;
            }
            dbl res = f(x);
            cache->insert(std::pair<dbl, dbl>(x, res));
            return res;
        };
    } else {
        F = f;
    }

    dbl h = (b - a) / (n - 1);
    dbl result = (F(a) + F(b)) / 2;
    for (size_t i = 1; i < n - 1; ++i) {
        result += F(a + h * i);
    }
    return h * result;
}


std::pair<dbl, dbl> adaptive_trapezium_helper(std::function<dbl(dbl)> f, dbl a, dbl b, dbl epsabs, dbl epsrel, size_t limit, size_t n=1, std::map<dbl, dbl> *cache=NULL) {
    dbl coarse = trapezium(f, a, b, 8, cache);
    dbl fine = trapezium(f, a, b, 15, cache);
    dbl error = fine - coarse;
    // printf("%e ± %e (%e)\n", fine, error, coarse);

    if (fabs(error) <= epsabs || fabs(error / fine) <= epsrel || n >= limit)
    {
        if (n >= limit) {
            printf("adaptive_trapezium did not converge");
        }
        return std::make_pair(fine, error);
    }

    dbl left_result, left_error;
    dbl right_result, right_error;

    std::tie(left_result, left_error) = adaptive_trapezium_helper(f, a, (a+b)/2, epsabs, epsrel, limit, n+1, cache);
    std::tie(right_result, right_error) = adaptive_trapezium_helper(f, (a+b)/2, b, epsabs, epsrel, limit, n+1, cache);
    // printf("%e ± %e\n", left_result + right_result, left_error + right_error);
    return std::pair<dbl, dbl>(left_result + right_result, left_error + right_error);
}


dbl adaptive_trapezium(std::function<dbl(dbl)> f, dbl a, dbl b, dbl &result, dbl &error, dbl epsabs, dbl epsrel, size_t limit) {
    std::map<dbl, dbl> cache;
    std::tie(result, error) = adaptive_trapezium_helper(f, a, b, epsabs, epsrel, limit, 1, &cache);
    // printf("%e ± %e\n", result, error);
    // cache.clear();
    return result;
}


dbl nonadaptive_trapezium(std::function<dbl(dbl)> f, dbl a, dbl b, dbl &result, dbl &error, dbl epsabs, dbl epsrel, size_t limit) {
    size_t n = 20;
    std::map<dbl, dbl> cache;

    dbl old_result = trapezium(f, a, b, n, &cache);
    bool converged = false;
    while (2 * n - 1 <= limit) {
        n += n - 1;
        result = trapezium(f, a, b, n, &cache);
        error = result - old_result;
        converged = fabs(error) <= epsabs || fabs(error / result) <= epsrel;

        if (converged) {
            break;
        }
        old_result = result;
    }

    if (converged) {
        return 0;
    }
    return n;
}


dbl integrand_1st_integration(
    dbl p2, const void *p
) {
    struct integration_params &params = *(struct integration_params *) p;
    dbl p0 = params.p0;
    dbl p1 = params.p1;
    return integrand_full(p0, p1, p2, *params.reaction, *params.Ms, params.kind);
}

dbl integrand_2nd_integration(
    dbl p1, const void *p
) {
    struct integration_params &params = *(struct integration_params *) p;
    dbl g = params.g;
    dbl h = params.h;
    params.p1 = p1;

    auto f = [params](dbl x) -> dbl {return integrand_1st_integration(x, &params);};
    // return trapezium(f, g, h, params.subdivisions);

    dbl result, error;
    size_t status;
    status = nonadaptive_trapezium(f, g, h, result, error, params.abseps, params.releps, params.subdivisions);
    if (status) {
        printf("1st: adaptive_trapezium: %e ± %e\n", result, error);
    }

    return result;
}