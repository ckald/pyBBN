#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#include <complex>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>


namespace py = pybind11;
using namespace pybind11::literals;
typedef double dbl;
typedef py::array_t<dbl> npdbl;


enum class Kinematics {
  DECAY = 2,
  SCATTERING = 0,
  CREATION = -2
};

struct M_t {
    M_t(std::array<int, 4> order, dbl K1, dbl K2)
        : order(order), K1(K1), K2(K2) {}
    std::array<int, 4> order;
    dbl K1;
    dbl K2;
};

struct grid_t {
    grid_t(std::vector<dbl> grid, std::vector<dbl> distribution)
        : grid(grid), distribution(distribution) {}
    std::vector<dbl> grid;
    std::vector<dbl> distribution;
};

struct particle_t {
    particle_t(int eta, dbl m, grid_t grid, int in_equilibrium, dbl T)
        : eta(eta), m(m), grid(grid), in_equilibrium(in_equilibrium), T(T) {}
    int eta;
    dbl m;
    grid_t grid;
    int in_equilibrium;
    dbl T;
};

struct reaction_t {
    reaction_t(particle_t specie, int side) : specie(specie), side(side) {}
    particle_t specie;
    int side;
};

dbl energy(dbl y, dbl mass);


std::pair<int, int> binary_find(const std::vector<dbl> &grid, dbl x);


dbl distribution_interpolation(const std::vector<dbl> &grid,
                               const std::vector<dbl> &distribution,
                               dbl p, dbl m, int eta, dbl T,
                               bool in_equilibrium);


dbl D1(dbl q1, dbl q2, dbl q3, dbl q4);
dbl D2(dbl q1, dbl q2, dbl q3, dbl q4);
dbl D3(dbl q1, dbl q2, dbl q3, dbl q4);

dbl D(const std::array<dbl, 4> &p, const std::array<dbl, 4> &E, const std::array<dbl, 4> &m,
      dbl K1, dbl K2,
      const std::array<int, 4> &order, const std::array<int, 4> &sides);


dbl Db1(dbl q2, dbl q3, dbl q4);
dbl Db2(dbl q2, dbl q3, dbl q4);

dbl Db(const std::array<dbl, 4> &p, const std::array<dbl, 4> &E, const std::array<dbl, 4> &m,
       dbl K1, dbl K2,
       const std::array<int, 4> &order, const std::array<int, 4> &sides);


dbl F_A(const std::vector<reaction_t> &reaction, const std::array<dbl, 4> &f, int skip_index);
dbl F_B(const std::vector<reaction_t> &reaction, const std::array<dbl, 4> &f, int skip_index);

dbl F_f(const std::vector<reaction_t> &reaction, const std::array<dbl, 4> &f);
dbl F_1(const std::vector<reaction_t> &reaction, const std::array<dbl, 4> &f);

dbl in_bounds(const std::array<dbl, 4> p, const std::array<dbl, 4> E, const std::array<dbl, 4> m);
