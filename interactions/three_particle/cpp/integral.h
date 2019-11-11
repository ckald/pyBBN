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

enum class Kinematics_3 {
  DECAY = 1,
  CREATION = -1
};

enum class CollisionIntegralKind_3 {
  Full = 0,
  F_1 = 1,
  F_f = 2,
  Full_vacuum_decay = 3,
  F_1_vacuum_decay = 4,
  F_f_vacuum_decay = 5,
  F_creation = 6,
  F_decay = 7
};

struct grid_t3 {
    grid_t3(std::vector<dbl> grid, std::vector<dbl> distribution)
        : grid(grid), distribution(distribution) {}
    std::vector<dbl> grid;
    std::vector<dbl> distribution;
};

struct particle_t3 {
    particle_t3(int eta, dbl m, grid_t3 grid, int in_equilibrium, dbl T)
        : eta(eta), m(m), grid(grid), in_equilibrium(in_equilibrium), T(T) {}
    int eta;
    dbl m;
    grid_t3 grid;
    int in_equilibrium;
    dbl T;
};

struct reaction_t3 {
    reaction_t3(particle_t3 specie, int side) : specie(specie), side(side) {}
    particle_t3 specie;
    int side;
};

dbl energy(dbl y, dbl mass);

std::pair<int, int> binary_find(const std::vector<dbl> &grid, dbl x);

dbl distribution_interpolation(const std::vector<dbl> &grid,
                               const std::vector<dbl> &distribution,
                               dbl p, dbl m, int eta, dbl T,
                               bool in_equilibrium);

dbl F_A(const std::vector<reaction_t3> &reaction, const std::array<dbl, 3> &f, int skip_index);
dbl F_B(const std::vector<reaction_t3> &reaction, const std::array<dbl, 3> &f, int skip_index);

dbl F_f(const std::vector<reaction_t3> &reaction, const std::array<dbl, 3> &f);
dbl F_1(const std::vector<reaction_t3> &reaction, const std::array<dbl, 3> &f);

int Sgn(double lam);
int in_bounds(const std::array<dbl, 3> p, const std::array<dbl, 3> E, const std::array<dbl, 3> m);
