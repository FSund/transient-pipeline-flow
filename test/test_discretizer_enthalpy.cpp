#include "debug.hpp"

#include "solver/discretizer/enthalpy.hpp"
#include "constants.hpp"

using arma::uword;
using arma::vec;
using arma::zeros;
using constants::pi;

TEST_CASE("Enthalpy discretizer")
{
    uword nGridPoints = 10;
    EnthalpyDiscretizer s(nGridPoints);
    vec v = zeros<vec>(nGridPoints) + 1; // dummy property
    vec p = v + 1; // pressure, to avoid inf
    vec gridPoints = arma::regspace(0, 1, nGridPoints - 1);
    s.discretizeFromPrimitives(1, v, v, gridPoints, v, v, p, v, v, p, v, v, v, v, v, v, v, v);

    // check grid point 0, all grid points will have same value
    // cont. and mom. eq. tests are the same as for internal energy version, only difference is in energy eq.

    // cont. eq.
    CHECK(s.term_i()  (0, 0, 0) == +4.0/pi); // m
    CHECK(s.term_ipp()(0, 0, 0) == -4.0/pi); // m
    CHECK(s.term_i()  (0, 0, 1) == 1.0/2.0); // p
    CHECK(s.term_ipp()(0, 0, 1) == 1.0/2.0); // p
    CHECK(s.term_i()  (0, 0, 2) == 2.0); // T
    CHECK(s.term_ipp()(0, 0, 2) == 2.0); // T
    CHECK(s.boundaryTerms()(0, 0) == 6.0);

    // mom. eq.
    CHECK(s.term_i()  (0, 1, 0) == 1.0/2.0 + 1.0/(2.0*pi) - 4.0/pi); // m
    CHECK(s.term_ipp()(0, 1, 0) == 1.0/2.0 + 1.0/(2.0*pi) + 4.0/pi); // m
    CHECK(s.term_i()  (0, 1, 1) == -(pi/4.0 + 1.0/pi)); // p
    CHECK(s.term_ipp()(0, 1, 1) == +(pi/4.0 + 1.0/pi)); // p
    CHECK(s.term_i()  (0, 1, 2) == -4.0/pi); // T
    CHECK(s.term_ipp()(0, 1, 2) == +4.0/pi); // T
    CHECK(s.boundaryTerms()(0, 1) == 1.0);

    // en. eq.
    CHECK(s.term_i()  (0, 2, 0) == +4.0/(5.0*pi) - 2.0/(5.0*pi*pi*pi)); // m
    CHECK(s.term_ipp()(0, 2, 0) == -4.0/(5.0*pi) - 2.0/(5.0*pi*pi*pi)); // m
    CHECK(s.term_i()  (0, 2, 1) == +2.0/(5.0*pi)); // p
    CHECK(s.term_ipp()(0, 2, 1) == -2.0/(5.0*pi)); // p
    CHECK(s.term_i()  (0, 2, 2) == 1.0/(2.0*1.0) + ((-1.0/2.0)*(-4.0/5.0))/2.0 - 2.0/(5.0*pi)); // T
    CHECK(s.term_ipp()(0, 2, 2) == 1.0/(2.0*1.0) + ((-1.0/2.0)*(-4.0/5.0))/2.0 + 2.0/(5.0*pi)); // T
    CHECK(s.boundaryTerms()(0, 2) == 1.0);
}
