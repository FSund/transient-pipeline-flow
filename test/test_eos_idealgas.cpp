#include "debug.hpp"
#include "equationofstate/idealgas.hpp"

using namespace arma;

TEST_CASE("IdealGas")
{
    IdealGas gas(Composition::defaultComposition);
    vec out = gas.evaluate(1, 1); // independent of everything

    // check ideal gas properties
    CHECK(out(0) == Approx(1)); // ideal gas
    CHECK(out(1) == Approx(0)); // derivatives are zero
    CHECK(out(2) == Approx(0)); // derivatives are zero
    CHECK(out(3) == Approx(0)); // derivatives are zero

    double cp = out(4);
    double cv = out(5);
    CHECK(cp == Approx(5/2*constants::gasConstant)); // c_p of ideal monoatomic gas
    CHECK(cv == Approx(3/2*constants::gasConstant)); // c_v of ideal monoatomic gas
    CHECK(cp - cv == Approx(constants::gasConstant)); // Mayer's relation for ideal gas
}
