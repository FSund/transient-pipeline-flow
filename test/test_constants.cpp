#include "debug.hpp"
#include "constants.hpp"
#include "composition.hpp"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062

TEST_SUITE_BEGIN("constants");

TEST_CASE("constants")
{
    CHECK(arma::sum(Composition::defaultComposition.vec()) == 1.0);
    CHECK(constants::pi == PI);
}

TEST_SUITE_END();
