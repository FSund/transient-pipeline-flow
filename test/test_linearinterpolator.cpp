#include "debug.hpp"

#include "utilities/linearinterpolator.hpp"

using namespace arma;
using namespace std;

TEST_CASE("LinearInterpolator")
{
    SUBCASE("No interpolation")
    {
        utils::LinearInterpolator interp(vec({0, 1}), vec({0, 1}));

        CHECK(interp.getValueAtPoint(-1) == 0);
        CHECK(interp.getValueAtPoint(0) == 0);
        CHECK(interp.getValueAtPoint(0.25) == 0);
        CHECK(interp.getValueAtPoint(0.5) == 0);
        CHECK(interp.getValueAtPoint(0.75) == 0);
//        CHECK(interp.getValueAtPoint(1.0) == 1);
//        CHECK(interp.getValueAtPoint(2.0) == 1);
    }

    SUBCASE("Linear interpolation")
    {
        utils::LinearInterpolator interp(vec({0, 1}), vec({0, 1}), 1);
        CHECK(interp.getValueAtPoint(-1) == 0);
        CHECK(interp.getValueAtPoint(0) == 0);
        CHECK(interp.getValueAtPoint(0.25) == 0.25);
        CHECK(interp.getValueAtPoint(0.5) == 0.5);
        CHECK(interp.getValueAtPoint(0.75) == 0.75);
        CHECK(interp.getValueAtPoint(1.0) == 1);
        CHECK(interp.getValueAtPoint(2.0) == 1);
    }
}
