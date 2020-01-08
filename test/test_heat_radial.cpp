#include "debug.hpp"
#include <iostream>
#include <iomanip>
#include <string>

#include "heattransfer/radial.hpp"
#include "heattransfer/steadystate.hpp"
#include "heattransfer/unsteady.hpp"

using doctest::Approx;
using namespace std;
using namespace arma;

TEST_CASE("makeState")
{
    SteadyStateHeatTransfer heat;
    CHECK(heat.makeState(1.0).heatFlux() == 1.0);
    CHECK(heat.makeState(1.0).hasTemperature() == true);
    CHECK(heat.makeState(1.0).temperature().size() == heat.size());

    CHECK(heat.makeState(1.0, 273, 273).heatFlux() == 1.0);
    CHECK(heat.makeState(1.0, 273, 273).hasTemperature() == true);
    CHECK(heat.makeState(1.0, 273, 273).temperature()(0) == 273);
    CHECK(heat.makeState(1.0, 273, 273).temperature().tail(1)(0) == 273);
    CHECK(heat.makeState(1.0).temperature().size() == heat.size());

    CHECK(heat.makeState(1.0, 200, 300).temperature()(0) == 200);
    CHECK(heat.makeState(1.0, 200, 300).temperature().tail(1)(0) == 300);
}
