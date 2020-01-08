#include "debug.hpp"
#include <iostream>
#include <iomanip>
#include <string>

#include "heattransfer/heattransferbase.hpp"
#include "heattransfer/steadystate.hpp"
#include "heattransfer/unsteady.hpp"
#include "heattransfer/pipewall.hpp"
#include "constants.hpp"

using doctest::Approx;
using namespace std;
using namespace arma;

TEST_SUITE_BEGIN("Examples");

TEST_CASE("no temperature diff")
{
    const double diameter = 0.9;
    const PipeWall wall = PipeWall::defaultPipeWall;
    const double burial = 1.2;
    const BurialMedium medium = BurialMedium::soil;
    const AmbientFluid fluid = AmbientFluid::seawater;

    const double pressure = 1e6;
    const double reyn = 1e5;
    const double cp = 2000;
    const double visc = 1e-5;

    SUBCASE("No temperature difference")
    {
        const double amb = 283.15;
        const double temperature = amb;

        UnsteadyHeatTransfer heat(diameter, wall, burial, medium, fluid);
        auto state = heat.thermalizeToSteadyState(amb, pressure, temperature, reyn, cp, visc);
        CHECK(state.heatFlux() == Approx(0));
        CHECK(arma::mean(arma::abs(state.temperature() - amb)) == Approx(0));
    }

    SUBCASE("No temperature difference 2")
    {
        const double amb = 283.15;
        const double temperature = amb;

        UnsteadyHeatTransfer heat(diameter, wall, burial, medium, fluid);

        // set up state manually
        HeatTransferState state(0, zeros<vec>(heat.size()) + amb);

        state = heat.evaluate(state, 60, amb, pressure, temperature, reyn, cp, visc);
        CHECK(state.heatFlux() == Approx(0));
        CHECK(arma::mean(arma::abs(state.temperature() - amb)) == Approx(0));
    }
}

TEST_CASE("10 degree difference")
{
    const double diameter = 1.0;
    PipeWall wall({PipeWall::Layer(0.5, 1000, 1000, 1000), PipeWall::Layer(0.5, 1000, 1000, 1000)});

    const double burial = -2*diameter; // exposed pipeline
    const BurialMedium medium = BurialMedium::soil;

    // reduce impact of outer film coeff by increasing velocity
    // the original tests had no outer film coeff
    const AmbientFluid fluid(1e5, 1e-5, Material::seawater);

    double ambient = 275;
    double temperature = 285;
    const double pressure = 1e6;
    const double reyn = 1e5;
    const double cp = 1000;
    const double rho = 10;
    const double visc = 1e-5;

    SUBCASE("10 degree difference")
    {
        UnsteadyHeatTransfer heat(diameter, wall, burial, medium, fluid);
        const double answer = -16.2854;
        double q;

        HeatTransferState state(0, zeros<vec>(heat.size()) + 280);
        q = heat.thermalizeToSteadyState(ambient, pressure, temperature, reyn, cp, visc).heatFlux();
        q = -4*q/(diameter*rho);
        CHECK(q == Approx(answer));

        // reverse difference
        ambient = temperature;
        temperature -= 10;

        q = heat.thermalizeToSteadyState(ambient, pressure, temperature, reyn, cp, visc).heatFlux();
        q = -4*q/(diameter*rho);
        CHECK(q == Approx(-answer));
    }
}

TEST_SUITE_END();

TEST_CASE("thermalize")
{
    const double dia = 0.9;
    const double burial = -2*dia;
    UnsteadyHeatTransfer heat1(dia, PipeWall::defaultPipeWall, burial, BurialMedium::soil, AmbientFluid::seawater);
    UnsteadyHeatTransfer heat2(dia, PipeWall::defaultPipeWall, burial, BurialMedium::soil, AmbientFluid::seawater);

    const double amb = 300;
    const double temperature = amb - 12;
    const double pressure = 1e6;
    const double reyn = 1e5;
    const double cp = 2000;
    const double visc = 1e-5;

    // check that manual thermalizing gives same result as thermalize method
    auto state1 = heat1.thermalizeToSteadyState(amb, pressure, temperature, reyn, cp, visc);
    HeatTransferState state2(0, zeros<vec>(heat2.size()));
    state2 = heat1.evaluate(state2, 1e300, amb, pressure, temperature, reyn, cp, visc);
    CHECK(state1.heatFlux() == state2.heatFlux());
    CHECK(arma::all(state1.temperature() == state2.temperature()));

    // check that heat flux after thermalization doesn't change with evaluations
    auto state3 = heat1.evaluate(state2, 60, amb, pressure, temperature, reyn, cp, visc);
    CHECK(state2.heatFlux() == state3.heatFlux());
    CHECK(arma::sum(arma::abs(state2.temperature() - state3.temperature())) == Approx(0));
}

TEST_CASE("Convergence to steady state")
{
    const double dia = 0.9;
    const double burial = 1.2;
    SteadyStateHeatTransfer steady(dia, PipeWall::defaultPipeWall, burial, BurialMedium::soil, AmbientFluid::seawater);
    UnsteadyHeatTransfer unsteady(dia, PipeWall::defaultPipeWall, burial, BurialMedium::soil, AmbientFluid::seawater);

    const double amb = 283.15;
    const double temperature = amb - 5;
    const double pressure = 1e6;
    const double reyn = 1e5;
    const double cp = 2000;
    const double visc = 1e-5;
    auto state = unsteady.thermalizeToSteadyState(amb, pressure, temperature, reyn, cp, visc);

    const arma::uword dt = 60;
    auto out1 = unsteady.evaluate(state, dt, amb, pressure, temperature, reyn, cp, visc);
    auto out2 = steady.evaluate(state, dt, amb, pressure, temperature, reyn, cp, visc);

    CHECK(out1.heatFlux() == doctest::Approx(out2.heatFlux()));

    // check that thermalized state doesn't depend on time step
    out1 = unsteady.evaluate(state, 200, amb, pressure, temperature, reyn, cp, visc);
    CHECK(out1.heatFlux() == doctest::Approx(out2.heatFlux()));
    out1 = unsteady.evaluate(state, 1.5, amb, pressure, temperature, reyn, cp, visc);
    CHECK(out1.heatFlux() == doctest::Approx(out2.heatFlux()));
}

TEST_CASE("Unsteady behaviour")
{
    const double gasDensity = 100;
    const double diameter = 1;
    const double burial = -2*diameter;
    const AmbientFluid fluid(1e5, 1e-5, Material::seawater); // remove effect of outer film coeff

    // defaults have changed, so use these layers
    const PipeWall::Layer steel(0.024, 50, 7800, 500);
    const PipeWall::Layer coating(0.007, 0.74, 1300, 1900);
    const PipeWall::Layer concrete(0.08, 2.9, 2500, 650);
    const PipeWall pipeWall = PipeWall({steel, coating, concrete});

    UnsteadyHeatTransfer heat(diameter, pipeWall, burial, BurialMedium::soil, fluid);

    const double pressure = 10e6;
    const double temperature = 280;
    const double reyn = 10e7;
    const double cp = 3000;
    const double visc = 1e-5;

    // first thermalize to 270 K
    double amb = 270;
    vec shellTemp = heat.thermalizeToSteadyState(amb, pressure, temperature, reyn, cp, visc).temperature();

    // then do a single timestep at 260 K, and check the result
    amb = 260;
    const auto check = heat.evaluateInternal(shellTemp, 60, amb, pressure, temperature, reyn, cp, visc);
    double q = -4*check.heatFlux()/(diameter*gasDensity);
    CHECK(q == Approx(-11.86275));

    // check that multiple evaluations doesn't change the result
    for (uword i = 0; i < 100; i++)
    {
        // DON'T update shell temperature here
        auto result = heat.evaluateInternal(shellTemp, 60, amb, pressure, temperature, reyn, cp, visc);
        CHECK(result.heatFlux() == check.heatFlux());
    }

    // get current temperature
    shellTemp = heat.evaluateInternal(shellTemp, 60, amb, pressure, temperature, reyn, cp, visc).temperature();

    // do another time step at 260 K, and check the new temperature
    auto result = heat.evaluateInternal(shellTemp, 60, amb, pressure, temperature, reyn, cp, visc);
    q = -4*result.heatFlux()/(diameter*gasDensity);
    CHECK(q == Approx(-12.11994));
}
