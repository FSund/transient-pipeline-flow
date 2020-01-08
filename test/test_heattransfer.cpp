#include "debug.hpp"
#include <iostream>
#include <iomanip>
#include <string>

#include "heattransfer/utils.hpp"
#include "heattransfer/heattransferbase.hpp"
#include "heattransfer/radial.hpp"
#include "heattransfer/fixedqvalue.hpp"
#include "heattransfer/fixeduvalue.hpp"
#include "heattransfer/steadystate.hpp"
#include "heattransfer/unsteady.hpp"
#include "heattransfer/pipewall.hpp"
#include "heattransfer/burialmedium.hpp"
#include "heattransfer/ambientfluid.hpp"

using doctest::Approx;
using namespace std;
using namespace arma;

TEST_CASE("Heat transfer state")
{
    CHECK(HeatTransferState().hasTemperature() == false);
    CHECK(HeatTransferState(0).hasTemperature() == false);
    CHECK(HeatTransferState(0, zeros<vec>(0)).hasTemperature() == true);
    CHECK(HeatTransferState(0, zeros<vec>(3)).hasTemperature() == true);
    CHECK(HeatTransferState(0, zeros<vec>(3)).temperature().size() == 3);
    CHECK_THROWS(HeatTransferState(0).temperature());
    CHECK_NOTHROW(HeatTransferState(0).setTemperature(zeros<vec>(3)));

    HeatTransferState state;
    CHECK(state.hasTemperature() == false);
    state.setTemperature(zeros<vec>(3));
    CHECK(state.hasTemperature() == true);
    CHECK(state.temperature().size() == 3);

    state.setTemperature({1, 2, 3});
    CHECK(state.temperature()(0) == 1);
    CHECK(state.temperature()(1) == 2);
    CHECK(state.temperature()(2) == 3);

    CHECK(HeatTransferState(2.5).heatFlux() == 2.5);
}

TEST_CASE("makeState")
{
    FixedQValue heat(1.0);
    CHECK(heat.makeState(2.0).heatFlux() == 2.0);
    CHECK(heat.makeState(2.0).hasTemperature() == false);
    CHECK(heat.makeState(2.0, 273, 273).hasTemperature() == false);
}

TEST_CASE("Heat transfer utils")
{
    SUBCASE("outer wall film coeff")
    {
        const double outerDiameter = 1.0;

        // properties of fluid outside pipe
        const double heatCapacityConstantPressure = 1000; // [J/kg K]
        const double viscosity = 1.0/1000.0; // [Pa s] = [kg/m*s]
        const double thermalConductivity = 0.1; // [W/m K]
        const double density = 1000; // [kg/m3]
        const double velocity = 0.1; // [m/s]
        double ho = utils::calcOuterWallFilmCoefficient(
            outerDiameter, heatCapacityConstantPressure, viscosity, thermalConductivity, density, velocity);
        CHECK(ho == Approx(61.6165102188));
    }

    SUBCASE("outer wall film coeff")
    {
        AmbientFluid fluid = AmbientFluid::seawater;
        const double dia = 0.9;
        const double ho = utils::calcOuterWallFilmCoefficient(
                    dia, fluid.heatCapacity(), fluid.viscosity(), fluid.conductivity(),
                    fluid.density(), fluid.velocity());
        CHECK(ho == utils::calcOuterWallFilmCoefficient(dia, fluid));
    }

    SUBCASE("inner wall film coeff")
    {
        const double innerDiameter = 1.0;
        const double pressure = 10e6;
        const double heatCapacityConstantPressure = 4000;
        const double viscosityOfGas = 1.0e-5;

        double ReynoldsNumber = 100;
        double hi = utils::calcInnerWallFilmCoefficient(
            innerDiameter, pressure, ReynoldsNumber, heatCapacityConstantPressure, viscosityOfGas);
        CHECK(hi == Approx(0)); // should be 0 when reynolds number is low

        ReynoldsNumber = 5000;
        hi = utils::calcInnerWallFilmCoefficient(
            innerDiameter, pressure, ReynoldsNumber, heatCapacityConstantPressure, viscosityOfGas);
        CHECK(hi == Approx(0.16104));

        ReynoldsNumber = 10e5;
        hi = utils::calcInnerWallFilmCoefficient(
            innerDiameter, pressure, ReynoldsNumber, heatCapacityConstantPressure, viscosityOfGas);
        CHECK(hi == Approx(61.4643705533));
    }

    SUBCASE("equivalent burial layer width")
    {
        // TODO
    }

    SUBCASE("equivalent burial layer radius")
    {
        // TODO
    }

    SUBCASE("log spaced shell widths")
    {
        // TODO
    }

    SUBCASE("equivalent bural layer width")
    {
        // TODO
    }
}

TEST_CASE("Fixed U and Q")
{
    SUBCASE("Fixed U")
    {
        CHECK(FixedUValue(0).evaluateInternal(200, 300).heatFlux() == 0);
        CHECK(FixedUValue(10).evaluateInternal(300, 300).heatFlux() == 0);

        const double U = 5;
        auto heat = FixedUValue(U);
        CHECK(heat.evaluateInternal(300, 300).heatFlux() == 0);
        CHECK(heat.evaluateInternal(300, 290).heatFlux() == 50);
        CHECK(heat.evaluateInternal(290, 300).heatFlux() == -50);
    }

    SUBCASE("Fixed Q")
    {
        CHECK(FixedQValue(0).evaluateInternal().heatFlux() == 0);
        CHECK(FixedQValue(10).evaluateInternal().heatFlux() == 10);
        CHECK(FixedQValue(-10).evaluateInternal().heatFlux() == -10);
        CHECK(FixedQValue(27).evaluateInternal().heatFlux() == 27);
        CHECK(FixedQValue(-27).evaluateInternal().heatFlux() == -27);

        const double Q = 5;
        auto heat = FixedQValue(Q);
        CHECK(heat.evaluateInternal().heatFlux() == 5);
    }
}

TEST_CASE("compare steady and unsteady heat")
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

//    HeatTransferState state(0, arma::zeros<arma::vec>(unsteady.size()) + amb);
    const arma::uword dt = 60;
    auto out1 = unsteady.evaluate(state, dt, amb, pressure, temperature, reyn, cp, visc);
    auto out2 = steady.evaluate(state, dt, amb, pressure, temperature, reyn, cp, visc);

    CHECK(out1.heatFlux() == doctest::Approx(out2.heatFlux()));
}
