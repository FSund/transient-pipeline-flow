#include "debug.hpp"

#include "physics.hpp"
#include "pipeline.hpp"
#include "solver/boundaryconditions.hpp"
#include "equationofstate/gerg04.hpp"
#include "heattransfer/unsteady.hpp"
#include "equationofstate/gerg04.hpp"
#include "equationofstate/equationofstate.hpp"

using namespace arma;
using namespace std;

TEST_SUITE_BEGIN("Physics");

TEST_CASE("Constructors")
{
    CHECK_THROWS(Physics(Pipeline(10), "garbage"));
    CHECK_THROWS(Physics(Pipeline(10), "garbage", "more garbage"));
    CHECK_THROWS(Physics(Pipeline(10), "BWRS", "garbage"));
    CHECK_THROWS(Physics(Pipeline(10), "garbage", "SteadyState"));

    CHECK_THROWS(Physics(Pipeline(10), "bwrss"));
    CHECK_THROWS(Physics(Pipeline(10), "BWRS", "Unsteadyy"));

    for (string eos : {"BWRS", "GERG04", "IdealGas"})
    {
        CHECK_NOTHROW(Physics(Pipeline(10), eos));
        for (string heat : {"SteadyState", "Unsteady"})
        {
            CHECK_NOTHROW(Physics(Pipeline(10), eos, heat));
        }
    }

    SUBCASE("From Config")
    {
        Config config;
        config.equationOfState = "GERG04";
        config.heatTransfer = "Unsteady";

        Physics physics(Pipeline(), config);

        // check type by casting
        CHECK(dynamic_cast<const GERG04*>(&physics.equationOfState().at(0)));
        CHECK(dynamic_cast<const UnsteadyHeatTransfer*>(&physics.heatTransfer().at(0)));
    }
}

TEST_CASE("update derived")
{
    Pipeline gas(10);
    Physics phys(gas);

    // should not touch pressure, temperature, flow or composition
    Pipeline newGas(gas);
    phys.updateDerivedProperties(newGas);
    CHECK(arma::all(newGas.pressure() == gas.pressure()));
    CHECK(arma::all(newGas.temperature() == gas.temperature()));
    CHECK(arma::all(newGas.flow() == gas.flow()));
    for (uword i = 0; i < gas.size(); i++)
    {
        CHECK(arma::all(vec(newGas.composition().at(i)) == vec(gas.composition().at(i))));
    }

    // or these other constants
    CHECK(arma::all(newGas.gridPoints() == gas.gridPoints()));
    CHECK(arma::all(newGas.diameter() == gas.diameter()));
    CHECK(arma::all(newGas.height() == gas.height()));
    CHECK(arma::all(newGas.burialDepth() == gas.burialDepth()));
    CHECK(arma::all(newGas.roughness() == gas.roughness()));
    CHECK(newGas.length() == gas.length());
}

TEST_CASE("initializeHeatTransferState")
{
    Pipeline pipeline;
    Physics physics(pipeline);
    physics.initializeHeatTransferState(pipeline);
    CHECK(pipeline.heatTransferIsInitialized());
    CHECK(pipeline.heatTransferState().size() == pipeline.size());
}

TEST_CASE("update derived props")
{
    Pipeline pipeline;
    Composition comp(arma::zeros<arma::vec>(10) + 1);
    pipeline.updateComposition(comp);
    Config config;
    config.equationOfState = "DummyGas";
    Physics physics(pipeline, config);

    // avoid zero friction, velocity and reynolds number
    pipeline.flow().fill(100);
    pipeline.roughness().fill(1e-5);

    // set these to zero now, and check that they are non-zero after update
    pipeline.specificGasConstant().fill(0);
    pipeline.density().fill(0);
    pipeline.viscosity().fill(0);
    pipeline.reynoldsNumber().fill(0);
    pipeline.velocity().fill(0);
    pipeline.frictionFactor().fill(0);

    physics.updateDerivedProperties(pipeline);

    CHECK(arma::all(pipeline.compressibilityFactor() == 1.0));
    CHECK(arma::all(pipeline.dZdtAtConstantPressure() == 2.0));
    CHECK(arma::all(pipeline.dZdpAtConstantTemperature() == 3.0));
    CHECK(arma::all(pipeline.dZdtAtConstantDensity() == 4.0));
    CHECK(arma::all(pipeline.heatCapacityConstantPressure() == 5.0));
    CHECK(arma::all(pipeline.heatCapacityConstantVolume() == 6.0));
    CHECK(arma::all(pipeline.composition().at(0).vec() == comp.vec()));
    CHECK(pipeline.molarMass()(0) == 1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9 + 10);

    CHECK(arma::all(pipeline.specificGasConstant() > 0));
    CHECK(arma::all(pipeline.density() > 0));
    CHECK(arma::all(pipeline.viscosity() > 0));
    CHECK(arma::all(pipeline.reynoldsNumber() > 0));
    CHECK(arma::all(pipeline.velocity() > 0));
    CHECK(arma::all(pipeline.frictionFactor() > 0));
}

//TEST_CASE("Simulator simulate")
//{
//    auto state = std::make_shared<Pipeline>(10);
//    state->flow().fill(100);
//    state->pressure() = arma::linspace(10e6, 5e6, state->size());
//    state->temperature().fill(273.15 + 5);
//    state->ambientTemperature().fill(273.15 + 10);

//    cout << state->heatTransferIsInitialized() << endl;
//    ClassWithNoName sim(state);
//    const arma::uword dt = 60;
//    cout << state->heatTransferIsInitialized() << endl;

//    vector<BoundaryConditions> bc(100, BoundaryConditions(*state));
//    sim.simulate(dt, bc);

//    cout << sim.state()->flow() << endl;
//    cout << state->flow() << endl;
//}


TEST_SUITE_END();
