#include "debug.hpp"

#include "simulator.hpp"
#include "physics.hpp"
#include "pipeline.hpp"
#include "solver/boundaryconditions.hpp"
#include "solver/discretizer/enthalpy.hpp"
#include "heattransfer/unsteady.hpp"
#include "equationofstate/gerg04.hpp"
#include "constants.hpp"
#include "equationofstate/equationofstate.hpp"

using namespace arma;
using namespace std;

TEST_SUITE_BEGIN("Simulator");

TEST_CASE("Constructors")
{
    SUBCASE("From Pipeline and Config")
    {
        Pipeline pipeline(53, 191e3);
        pipeline.flow().fill(constants::pi);
        Simulator sim(pipeline);
        CHECK(sim.state().flow()(0) == constants::pi);
        CHECK(sim.pipeline().length() == 191e3);
        CHECK(sim.size() == 53);
        CHECK(sim.pipeline().size() == 53);
        CHECK(sim.state().pressure().n_elem == 53);
    }
    SUBCASE("Config checks")
    {
        Config config;
        config.relaxationFactors = {0.3, 0.4, 0.5};
        config.tolerances = {0.001, 0.002, 0.003};
        config.toleranceType = "absolute";
        config.heatTransfer = "Unsteady";
        config.equationOfState = "GERG04";

        Simulator sim(Pipeline(), config);

        CHECK(sim.solver().relaxationFactors()(0) == 0.3);
        CHECK(sim.solver().relaxationFactors()(1) == 0.4);
        CHECK(sim.solver().relaxationFactors()(2) == 0.5);

        CHECK(sim.solver().tolerances()(0) == 0.001);
        CHECK(sim.solver().tolerances()(1) == 0.002);
        CHECK(sim.solver().tolerances()(2) == 0.003);

        CHECK(sim.solver().toleranceType() == "absolute");

        CHECK(dynamic_cast<const GERG04*>(&sim.physics().equationOfState().at(0)));
        CHECK(dynamic_cast<const UnsteadyHeatTransfer*>(&sim.physics().heatTransfer().at(0)));
    }
    SUBCASE("From Physics and Solver")
    {
        Pipeline pipeline(53, 191e3);
        pipeline.flow().fill(constants::pi);

        auto physics = std::make_unique<Physics>(pipeline, "GERG04", "Unsteady");
        std::unique_ptr<Solver> solver = std::make_unique<Solver>(pipeline.size());

        Simulator sim(pipeline, physics, solver);

        CHECK(sim.state().flow()(0) == constants::pi);
        CHECK(sim.pipeline().length() == 191e3);
        CHECK(sim.size() == 53);
        CHECK(sim.pipeline().size() == 53);
        CHECK(sim.state().pressure().n_elem == 53);
        CHECK(sim.physics().size() == 53);

        CHECK(dynamic_cast<const GERG04*>(&sim.physics().equationOfState().at(0)));
        CHECK(dynamic_cast<const UnsteadyHeatTransfer*>(&sim.physics().heatTransfer().at(0)));
    }
    SUBCASE("empty optional sampler")
    {
        Config config;
        config.outputPath = "";

        Simulator sim(Pipeline(), config);
        CHECK_THROWS(sim.sampler());
    }
    SUBCASE("with sampler")
    {
        Config config;
        config.outputPath = "./output/";
        Simulator sim(Pipeline(), config);
        CHECK_NOTHROW(sim.sampler());
        CHECK(sim.sampler().outputDir() == std::filesystem::absolute("./output/"));
    }
}

TEST_CASE("Setters")
{
    SUBCASE("enableBatchTracking()")
    {
        Simulator sim;
        sim.enableBatchTracking();
        CHECK(sim.pipeline().constantComposition() == false);
        CHECK(sim.pipeline().batchTrackingIsInitialized() == true);
    }
}

TEST_CASE("New simulator simulate")
{
    // set up pipeline
    Pipeline pipeline;
    pipeline.flow().fill(100);
    pipeline.pressure() = arma::linspace(10e6, 9.9e6, pipeline.size());
    pipeline.temperature().fill(273.15 + 5);
    pipeline.roughness().fill(5e-7);
    pipeline.ambientTemperature().fill(273.15 + 10);

    // make simulator
    Config config;
    Simulator sim(pipeline, config);

    // set up boundary conditions
    const arma::uword dt = 60;
    const arma::uword nSteps = 100;
    BoundaryConditions bc(pipeline);
    TimeSeries ts(dt, vector<BoundaryConditions>(nSteps, bc));
    ts.setBoundarySettings({"inlet", "outlet", "inlet"});
    sim.simulate(ts);

    CHECK(sim.state().flow()(0) == pipeline.flow()(0));
    CHECK(sim.state().pressure().tail(1)(0) == pipeline.pressure().tail(1)(0));
    CHECK(sim.state().temperature()(0) == pipeline.temperature()(0));
}

TEST_SUITE_END();
