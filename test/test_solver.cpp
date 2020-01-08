#include "debug.hpp"

#include "solver/solver.hpp"
#include "solver/boundaryconditions.hpp"
#include "physics.hpp"
#include "solver/discretizer/enthalpy.hpp"
#include "solver/discretizer/internalenergy.hpp"
#include "pipeline.hpp"
#include "solver/governingequationsolver.hpp"
#include "heattransfer/unsteady.hpp"

using namespace arma;
using std::cout;
using std::endl;

TEST_CASE("Constructors")
{
    SUBCASE("From config")
    {
        Config config;
        config.discretizer = "Enthalpy";
        config.relaxationFactors = {0.3, 0.4, 0.5};
        config.tolerances = {0.001, 0.002, 0.003};
        config.toleranceType = "absolute";

        Solver solver(10, config);

        CHECK(solver.relaxationFactors()(0) == 0.3);
        CHECK(solver.relaxationFactors()(1) == 0.4);
        CHECK(solver.relaxationFactors()(2) == 0.5);

        CHECK(solver.tolerances()(0) == 0.001);
        CHECK(solver.tolerances()(1) == 0.002);
        CHECK(solver.tolerances()(2) == 0.003);

        CHECK(solver.toleranceType() == "absolute");

        CHECK(dynamic_cast<const GoverningEquationSolver<EnthalpyDiscretizer>*>(&solver.governingEquationSolver()));
        CHECK_FALSE(dynamic_cast<const GoverningEquationSolver<InternalEnergyDiscretizer>*>(&solver.governingEquationSolver()));
    }
}

TEST_CASE("batch tracking evaluation correctly implemented")
{
    // check that the iteration stuff handles batch tracking correctly
    const arma::uword n = 11;
    const double length = 1e3;

    Pipeline pipeline(n, length);
//    pipeline.pressure().fill(1e6);
    pipeline.pressure() = arma::linspace(1000100, 1e6, n);
    pipeline.flow().fill(77); // approx 10 m/s
    pipeline.roughness().fill(0); // remove friction -- no pressure drop
    pipeline.constantComposition() = false;

    pipeline.initializeBatchTracking();

    Physics physics(pipeline, "BWRS", "FixedQValue");

    physics.updateDerivedProperties(pipeline);
    physics.initializeHeatTransferState(pipeline);
    physics.thermalizeHeatTransfer(pipeline);

    // use low tolerances to force multiple iterations
    const arma::vec tolerances {0.00001, 0.00001, 0.00001};
    Solver solver(n, "InternalEnergy", {1, 1, 2/3.0}, "relative", tolerances);

    const arma::uword dt = 10;
    const Pipeline current = pipeline; // make copy
    BoundaryConditions bc(pipeline);
    bc.inletComposition() = Composition({94, 7, 0.5, 0.03, 0.02, 0.002, 0.001, 0.01, 0.7, 2}).normalized();
//    cout << "C1 @ " << current.gridPoints()(0) << " = " << bc.inletComposition()(0) << endl;

    Pipeline out = solver.solveWithIterations(dt, current, bc, physics);
//    cout << "velocity = " << out.velocity().t() << endl;
//    for (arma::uword i = 0; i < current.gridPoints().size(); i++)
//    {
//        cout << "C1 @ " << current.gridPoints()(i) << " = " << out.composition().at(i)(0) << endl;
//    }

    // have 1000 m pipeline, 11 grid points, so 100 m between each grid point
    // at just over m/s for 10 seconds this should mean that the second grid
    // point has the inlet composition, while the third grid point (and all
    // the ones after that) still has the original composition
    CHECK(out.composition().at(0) == bc.inletComposition());
    CHECK(out.composition().at(1) == bc.inletComposition());
    for (uword i = 2; i < out.size(); i++)
    {
        CHECK(out.composition().at(i) == pipeline.composition().at(0));
    }

    // can also chech that the second batch has not moved further than approx.
    // 100 m
    CHECK(out.batchTrackingState().batches().at(1).position() < 110);
}

TEST_CASE("heat transfer evaluation correctly implemented")
{
    // when performing iterations it's important that the heat transfer is not
    // evaluated repeatedly (at least when using unsteady heat transfer),
    // since this actually moves the heat transfer state forward in time
    // when it shouldn't be

    // so either reset heat transfer state at the end of each iteration, or
    // make sure to use "current.heatTransferState()" to calculate the new guess

    // but how to test for this..?

    const arma::uword n = 11;
    const double length = 1e3;

    Pipeline pipeline(n, length);
    pipeline.pressure().fill(1e6);
    pipeline.flow().fill(100);
    pipeline.roughness().fill(0); // remove friction -- no pressure drop
    pipeline.constantComposition() = true;
    pipeline.ambientTemperature().fill(273.15 + 4);

    Config config;
    config.equationOfState = "BWRS";
    config.heatTransfer = "Unsteady";
    config.bruteForce = true;
    config.maxIterations = 10; // force 10 iterations

    Physics physics(pipeline, config);

    physics.updateDerivedProperties(pipeline);
    physics.initializeHeatTransferState(pipeline);

    pipeline.heatTransferState().back().setTemperature({280, 290, 300});

    Solver solver(n, config);

    const arma::uword dt = 60;
    const Pipeline current = pipeline; // make copy
    BoundaryConditions bc(pipeline);

    Pipeline bruteForce = solver.solveWithIterations(dt, current, bc, physics);

    // the heat transfer state/temperature should be the same after n iterations
    // as performing a single heat transfer evaluation of current.heatTransferState
    // at the new pipeline state

    // heat transfer does not depend on state, so can just use pointer to the
    // current heat transfer object, and call using the state we want
    const HeatTransferBase& heat = physics.heatTransfer().at(n - 1);
    auto check =
            heat.evaluate(current.heatTransferState().at(n - 1), dt, bruteForce.ambientTemperature()(n - 1),
                  bruteForce.pressure()(n - 1), bruteForce.temperature()(n - 1),
                  bruteForce.reynoldsNumber()(n - 1), bruteForce.heatCapacityConstantPressure()(n - 1),
                  bruteForce.viscosity()(n - 1));

    auto heatState = bruteForce.heatTransferState().at(n - 1);
    CHECK(check.heatFlux() == heatState.heatFlux());
    CHECK(equal(check.temperature(), heatState.temperature()));
}

TEST_CASE("Tolerances")
{
    SUBCASE("Same")
    {
        const uword N = 10;
        Pipeline previous(N);
        previous.flow().fill(100);
        previous.temperature().fill(100);
        previous.pressure().fill(100);

        Pipeline guess = previous;
        guess.flow() += 1;
        guess.temperature() += 1;
        guess.pressure() += 1;

        vec tolerances = {0.01, 0.01, 0.01};
        std::string toleranceType = "relative";
        vec relaxationFactors = {1, 1, 1};
        CHECK(Solver::differencesWithinTolerance(guess, previous, tolerances, toleranceType, relaxationFactors));

        tolerances *= 0.9999999;
        CHECK(!Solver::differencesWithinTolerance(guess, previous, tolerances, toleranceType, relaxationFactors));

        tolerances = {0.02, 0.02, 0.02};
        CHECK(Solver::differencesWithinTolerance(guess, previous, tolerances, toleranceType, relaxationFactors));

        tolerances = {0.01, 0.01, 0.01};
        tolerances *= 2;
        relaxationFactors /= 2;
        CHECK(Solver::differencesWithinTolerance(guess, previous, tolerances, toleranceType, relaxationFactors));

        tolerances = {0.01, 0.01, 0.01};
        relaxationFactors /= 2;
        CHECK(!Solver::differencesWithinTolerance(guess, previous, tolerances, toleranceType, relaxationFactors));
    }
    SUBCASE("Different")
    {
        const uword N = 10;
        Pipeline previous(N);
        previous.flow().fill(1e2);
        previous.temperature().fill(1e3);
        previous.pressure().fill(1e4);

        Pipeline guess = previous;
        guess.flow() += 1;
        guess.temperature() += 1;
        guess.pressure() += 1;

        vec relaxationFactors = {1, 1, 1};
        std::string toleranceType = "relative";

        vec tolerances = {0.01, 0.001, 0.001};
        CHECK(Solver::differencesWithinTolerance(guess, previous, tolerances, toleranceType, relaxationFactors));

        tolerances *= 0.99999999;
        CHECK(!Solver::differencesWithinTolerance(guess, previous, tolerances, toleranceType, relaxationFactors));
    }
}

TEST_CASE("Uniform test")
{
    const uword nGridPoints = 10;
    const vec pressure = zeros<vec>(nGridPoints) + 1e6;
    const vec temperature = zeros<vec>(nGridPoints) + 273.15;
    const vec flow = zeros<vec>(nGridPoints) + 100;

    Pipeline gas(nGridPoints);
    gas.pressure() = pressure;
    gas.temperature() = temperature;
    gas.flow() = flow;

    // uniform test, so remove source terms like friction, heat transfer, elevation changes etc.
    gas.roughness().fill(0); // disable friction
    gas.height().fill(0); // disable elevation changes
    gas.diameter().fill(1); // uniform diameter

    Physics physics(gas, "BWRS", "FixedQValue"); // no heat transfer
    physics.updateDerivedProperties(gas); // initialize all derived properties
    physics.initializeHeatTransferState(gas);

    const BoundaryConditions boundaryConditions(gas);
    Solver solver(nGridPoints);
    const arma::uword dt = 3600;

    Pipeline output = solver.solveWithIterations(dt, gas, boundaryConditions, physics);

    CHECK(output.flow()(0) == doctest::Approx(output.flow().tail(1)(0)));
    CHECK(output.temperature()(0) == doctest::Approx(output.temperature().tail(1)(0)));
    CHECK(output.pressure()(0) == doctest::Approx(output.pressure().tail(1)(0)));
}

TEST_CASE("No flow, heat transfer test")
{
//    const uword nGridPoints = 10;
//    const vec pressure = zeros<vec>(nGridPoints) + 1e6;
//    const vec temperature = zeros<vec>(nGridPoints) + 273.15;

//    Pipeline gas(pressure, temperature);

//    gas.roughness().fill(0); // disable friction
//    gas.ambientTemperature() = gas.temperature() + 5; // ENABLE heat transfer
//    gas.height().fill(0); // disable elevation changes
//    gas.diameter().fill(1); // uniform diameter

//    NaturalGas naturalGas(gas);
//    gas = naturalGas.evaluate(gas); // update derived properties

//    Pipeline pipeline(gas);
//    gas = pipeline.evaluate(60, gas); // update friction, heat transfer etc.

//    const BoundaryConditions boundaryConditions(gas);
//    Solver solver(nGridPoints);
//    Pipeline output = solver.solve(3600, gas, boundaryConditions, naturalGas, pipeline);

    // inner wall heat transfer is hardcoded as zero for Re < 4000, so won't
    // get any heat transfer here...
}

TEST_CASE("Heat transfer test" * doctest::skip(true))
{
    const uword nGridPoints = 10;
    const vec pressure = zeros<vec>(nGridPoints) + 1e6;
    const vec temperature = zeros<vec>(nGridPoints) + 273.15;
    const vec flow = zeros<vec>(nGridPoints) + 100;

    Pipeline gas(nGridPoints);
    gas.pressure() = pressure;
    gas.temperature() = temperature;
    gas.flow() = flow;

    gas.roughness().fill(0); // disable friction
    gas.ambientTemperature() = gas.temperature() + 5; // ENABLE heat transfer
    gas.height().fill(0); // disable elevation changes
    gas.diameter().fill(1); // uniform diameter

    Physics physics(gas);
//    physics.initializeHeatTransferState(gas);
    physics.updateDerivedProperties(gas);

    const BoundaryConditions boundaryConditions(gas);
    Solver solver(nGridPoints);
    Pipeline output = solver.solve(360, gas, boundaryConditions, physics);
}

TEST_CASE("Hydrostatic pressure")
{
    SUBCASE("default composition")
    {
        const uword nGridPoints = 10;
        const vec temperature = zeros<vec>(nGridPoints) + 273.15;
        const vec flow = zeros<vec>(nGridPoints);

        // start as close as possible to thermalized state, to avoid issues with
        // pressure and temperature/energy effects
        const vec pressure = linspace(1e6, 1e6 + 7781.79, nGridPoints);

        Pipeline gas(nGridPoints);
        gas.pressure() = pressure;
        gas.temperature() = temperature;
        gas.flow() = flow;
        gas.setLength(100);

        gas.roughness().fill(0); // disable friction
        gas.ambientTemperature() = gas.temperature();
        gas.diameter().fill(1); // uniform diameter

        // hydrostatic pressure
        const double height = -100;
        gas.height() = arma::linspace(0, height, nGridPoints);

        // use ideal gas and no heat transfer
        Physics physics(gas, "IdealGas", "FixedQValue");
        physics.updateDerivedProperties(gas);
        physics.initializeHeatTransferState(gas);

        BoundaryConditions boundaryConditions(gas);
        boundaryConditions.setBoundarySettings({"both", "inlet", "inlet"});

        Solver solver(nGridPoints);

        Pipeline output = gas; // make copy to iterate over
        for (uword i = 0; i < 100; i++)
        {
            output = solver.solve(360, output, boundaryConditions, physics);
        }

        // pressure difference should be approximately equal to the hydrostatic
        // pressure, since we are using ideal gas
        const double dp = output.pressure().tail(1)(0) - output.pressure()(0);
        const double hydrostatic = -mean(output.density())*9.81*height;

        CHECK(dp == doctest::Approx(hydrostatic).epsilon(1e-6));

        // flow should be zero
        CHECK(arma::sum(output.flow()) == Approx(0));
    }

    SUBCASE("only C1")
    {
        const uword nGridPoints = 10;
        const vec temperature = zeros<vec>(nGridPoints) + 273.15;
        const vec flow = zeros<vec>(nGridPoints);
        const Composition x(vec({1, 0, 0, 0, 0, 0, 0, 0, 0, 0}));

        // start as close as possible to thermalized state, to avoid issues with
        // pressure and temperature/energy effects
        const vec pressure = linspace(1e6, 1e6 + 6952.57, nGridPoints);

        Pipeline gas(nGridPoints);
        gas.pressure() = pressure;
        gas.temperature() = temperature;
        gas.flow() = flow;
        gas.updateComposition(x);
        gas.setLength(100);

        gas.roughness().fill(0); // disable friction
        gas.ambientTemperature() = gas.temperature();
        gas.diameter().fill(1); // uniform diameter

        // hydrostatic pressure
        const double height = -100;
        gas.height() = arma::linspace(0, height, nGridPoints);

        // use ideal gas and no heat transfer
        Physics physics(gas, "IdealGas", "FixedQValue");
        physics.updateDerivedProperties(gas);
        physics.initializeHeatTransferState(gas);

        BoundaryConditions boundaryConditions(gas);
        boundaryConditions.setBoundarySettings({"both", "inlet", "inlet"});

        Solver solver(nGridPoints);

        Pipeline output = gas; // make copy to iterate over
        for (uword i = 0; i < 100; i++)
        {
            output = solver.solve(360, output, boundaryConditions, physics);
        }

        // pressure difference should be approximately equal to the hydrostatic
        // pressure, since we are using ideal gas
        const double dp = output.pressure().tail(1)(0) - output.pressure()(0);
        const double hydrostatic = -mean(output.density())*9.81*height;

        CHECK(dp == doctest::Approx(hydrostatic).epsilon(1e-6));

        // flow should be zero
        CHECK(arma::sum(output.flow()) == Approx(0));
    }
}
