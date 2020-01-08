#include "debug.hpp"

#include "pipeline.hpp"
#include "solver/boundaryconditions.hpp"

using namespace arma;

TEST_CASE("Pipeline class")
{
    SUBCASE("constructor 1")
    {
        // check that members are initialized to the correct size
        CHECK(Pipeline(200).gridPoints().n_elem == 200);
        CHECK(Pipeline(200).diameter().n_elem == 200);
        CHECK(Pipeline(200).height().n_elem == 200);
        CHECK(Pipeline(200).roughness().n_elem == 200);

        CHECK(Pipeline(200).burialDepth().n_elem == 200);
        CHECK(Pipeline(200).pipeWall().size() == 200);
        CHECK(Pipeline(200).burialMedium().size() == 200);
        CHECK(Pipeline(200).ambientFluid().size() == 200);

        // check length argument
        CHECK(Pipeline(200, 150).length() == 150);
        CHECK(Pipeline(200, 150).gridPoints()(0) == 0);
        CHECK(Pipeline(200, 150).gridPoints().tail(1)(0) == 150);
    }

    SUBCASE("heat transfer state and copy constructor")
    {
        Pipeline state(10);
        state.heatTransferState().at(0) = HeatTransferState(7);

        Pipeline state2 = state; // copy constructor
        CHECK(state2.heatTransferState().at(0).heatFlux() == doctest::Approx(7));

        state2 = Pipeline(state); // copy constructor
        CHECK(state.heatTransferState().at(0).heatFlux() == doctest::Approx(7));
    }
}

TEST_CASE("setComposition")
{
    Pipeline pipeline(10);
    Composition comp(vec({1, 0, 0, 0, 0, 0, 0, 0, 0, 0}));
    pipeline.updateComposition(comp);

    // check m_composition
    CHECK(pipeline.composition().size() == pipeline.size());
    CHECK(arma::all(pipeline.composition().at(0).vec() == comp.vec()));
    CHECK(arma::all(pipeline.composition().back().vec() == comp.vec()));

    // check batch tracking state
    CHECK(arma::all(
              pipeline.batchTrackingState().batches().front().concentration()
              == comp.vec())
          );

    CHECK(arma::all(
              pipeline.batchTrackingState().sample().front().vec()
              == comp.vec())
          );
}

TEST_CASE("setlength")
{
    Pipeline pipeline(100);
    pipeline.setLength(150);
    CHECK(pipeline.length() == 150);
    CHECK(pipeline.gridPoints()(0) == 0);
    CHECK(pipeline.gridPoints().tail(1)(0) == 150);

    CHECK(pipeline.batchTrackingIsInitialized());
    CHECK_THROWS(pipeline.batchTrackingState().sample(vec({-1})));
    CHECK_THROWS(pipeline.batchTrackingState().sample(vec({151})));
    CHECK_NOTHROW(pipeline.batchTrackingState().sample(vec({0})));
    CHECK_NOTHROW(pipeline.batchTrackingState().sample(vec({150})));
}

TEST_CASE("enable batch tracking")
{
    Pipeline pipeline;
    pipeline.enableBatchTracking();
    CHECK(pipeline.constantComposition() == false);
    CHECK(pipeline.batchTrackingIsInitialized() == true);
}

TEST_CASE("initialize batch tracking")
{
    Pipeline pipeline;
    pipeline.initializeBatchTracking();
    CHECK(pipeline.batchTrackingIsInitialized() == true);
    CHECK(pipeline.batchTrackingState().batches().size() > 0);
    CHECK(arma::all(
              pipeline.batchTrackingState().batches().front().concentration() ==
              pipeline.composition().front().vec())
          );
}

TEST_CASE("get boundary conditions")
{
    Pipeline pipeline;
    auto bc = pipeline.getBoundaryConditions();

    CHECK(bc.inletFlow() == pipeline.flow()(0));
    CHECK(bc.outletFlow() == pipeline.flow().tail(1)(0));
    CHECK(bc.inletPressure() == pipeline.pressure()(0));
    CHECK(bc.outletPressure() == pipeline.pressure().tail(1)(0));
    CHECK(bc.inletTemperature() == pipeline.temperature()(0));
    CHECK(bc.outletTemperature() == pipeline.temperature().tail(1)(0));

    CHECK(arma::all(
              bc.inletComposition().vec()
              == pipeline.composition().front().vec()));
    CHECK(arma::all(
              bc.outletComposition().vec()
              == pipeline.composition().back().vec()));
}

TEST_CASE("setters and getters")
{
    const arma::uword n = 10;
    Pipeline pipeline(n);

    pipeline.gridPoints() = arma::zeros<vec>(n) + 1;
    CHECK(arma::all(pipeline.gridPoints() == arma::zeros<vec>(n) + 1));
    pipeline.diameter() = arma::zeros<vec>(n) + 2;
    CHECK(arma::all(pipeline.diameter() == arma::zeros<vec>(n) + 2));
    pipeline.height() = arma::zeros<vec>(n) + 3;
    CHECK(arma::all(pipeline.height() == arma::zeros<vec>(n) + 3));
    pipeline.roughness() = arma::zeros<vec>(n) + 4;
    CHECK(arma::all(pipeline.roughness() == arma::zeros<vec>(n) + 4));

    pipeline.burialDepth() = arma::zeros<vec>(n) + 5;
    CHECK(arma::all(pipeline.burialDepth() == arma::zeros<vec>(n) + 5));

    // pipeWall
    const PipeWall pipeWall(
        std::vector<PipeWall::Layer>(
            {
                PipeWall::Layer(0.55, Material::soil)
            }
        )
    );
    pipeline.pipeWall() = std::vector<PipeWall>(n, pipeWall);
    CHECK(pipeline.pipeWall().front().layers().size() == 1);
    CHECK(pipeline.pipeWall().front().layer(0).thickness() == 0.55);
    CHECK(pipeline.pipeWall().back().layers().size() == 1);
    CHECK(pipeline.pipeWall().back().layer(0).thickness() == 0.55);

    // burialMedium
    const BurialMedium burialMedium(1, 2, 3);
    pipeline.burialMedium() = std::vector<BurialMedium>(n, burialMedium);
    CHECK(pipeline.burialMedium().front().conductivity() == 1);
    CHECK(pipeline.burialMedium().front().density() == 2);
    CHECK(pipeline.burialMedium().front().heatCapacity() == 3);
    CHECK(pipeline.burialMedium().back().conductivity() == 1);
    CHECK(pipeline.burialMedium().back().density() == 2);
    CHECK(pipeline.burialMedium().back().heatCapacity() == 3);

    // ambientFluid
    const AmbientFluid fluid(1, 2, 3, 4, 5);
    pipeline.ambientFluid() = std::vector<AmbientFluid>(n, fluid);
    CHECK(pipeline.ambientFluid().front().velocity() == 1);
    CHECK(pipeline.ambientFluid().front().viscosity() == 2);
    CHECK(pipeline.ambientFluid().front().conductivity() == 3);
    CHECK(pipeline.ambientFluid().front().density() == 4);
    CHECK(pipeline.ambientFluid().front().heatCapacity() == 5);

    pipeline.flow() = arma::zeros<vec>(n) + 1;
    CHECK(arma::all(pipeline.flow() == arma::zeros<vec>(n) + 1));
    pipeline.pressure() = arma::zeros<vec>(n) + 2;
    CHECK(arma::all(pipeline.pressure() == arma::zeros<vec>(n) + 2));
    pipeline.temperature() = arma::zeros<vec>(n) + 3;
    CHECK(arma::all(pipeline.temperature() == arma::zeros<vec>(n) + 3));
    pipeline.updateComposition(std::vector<Composition>(n, Composition(arma::zeros<vec>(10) + 4)));
    CHECK(arma::all(pipeline.composition().front().vec() == arma::zeros<vec>(10) + 4));
    CHECK(arma::all(pipeline.composition().back().vec() == arma::zeros<vec>(10) + 4));

    pipeline.heatCapacityConstantVolume() = arma::zeros<vec>(n) + 5;
    CHECK(arma::all(pipeline.heatCapacityConstantVolume() == arma::zeros<vec>(n) + 5));
    pipeline.heatCapacityConstantPressure() = arma::zeros<vec>(n) + 6;
    CHECK(arma::all(pipeline.heatCapacityConstantPressure() == arma::zeros<vec>(n) + 6));
    pipeline.density() = arma::zeros<vec>(n) + 7;
    CHECK(arma::all(pipeline.density() == arma::zeros<vec>(n) + 7));
    pipeline.viscosity() = arma::zeros<vec>(n) + 8;
    CHECK(arma::all(pipeline.viscosity() == arma::zeros<vec>(n) + 8));
    pipeline.specificGasConstant() = arma::zeros<vec>(n) + 9;
    CHECK(arma::all(pipeline.specificGasConstant() == arma::zeros<vec>(n) + 9));
    pipeline.molarMass() = arma::zeros<vec>(n) + 1;
    CHECK(arma::all(pipeline.molarMass() == arma::zeros<vec>(n) + 1));

    pipeline.compressibilityFactor() = arma::zeros<vec>(n) + 2;
    CHECK(arma::all(pipeline.compressibilityFactor() == arma::zeros<vec>(n) + 2));
    pipeline.dZdtAtConstantPressure() = arma::zeros<vec>(n) + 3;
    CHECK(arma::all(pipeline.dZdtAtConstantPressure() == arma::zeros<vec>(n) + 3));
    pipeline.dZdpAtConstantTemperature() = arma::zeros<vec>(n) + 4;
    CHECK(arma::all(pipeline.dZdpAtConstantTemperature() == arma::zeros<vec>(n) + 4));
    pipeline.dZdtAtConstantDensity() = arma::zeros<vec>(n) + 5;
    CHECK(arma::all(pipeline.dZdtAtConstantDensity() == arma::zeros<vec>(n) + 5));

    pipeline.velocity() = arma::zeros<vec>(n) + 6;
    CHECK(arma::all(pipeline.velocity() == arma::zeros<vec>(n) + 6));
    pipeline.frictionFactor() = arma::zeros<vec>(n) + 7;
    CHECK(arma::all(pipeline.frictionFactor() == arma::zeros<vec>(n) + 7));
    pipeline.reynoldsNumber() = arma::zeros<vec>(n) + 8;
    CHECK(arma::all(pipeline.reynoldsNumber() == arma::zeros<vec>(n) + 8));

    pipeline.ambientTemperature() = arma::zeros<vec>(n) + 9;
    CHECK(arma::all(pipeline.ambientTemperature() == arma::zeros<vec>(n) + 9));
    pipeline.heatFlow() = arma::zeros<vec>(n) + 1;
    CHECK(arma::all(pipeline.heatFlow() == arma::zeros<vec>(n) + 1));

    pipeline.constantComposition() = true;
    CHECK(pipeline.constantComposition() == true);
    pipeline.constantComposition() = false;
    CHECK(pipeline.constantComposition() == false);

    pipeline.heatTransferIsInitialized() = true;
    CHECK(pipeline.heatTransferIsInitialized() == true);
    pipeline.heatTransferIsInitialized() = false;
    CHECK(pipeline.heatTransferIsInitialized() == false);

    pipeline.batchTrackingIsInitialized() = true;
    CHECK(pipeline.batchTrackingIsInitialized() == true);
    pipeline.batchTrackingIsInitialized() = false;
    CHECK(pipeline.batchTrackingIsInitialized() == false);

    const arma::vec comp {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    const BatchTrackingState newState(pipeline.gridPoints(), comp);
    pipeline.batchTrackingState() = newState;
    CHECK(arma::all(pipeline.batchTrackingState().batches().at(0).concentration() == comp));

    // TODO: test heatTransferState and batchTrackingState
}
