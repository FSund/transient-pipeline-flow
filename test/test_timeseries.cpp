#include "debug.hpp"

#include <vector>

#include "timeseries.hpp"
#include "pipeline.hpp"
#include "solver/boundaryconditions.hpp"

using arma::mat;
using arma::vec;
using arma::uvec;
using arma::uword;
using arma::zeros;
using arma::linspace;
using std::vector;
using doctest::Approx;

TEST_CASE("TimeSeries constructors")
{
    SUBCASE("Simple constructor")
    {
        CHECK(TimeSeries(10).size() == 10);
        CHECK(TimeSeries(10).timestamps().size() == 10);

        CHECK(TimeSeries(10).inletFlow().size() == 10);
        CHECK(TimeSeries(10).outletFlow().size() == 10);
        CHECK(TimeSeries(10).inletPressure().size() == 10);
        CHECK(TimeSeries(10).outletPressure().size() == 10);
        CHECK(TimeSeries(10).inletTemperature().size() == 10);
        CHECK(TimeSeries(10).outletTemperature().size() == 10);

        CHECK(TimeSeries(10).inletComposition().size() == 10);
        CHECK(TimeSeries(10).outletComposition().size() == 10);

        TimeSeries ts(10, 120);
        CHECK(ts.timestamps()(1) - ts.timestamps()(0) == 120);
    }

    SUBCASE("from Pipeline")
    {
        const uword n = 10;
        Pipeline state(n);
        const uvec timestamps = linspace<uvec>(0, 60*(n-1), n);
        state.flow().front() = 10;
        state.flow().back() = 20;
        state.pressure().front() = 100;
        state.pressure().back() = 200;
        state.temperature().front() = 1000;
        state.temperature().back() = 2000;
        std::vector<Composition> comp = state.composition();
        comp.front() = Composition(zeros<vec>(10) + 1);
        comp.back() = Composition(zeros<vec>(10) + 2);
        state.setCompositionUnsafe(comp);
        TimeSeries ts(state, 60, n);

        CHECK(ts.inletFlow()(0) == 10);
        CHECK(ts.outletFlow()(0) == 20);
        CHECK(ts.inletPressure()(0) == 100);
        CHECK(ts.outletPressure()(0) == 200);
        CHECK(ts.inletTemperature()(0) == 1000);
        CHECK(ts.outletTemperature()(0) == 2000);
        CHECK(arma::all(ts.inletComposition().at(0).vec() == zeros<vec>(10) + 1));
        CHECK(arma::all(ts.outletComposition().at(0).vec() == zeros<vec>(10) + 2));
    }

    SUBCASE("From dt/timestamps and vector<BoundaryConditions>")
    {
        vector<BoundaryConditions> bc(10, BoundaryConditions(zeros<mat>(3, 2)));
        CHECK_NOTHROW(TimeSeries(60, bc));
        CHECK(TimeSeries(60, bc).size() == 10);

        arma::uvec timestamps = linspace<uvec>(0, 60*9, 10);
        CHECK_NOTHROW(TimeSeries(timestamps, bc));
        CHECK(TimeSeries(timestamps, bc).timestamps().size() == 10);

        timestamps = linspace<uvec>(0, 60*9, 11); // incompatible with bc size
        CHECK_THROWS(TimeSeries(timestamps, bc));
        timestamps = linspace<uvec>(0, 60*9, 9); // incompatible with bc size
        CHECK_THROWS(TimeSeries(timestamps, bc));
    }

    SUBCASE("from csv-file")
    {   
        // no composition
        CHECK_NOTHROW(TimeSeries(std::string(TRANSFLOW_RESOURCE_PATH) + "/examples/bc-no_composition.csv"));
        TimeSeries ts(std::string(TRANSFLOW_RESOURCE_PATH) + "/examples/bc-no_composition.csv");
        CHECK(ts.size() == 3000);
        CHECK(ts.inletFlow().size() == 3000);
        CHECK(ts.inletComposition().size() == 3000);
        CHECK(ts.outletComposition().size() == 3000);
        CHECK(ts.timestamps().n_elem == 3000);

        CHECK(ts.inletFlow()(0) == 384.96);
        CHECK(ts.outletFlow()(0) == 384.94);
        CHECK(ts.inletPressure()(0) == 14190780);
        CHECK(ts.outletPressure()(0) == 8050544);
        CHECK(ts.inletTemperature()(0) == 277.47);
        CHECK(ts.outletTemperature()(0) == 283.20);

        // with composition
        CHECK_NOTHROW(TimeSeries(std::string(TRANSFLOW_RESOURCE_PATH) + "/examples/bc-with_composition.csv"));
        ts = (std::string(TRANSFLOW_RESOURCE_PATH) + "/examples/bc-with_composition.csv");
        CHECK(ts.size() == 3000);
        CHECK(ts.inletFlow().size() == 3000);
        CHECK(ts.inletComposition().size() == 3000);
        CHECK(ts.outletComposition().size() == 3000);
        CHECK(ts.timestamps().n_elem == 3000);

        CHECK(ts.inletFlow()(0) == 384.96);
        CHECK(ts.outletFlow()(0) == 384.94);
        CHECK(ts.inletPressure()(0) == 14190780);
        CHECK(ts.outletPressure()(0) == 8050544);
        CHECK(ts.inletTemperature()(0) == 277.47);
        CHECK(ts.outletTemperature()(0) == 283.20);

        const auto& inlet = ts.inletComposition().at(0);
        CHECK(inlet(0) == Approx(0.922225).epsilon(1e-6));
        CHECK(inlet(1) == Approx(0.044388).epsilon(1e-6));
        CHECK(inlet(2) == Approx(0.004359).epsilon(1e-6));
        CHECK(inlet(3) == Approx(0.002255).epsilon(1e-6));
        CHECK(inlet(4) == Approx(0.000666).epsilon(1e-6));
        CHECK(inlet(5) == Approx(0.000426).epsilon(1e-6));
        CHECK(inlet(6) == Approx(0.000197).epsilon(1e-6));
        CHECK(inlet(7) == Approx(0.000000).epsilon(1e-6));
        CHECK(inlet(8) == Approx(0.013443).epsilon(1e-6));
        CHECK(inlet(9) == Approx(0.012040).epsilon(1e-6));

        const auto& outlet = ts.outletComposition().at(0);
        CHECK(outlet(0) == Approx(0.922225).epsilon(1e-6));
        CHECK(outlet(1) == Approx(0.044388).epsilon(1e-6));
        CHECK(outlet(2) == Approx(0.004359).epsilon(1e-6));
        CHECK(outlet(3) == Approx(0.002255).epsilon(1e-6));
        CHECK(outlet(4) == Approx(0.000666).epsilon(1e-6));
        CHECK(outlet(5) == Approx(0.000426).epsilon(1e-6));
        CHECK(outlet(6) == Approx(0.000197).epsilon(1e-6));
        CHECK(outlet(7) == Approx(0.000000).epsilon(1e-6));
        CHECK(outlet(8) == Approx(0.013443).epsilon(1e-6));
        CHECK(outlet(9) == Approx(0.012040).epsilon(1e-6));
    }

    SUBCASE("from csv-file with boundary settings")
    {
        TimeSeries ts(std::string(TRANSFLOW_RESOURCE_PATH) + "/examples/bc-no_composition.csv", {"outlet", "both", "outlet"});
        CHECK(ts.inletFlow().isActive() == false);
        CHECK(ts.outletFlow().isActive() == true);
        CHECK(ts.inletPressure().isActive() == true);
        CHECK(ts.outletPressure().isActive() == true);
        CHECK(ts.inletTemperature().isActive() == false);
        CHECK(ts.outletTemperature().isActive() == true);
    }

    SUBCASE("From CSV-file with selected rows")
    {
        CHECK_NOTHROW(TimeSeries(std::string(TRANSFLOW_RESOURCE_PATH) + "/examples/bc-with_composition.csv", 0));
        CHECK_NOTHROW(TimeSeries(std::string(TRANSFLOW_RESOURCE_PATH) + "/examples/bc-with_composition.csv", 100));
        CHECK_NOTHROW(TimeSeries(std::string(TRANSFLOW_RESOURCE_PATH) + "/examples/bc-with_composition.csv", 2999));
        CHECK_THROWS(TimeSeries(std::string(TRANSFLOW_RESOURCE_PATH) + "/examples/bc-with_composition.csv", 3000));
        CHECK_NOTHROW(TimeSeries(std::string(TRANSFLOW_RESOURCE_PATH) + "/examples/bc-with_composition.csv", 0, 100));
        CHECK_NOTHROW(TimeSeries(std::string(TRANSFLOW_RESOURCE_PATH) + "/examples/bc-with_composition.csv", 0, 2999));
        CHECK_THROWS(TimeSeries(std::string(TRANSFLOW_RESOURCE_PATH) + "/examples/bc-with_composition.csv", 0, 3000));

        CHECK(TimeSeries(std::string(TRANSFLOW_RESOURCE_PATH) + "/examples/bc-with_composition.csv", 0, 2999).size() == 3000);
        CHECK(TimeSeries(std::string(TRANSFLOW_RESOURCE_PATH) + "/examples/bc-with_composition.csv", 100).size() == 101);
        CHECK(TimeSeries(std::string(TRANSFLOW_RESOURCE_PATH) + "/examples/bc-with_composition.csv", 0, 100).size() == 101);
    }
}

TEST_CASE("user-defined conversion to vector<TimeStep>")
{
    TimeSeries ts(10);
    const arma::vec flow = arma::linspace<vec>(0, 9, 10);
    ts.inletFlow() = flow;
    std::vector<TimeStep> steps(ts);
    for (std::size_t i = 0; i < 10; i++)
    {
        CHECK(steps.at(i).inletFlow() == flow(i));
    }
}

TEST_CASE("setters and getters")
{
    TimeSeries ts(10);

    ts.inletFlow().fill(5);
    CHECK(ts.inletFlow()(0) == 5);

    ts.outletFlow().fill(6);
    CHECK(ts.outletFlow()(0) == 6);

    ts.inletPressure().fill(7);
    CHECK(ts.inletPressure()(0) == 7);

    ts.outletPressure().fill(7);
    CHECK(ts.outletPressure()(0) == 7);

    ts.inletTemperature().fill(8);
    CHECK(ts.inletTemperature()(0) == 8);

    ts.outletTemperature().fill(9);
    CHECK(ts.outletTemperature()(0) == 9);

    ts.inletComposition().at(0) = Composition(zeros<vec>(10) + 1);
    CHECK(ts.inletComposition().at(0)(0) == 1);

    ts.outletComposition().at(0) = Composition(zeros<vec>(10) + 2);
    CHECK(ts.outletComposition().at(0)(0) == 2);
}

typedef TimeSeries::Series Series;
TEST_CASE("TimeSeries::Series")
{
    SUBCASE("constructors")
    {
        // Series(const bool)
        CHECK(Series(true).isActive() == true);
        CHECK(Series(false).isActive() == false);

        // Series(const arma::vec&)
        CHECK(Series(arma::vec({1, 2, 3})).size() == 3);
        CHECK(Series(arma::vec({1, 2, 3}))(0) == 1);
        CHECK(Series(arma::vec({1, 2, 3}))(1) == 2);
        CHECK(Series(arma::vec({1, 2, 3}))(2) == 3);
        CHECK(Series(arma::vec({1, 2, 3})).isActive() == true); // active by default

        // Series(const arma::vec&, const bool)
        CHECK(Series(arma::vec({1, 2, 3}), true).isActive() == true);
        CHECK(Series(arma::vec({1, 2, 3}), false).isActive() == false);
        CHECK(Series(arma::vec({1, 2, 3}), false)(0) == 1);
        CHECK(Series(arma::vec({1, 2, 3}), false)(1) == 2);
        CHECK(Series(arma::vec({1, 2, 3}), false)(2) == 3);
    }

    SUBCASE("setters")
    {
        Series s;
        CHECK(s.isActive() == false);

        s.set(arma::vec({1, 2, 3}));
        CHECK(s.size() == 3);
        CHECK(s.isActive() == true);
        CHECK(s(0) == 1);
        CHECK(s(1) == 2);
        CHECK(s(2) == 3);

        s.set(arma::vec({1, 2, 3}), true);
        CHECK(s.size() == 3);
        CHECK(s.isActive() == true);
        s.set(arma::vec({1, 2, 3}), false);
        CHECK(s.isActive() == false);
        CHECK(s.size() == 3);

        s.fill(4);
        CHECK(s.isActive() == true);
        CHECK(s.size() == 3);
        CHECK(s(0) == 4);
        CHECK(s(1) == 4);
        CHECK(s(2) == 4);

        s.setActive(false);
        CHECK(s.isActive() == false);
        s.setActive(true);
        CHECK(s.isActive() == true);

        s.setActive(false);
        CHECK(s.isActive() == false);
        s = arma::vec({1, 2, 3, 4});
        CHECK(s.isActive() == true);
        CHECK(s.size() == 4);
        CHECK(s(0) == 1);
        CHECK(s(1) == 2);
        CHECK(s(2) == 3);
        CHECK(s(3) == 4);
    }
}
