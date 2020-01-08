#include "debug.hpp"

#include "advection/batchtracking.hpp"
#include "advection/batchtrackingstate.hpp"
#include "pipeline.hpp"

using namespace arma;
using namespace std;

TEST_SUITE_BEGIN("BatchTracking");

TEST_CASE("Constructors")
{
    const uword n = 10;
    const vec gridPoints = arma::linspace(100, 1000, n); // test with starting point != 0
    const uword nBatches = n - 1;
    const vec x = {1, 2, 3, 4, 5, 6, 7, 8, 9};
//    SUBCASE("1")
//    {
//        BatchTrackingState state(gridPoints, nBatches, x);

//        CHECK(state.batches().at(0).position() == gridPoints(0));
//        CHECK(state.batches().back().position() == gridPoints.tail(2)(0));
//        for (const auto& batch : state.batches())
//        {
//            CHECK(arma::all(batch.concentration() == x));
//        }
//        CHECK(state.batches().size() == nBatches);
//    }
    SUBCASE("2")
    {
        BatchTrackingState state(gridPoints, x, nBatches);

        CHECK(state.batches().at(0).position() == gridPoints(0));
        CHECK(state.batches().back().position() == gridPoints.tail(2)(0));
        for (const auto& batch : state.batches())
        {
            CHECK(arma::all(batch.concentration() == x));
        }
        CHECK(state.batches().size() == nBatches);
    }
    SUBCASE("3")
    {
        Composition c({0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
        vector<Composition> comp(n, c);
        BatchTrackingState state(gridPoints, comp);

        CHECK(state.batches().at(0).position() == gridPoints(0));
        CHECK(state.batches().back().position() == gridPoints.tail(2)(0));
        for (const auto& batch : state.batches())
        {
            CHECK(arma::all(batch.concentration() == c.vec()));
        }
        CHECK(state.batches().size() == nBatches);
    }
}

TEST_CASE("Long time")
{
    const vec gridPoints = {0, 10, 20, 30, 40, 50, 75, 100, 200, 300};
//    BatchTracking bt(gridPoints, 0, zeros<vec>(1) + 1);
    BatchTrackingState state(gridPoints, zeros<vec>(1) + 1, 5);

    const vec velocity = {1,1,1,1,1,1,1,1,1};

    mat compBC;
    compBC << 0.5 << 0.5;

    state = BatchTracking::advect(state, 60, compBC, velocity);

    compBC << 1.0 << 1.0;
    state = BatchTracking::advect(state, 60, compBC, velocity);

    state = BatchTracking::advect(state, 1000, compBC, velocity);
    vector<vec> out = state.sampleToVec();

    for (auto v : out)
    {
        for (auto value : v)
        {
            CHECK(value == 1.0);
        }
    }
}

TEST_CASE("Exact time")
{
    const uword nGridPoints = 11;
    const vec gridPoints = arma::linspace(0, 100, nGridPoints);
//    BatchTracking bt(gridPoints, 0, zeros<vec>(1) + 1);
    BatchTrackingState state(gridPoints, zeros<vec>(1) + 1, 5);

    const vec velocity = arma::zeros<vec>(nGridPoints - 1) + 1;

    mat compBC;
    compBC << 0.5 << 0.5;

    state = BatchTracking::advect(state, 10.0, compBC, velocity);
    auto out = state.sampleToVec();
    CHECK(out.at(0)(0) == 0.5);
    CHECK(out.at(1)(0) == 1.0);
    CHECK(out.at(2)(0) == 1.0);

    compBC << 1.0 << 0.5;
    state = BatchTracking::advect(state, 10.0, compBC, velocity);
    out = state.sampleToVec();
    CHECK(out.at(0)(0) == 1.0);
    CHECK(out.at(1)(0) == 0.5);
    CHECK(out.at(2)(0) == 1.0);

    state = BatchTracking::advect(state, 79.0, compBC, velocity);
    out = state.sampleToVec();
    CHECK(out.at(nGridPoints - 2)(0) == 0.5);
    CHECK(out.at(nGridPoints - 1)(0) == 1.0);

    state = BatchTracking::advect(state, 10.0, compBC, velocity);
    out = state.sampleToVec();
    CHECK(out.at(nGridPoints - 2)(0) == 1.0);
    CHECK(out.at(nGridPoints - 1)(0) == 0.5);
}

TEST_CASE("throws")
{
    const uword nGridPoints = 11;
    const vec gridPoints = arma::linspace(55, 155, nGridPoints);
    BatchTrackingState state(gridPoints);

    CHECK_THROWS(state.sample(vec({54.9})));
    CHECK_THROWS(state.sample(vec({155.1})));
    CHECK_THROWS(state.sampleInternal(vec({54.9})));
    CHECK_THROWS(state.sampleInternal(vec({155.1})));
    CHECK_THROWS(state.sampleInternal(vec({0})));
    CHECK_THROWS(state.sampleInternal(vec({-50})));
    CHECK_THROWS(state.sampleInternal(vec({200})));

    CHECK_NOTHROW(state.sampleInternal(vec({55})));
    CHECK_NOTHROW(state.sampleInternal(vec({100})));
    CHECK_NOTHROW(state.sampleInternal(vec({155})));
}

TEST_CASE("gridpoints start at != 0")
{
    const uword nGridPoints = 11;
    const vec gridPoints = arma::linspace(55, 155, nGridPoints);

    const vec x = zeros<vec>(1) + 1; // initial concentration
    BatchTrackingState state(gridPoints, x, 5);

    const vec velocity = arma::zeros<vec>(nGridPoints - 1) + 2;

    mat compBC;
    compBC << 2 << 2;

    state = BatchTracking::advect(state, 5.0, compBC, velocity);
    auto out = state.sampleToVec();
    CHECK(out.at(0)(0) == 2.0);
    CHECK(out.at(1)(0) == 1.0);
    CHECK(out.at(2)(0) == 1.0);

    out = state.sampleToVec(vec({55, 64.9, 64.9999999, 65, 65.00000001, 155}));
    CHECK(out.at(0)(0) == 2.0);
    CHECK(out.at(1)(0) == 2.0);
    CHECK(out.at(2)(0) == 2.0);
    CHECK(out.at(3)(0) == 1.0);
    CHECK(out.at(4)(0) == 1.0);
    CHECK(out.at(5)(0) == 1.0);
}

TEST_CASE("zero velocity at inlet")
{
    const uword nGridPoints = 11;
    const vec gridPoints = arma::linspace(55, 155, nGridPoints);

    const vec x = zeros<vec>(1) + 1; // initial concentration
    BatchTrackingState state(gridPoints, x, 5);

    const vec velocity = arma::zeros<vec>(nGridPoints - 1);

    mat compBC;
    compBC << 2 << 2;

    state = BatchTracking::advect(state, 5.0, compBC, velocity);
    vector<vec> out = state.sampleToVec();

    CHECK(state.batches().at(0).position() == gridPoints(0));
    CHECK(state.batches().at(0).position() < state.batches().at(1).position());
}

TEST_SUITE_END();
