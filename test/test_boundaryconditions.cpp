#include "debug.hpp"

#include <string>
#include <map>

#include "solver/boundaryconditions.hpp"
#include "pipeline.hpp"

using arma::mat;
using arma::vec;

typedef BoundaryConditions::SingleCondition SingleCondition;

TEST_CASE("BoundaryConditions")
{
    CHECK_NOTHROW(BoundaryConditions(arma::zeros(3, 2)));

    CHECK_THROWS(BoundaryConditions(arma::zeros(2, 2)));
    CHECK_THROWS(BoundaryConditions(arma::zeros(4, 2)));
    CHECK_THROWS(BoundaryConditions(arma::zeros(3, 0)));
    CHECK_THROWS(BoundaryConditions(arma::zeros(3, 1)));
    CHECK_THROWS(BoundaryConditions(arma::zeros(3, 3)));

    // removed this constructor
//    CHECK_NOTHROW(BoundaryConditions(arma::zeros(13, 2)));

//    CHECK_THROWS(BoundaryConditions(arma::zeros(12, 2)));
//    CHECK_THROWS(BoundaryConditions(arma::zeros(14, 2)));
//    CHECK_THROWS(BoundaryConditions(arma::zeros(13, 0)));
//    CHECK_THROWS(BoundaryConditions(arma::zeros(13, 1)));
//    CHECK_THROWS(BoundaryConditions(arma::zeros(13, 3)));


    CHECK_NOTHROW(BoundaryConditions(arma::zeros<mat>(3, 2), Composition(arma::zeros<vec>(10)), Composition(arma::zeros<vec>(10))));

    CHECK_THROWS(BoundaryConditions(arma::zeros<mat>(4, 2), Composition(arma::zeros<vec>(10)), Composition(arma::zeros<vec>(10))));
    CHECK_THROWS(BoundaryConditions(arma::zeros<mat>(2, 2), Composition(arma::zeros<vec>(10)), Composition(arma::zeros<vec>(10))));
    CHECK_THROWS(BoundaryConditions(arma::zeros<mat>(3, 1), Composition(arma::zeros<vec>(10)), Composition(arma::zeros<vec>(10))));
    CHECK_THROWS(BoundaryConditions(arma::zeros<mat>(4, 1), Composition(arma::zeros<vec>(10)), Composition(arma::zeros<vec>(10))));

    CHECK_THROWS(BoundaryConditions(arma::zeros<mat>(3, 2), Composition(arma::zeros<vec>(11)), Composition(arma::zeros<vec>(10))));
    CHECK_THROWS(BoundaryConditions(arma::zeros<mat>(3, 2), Composition(arma::zeros<vec>(9)), Composition(arma::zeros<vec>(10))));
    CHECK_THROWS(BoundaryConditions(arma::zeros<mat>(3, 2), Composition(arma::zeros<vec>(10)), Composition(arma::zeros<vec>(11))));
    CHECK_THROWS(BoundaryConditions(arma::zeros<mat>(3, 2), Composition(arma::zeros<vec>(10)), Composition(arma::zeros<vec>(9))));

//    mat a = arma::zeros<mat>(3 + 10, 2);
//    a.col(0).rows(3, 12).fill(0.1);
//    a.col(1).rows(3, 12).fill(0.5);
//    a.col(0).tail(1).fill(1);
//    a.col(1).tail(1).fill(2);
//    a(3, 0) = 3;
//    a(3, 1) = 6;
//    std::cout << a << std::endl;

//    std::cout << a.col(0).rows(3, 12) << std::endl;
//    std::cout << a.rows(0, 2) << std::endl;

//    BoundaryConditions bc(a);

//    std::cout << vec(bc.inletComposition()) << std::endl;
//    std::cout << bc.outletComposition() << std::endl;
//    std::cout << bc(0, 1) << std::endl;
//    std::cout << bc(2, 1) << std::endl;

//    std::cout << bc << std::endl;
}

TEST_CASE("BC from Pipeline")
{
    SUBCASE("getBoundaryConditions")
    {
        Pipeline state(10);
        state.flow() = arma::linspace(1, 10, 10);
        state.pressure() = arma::linspace(10, 100, 10);
        state.temperature() = arma::linspace(100, 1000, 10);

        auto comp = state.composition();
        comp.front() = Composition(arma::linspace(0.1, 1, 10));
        comp.back() = Composition(arma::linspace(1, 0.1, 10));
        state.updateComposition(comp);

        BoundaryConditions bc = state.getBoundaryConditions();

        CHECK(bc.inletFlow() == 1);
        CHECK(bc.outletFlow() == 10);
        CHECK(bc.inletPressure() == 10);
        CHECK(bc.outletPressure() == 100);
        CHECK(bc.inletTemperature() == 100);
        CHECK(bc.outletTemperature() == 1000);

        CHECK(bc.inletComposition()(0) == 0.1);
        CHECK(bc.inletComposition()(9) == 1);
        CHECK(bc.outletComposition()(0) == 1);
        CHECK(bc.outletComposition()(9) == 0.1);
    }
}

TEST_CASE("Constructors")
{
    SUBCASE("empty")
    {}

    SUBCASE("from matrix and compositions")
    {
        CHECK_THROWS(BoundaryConditions(arma::zeros<mat>(3, 1))); // too few columns
        CHECK_THROWS(BoundaryConditions(arma::zeros<mat>(3, 3))); // too many columns
        CHECK_THROWS(BoundaryConditions(arma::zeros<mat>(2, 2))); // too few rows
        CHECK_THROWS(BoundaryConditions(arma::zeros<mat>(4, 2))); // too many rows

        CHECK_NOTHROW(BoundaryConditions bc(arma::zeros<mat>(3, 2)));

        BoundaryConditions bc(
            arma::mat({{1, 2}, {3, 4}, {5, 6}}),
            Composition({1, 2, 3, 4, 5, 6, 7, 8, 9, 10}),
            Composition({2, 3, 4, 5, 6, 7, 8, 9, 10, 11})
        );

        CHECK(bc.inletFlow() == 1);
        CHECK(bc.outletFlow() == 2);

        CHECK(bc.inletPressure() == 3);
        CHECK(bc.outletPressure() == 4);

        CHECK(bc.inletTemperature() == 5);
        CHECK(bc.outletTemperature() == 6);

        CHECK(bc.inletComposition()(0) == 1);
        CHECK(bc.inletComposition()(9) == 10);

        CHECK(bc.outletComposition()(0) == 2);
        CHECK(bc.outletComposition()(9) == 11);

        // defaults
        CHECK(bc.inletFlow().isActive() == true);
        CHECK(bc.outletFlow().isActive() == false);
        CHECK(bc.inletPressure().isActive() == false);
        CHECK(bc.outletPressure().isActive() == true);
        CHECK(bc.inletTemperature().isActive() == true);
        CHECK(bc.outletTemperature().isActive() == false);
    }

    SUBCASE("from doubles and composition")
    {
        BoundaryConditions bc(
            1, 2, 3, 4, 5, 6,
            Composition({1, 2, 3, 4, 5, 6, 7, 8, 9, 10}),
            Composition({2, 3, 4, 5, 6, 7, 8, 9, 10, 11})
        );

        CHECK(bc.inletFlow() == 1);
        CHECK(bc.outletFlow() == 2);

        CHECK(bc.inletPressure() == 3);
        CHECK(bc.outletPressure() == 4);

        CHECK(bc.inletTemperature() == 5);
        CHECK(bc.outletTemperature() == 6);

        CHECK(bc.inletComposition()(0) == 1);
        CHECK(bc.inletComposition()(9) == 10);

        CHECK(bc.outletComposition()(0) == 2);
        CHECK(bc.outletComposition()(9) == 11);

        // defaults
        CHECK(bc.inletFlow().isActive() == true);
        CHECK(bc.outletFlow().isActive() == false);
        CHECK(bc.inletPressure().isActive() == false);
        CHECK(bc.outletPressure().isActive() == true);
        CHECK(bc.inletTemperature().isActive() == true);
        CHECK(bc.outletTemperature().isActive() == false);
    }

    SUBCASE("from SingleCondition and composition")
    {
        BoundaryConditions bc(
            SingleCondition(1, false),
            SingleCondition(2, false),
            SingleCondition(3, false),
            SingleCondition(4, false),
            SingleCondition(5, false),
            SingleCondition(6, false),
            Composition({1, 2, 3, 4, 5, 6, 7, 8, 9, 10}),
            Composition({2, 3, 4, 5, 6, 7, 8, 9, 10, 11})
        );

        CHECK(bc.inletFlow() == 1);
        CHECK(bc.inletFlow().isActive() == false);
        CHECK(bc.outletFlow() == 2);
        CHECK(bc.outletFlow().isActive() == false);

        CHECK(bc.inletPressure() == 3);
        CHECK(bc.inletPressure().isActive() == false);
        CHECK(bc.outletPressure() == 4);
        CHECK(bc.outletPressure().isActive() == false);

        CHECK(bc.inletTemperature() == 5);
        CHECK(bc.inletTemperature().isActive() == false);
        CHECK(bc.outletTemperature() == 6);
        CHECK(bc.outletTemperature().isActive() == false);

        CHECK(bc.inletComposition()(0) == 1);
        CHECK(bc.inletComposition()(9) == 10);

        CHECK(bc.outletComposition()(0) == 2);
        CHECK(bc.outletComposition()(9) == 11);
    }

    SUBCASE("From Pipeline")
    {
        Pipeline state(10);
        state.flow() = arma::linspace(1, 10, 10);
        state.pressure() = arma::linspace(10, 100, 10);
        state.temperature() = arma::linspace(100, 1000, 10);

        auto comp = state.composition();
        comp.front() = Composition(arma::linspace(0.1, 1, 10));
        comp.back() = Composition(arma::linspace(1, 0.1, 10));
        state.updateComposition(comp);

        BoundaryConditions bc(state);

        CHECK(bc.inletFlow() == 1);
        CHECK(bc.outletFlow() == 10);
        CHECK(bc.inletPressure() == 10);
        CHECK(bc.outletPressure() == 100);
        CHECK(bc.inletTemperature() == 100);
        CHECK(bc.outletTemperature() == 1000);

        CHECK(bc.inletComposition()(0) == 0.1);
        CHECK(bc.inletComposition()(9) == 1);
        CHECK(bc.outletComposition()(0) == 1);
        CHECK(bc.outletComposition()(9) == 0.1);
    }
}

TEST_CASE("setBoundarySettings")
{
    BoundaryConditions bc(1, 2, 3, 4, 5, 6);
    bc.setBoundarySettings({"none", "none", "none"});

    // flow
    bc.setBoundarySettings({"none", "none", "none"});
    CHECK(bc.inletFlow().isActive() == false);
    CHECK(bc.outletFlow().isActive() == false);

    bc.setBoundarySettings({"inlet", "none", "none"});
    CHECK(bc.inletFlow().isActive() == true);
    CHECK(bc.outletFlow().isActive() == false);

    bc.setBoundarySettings({"outlet", "none", "none"});
    CHECK(bc.inletFlow().isActive() == false);
    CHECK(bc.outletFlow().isActive() == true);

    bc.setBoundarySettings({"both", "none", "none"});
    CHECK(bc.inletFlow().isActive() == true);
    CHECK(bc.outletFlow().isActive() == true);

    // pressure
    bc.setBoundarySettings({"none", "none", "none"});
    CHECK(bc.inletPressure().isActive() == false);
    CHECK(bc.outletPressure().isActive() == false);

    bc.setBoundarySettings({"none", "inlet", "none"});
    CHECK(bc.inletPressure().isActive() == true);
    CHECK(bc.outletPressure().isActive() == false);

    bc.setBoundarySettings({"none", "outlet", "none"});
    CHECK(bc.inletPressure().isActive() == false);
    CHECK(bc.outletPressure().isActive() == true);

    bc.setBoundarySettings({"none", "both", "none"});
    CHECK(bc.inletPressure().isActive() == true);
    CHECK(bc.outletPressure().isActive() == true);

    // temperature
    bc.setBoundarySettings({"none", "none", "none"});
    CHECK(bc.inletTemperature().isActive() == false);
    CHECK(bc.outletTemperature().isActive() == false);

    bc.setBoundarySettings({"none", "none", "inlet"});
    CHECK(bc.inletTemperature().isActive() == true);
    CHECK(bc.outletTemperature().isActive() == false);

    bc.setBoundarySettings({"none", "none", "outlet"});
    CHECK(bc.inletTemperature().isActive() == false);
    CHECK(bc.outletTemperature().isActive() == true);

    bc.setBoundarySettings({"none", "none", "both"});
    CHECK(bc.inletTemperature().isActive() == true);
    CHECK(bc.outletTemperature().isActive() == true);
}

TEST_CASE("nActiveBoundaryConditions")
{
    BoundaryConditions bc(1, 2, 3, 4, 5, 6);

    std::map<std::string, arma::uword> mapping {
        {"none", 0},
        {"outlet", 1},
        {"inlet", 1},
        {"both", 2}
    };

    for (std::string flow : {"none", "outlet", "inlet", "both"})
    {
        for (std::string pressure: {"none", "outlet", "inlet", "both"})
        {
            for (std::string temperature : {"none", "outlet", "inlet", "both"})
            {
                bc.setBoundarySettings({flow, pressure, temperature});
                arma::uword count = mapping.at(flow) + mapping.at(pressure) + mapping.at(temperature);
                CHECK(bc.nActiveBoundaryConditions() == count);
            }
        }
    }
}

TEST_CASE("inlet/outlet get by index")
{
    const BoundaryConditions bc(1, 4, 2, 5, 3, 6);

    CHECK(bc.inlet(0) == 1);
    CHECK(bc.inlet(0) == bc.inletFlow());
    CHECK(bc.outlet(0) == 4);
    CHECK(bc.outlet(0) == bc.outletFlow());

    CHECK(bc.inlet(1) == 2);
    CHECK(bc.inlet(1) == bc.inletPressure());
    CHECK(bc.outlet(1) == 5);
    CHECK(bc.outlet(1) == bc.outletPressure());

    CHECK(bc.inlet(2) == 3);
    CHECK(bc.inlet(2) == bc.inletTemperature());
    CHECK(bc.outlet(2) == 6);
    CHECK(bc.outlet(2) == bc.outletTemperature());
}

typedef BoundaryConditions::SingleCondition SingleCondition;
TEST_CASE("SingleCondition")
{
    SUBCASE("Constructors")
    {
        SingleCondition s(1, true);
        CHECK(s.isActive() == true);
        CHECK(s.value() == 1);
        s = SingleCondition(2, false);
        CHECK(s.isActive() == false);
        CHECK(s.value() == 2);
    }
}
