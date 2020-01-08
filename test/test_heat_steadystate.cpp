#include "debug.hpp"
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <vector>

#include "heattransfer/heattransferbase.hpp"
#include "heattransfer/steadystate.hpp"
#include "heattransfer/unsteady.hpp"
#include "heattransfer/utils.hpp"
#include "pipeline.hpp"

using namespace std;
using namespace arma;

TEST_SUITE_BEGIN("Steady state heat transfer");

TEST_CASE("Constructor 1")
{
    SUBCASE("Test")
    {
        const double diameter = 0.9;
        const PipeWall pipeWall = PipeWall(
                        {
                            PipeWall::Layer(0.024, Material::steel),
                            PipeWall::Layer(0.007, Material::coating),
                            PipeWall::Layer(0.08, Material::concrete)
                        }
                    );
        const BurialMedium medium = BurialMedium::soil;
        const double burial = -2*diameter;
        const AmbientFluid fluid = AmbientFluid::seawater;
        SteadyStateHeatTransfer heat(diameter, pipeWall, burial, medium, fluid);
        CHECK(heat.size() == 3);
    }

    SUBCASE("Test")
    {
        const double diameter = 0.9;
        const PipeWall pipeWall = PipeWall(
                        {
                            PipeWall::Layer(0.024, Material::steel),
                            PipeWall::Layer(0.007, Material::coating),
                            PipeWall::Layer(0.08, Material::concrete)
                        }
                    );
        const BurialMedium medium = BurialMedium::soil;
        const double burial = 1.2;
        const AmbientFluid fluid = AmbientFluid::seawater;
        SteadyStateHeatTransfer heat(diameter, pipeWall, burial, medium, fluid);
        CHECK(heat.size() == 13); // constant 10 burial layers...

        const double wallThickness = pipeWall.layer(0).thickness() + pipeWall.layer(1).thickness() + pipeWall.layer(2).thickness();
        utils::calcEquivalentBurialLayerWidths(diameter, wallThickness, burial, medium.conductivity());
    }
}

TEST_CASE("thermal resistance calculation")
{
    SUBCASE("1 layer")
    {
        const double dia = 1;
        const double burial = 0; // disable burial
        const PipeWall wall = PipeWall({PipeWall::Layer(0.1, Material::steel)});
        const AmbientFluid fluid = AmbientFluid::seawater;
        SteadyStateHeatTransfer heat(dia, wall, burial, BurialMedium::soil, fluid);

        const double ri = dia/2.0;
        const double ro = dia/2.0 + wall.layer(0).thickness();
        const double lam = wall.layer(0).conductivity();
        const double outerDia = 2*ro;
        const double ho = utils::calcOuterWallFilmCoefficient(outerDia, fluid);
        const double ans = 1/(ri*log(ro/ri)/lam + ri/(ro*ho));
        const double U = heat.getOverallHeatTransferCoefficient();

        CHECK(U == ans);
    }

    SUBCASE("3 layers")
    {
        const double dia = 1;
        const double burial = 0; // disable burial
        const AmbientFluid fluid(1e4, 1, 1, 1, 1);
        const PipeWall wall = PipeWall::defaultPipeWall;
        SteadyStateHeatTransfer heat(dia, wall, burial, BurialMedium::soil, fluid);

        vec r;
        r << dia/2.0
          << dia/2.0 + wall.layer(0).thickness()
          << dia/2.0 + wall.layer(0).thickness() + wall.layer(1).thickness()
          << dia/2.0 + wall.layer(0).thickness() + wall.layer(1).thickness() + wall.layer(2).thickness();

        vec lam = {
            wall.layer(0).conductivity(),
            wall.layer(1).conductivity(),
            wall.layer(2).conductivity()
        };

        const double outerDia = 2*r.tail(1)(0);
        const double ho = utils::calcOuterWallFilmCoefficient(outerDia, fluid);
        const double ans = 1/(r(0)*log(r(1)/r(0))/lam(0) + r(0)*log(r(2)/r(1))/lam(1) + r(0)*log(r(3)/r(2))/lam(2) + r(0)/(r(3)*ho));
        const double U = heat.getOverallHeatTransferCoefficient();

        CHECK(U == ans);
    }
}

TEST_CASE("examples")
{
    const double diameter = 1;
    const double burial = 1;

    const double gasPressure = 1e6;
    const double gasReynoldsNumber = 1e5;
    const double gasHeatCapacity = 1000;
    const double gasViscosity = 1e-5;
    const double gasDensity = 10;
    double ambientTemperature = 300;
    double gasTemperature = 300;

    // no temperature difference
    CHECK(SteadyStateHeatTransfer(diameter, PipeWall::defaultPipeWall, burial, BurialMedium::soil, AmbientFluid::seawater).evaluateInternal(ambientTemperature, gasPressure, gasTemperature, gasReynoldsNumber, gasHeatCapacity, gasViscosity).heatFlux() == Approx(0));

    // no heat transfer at zero reynolds number
    double reyn = 0;
    CHECK(SteadyStateHeatTransfer(diameter, PipeWall::defaultPipeWall, burial, BurialMedium::soil, AmbientFluid::seawater).evaluateInternal(ambientTemperature, gasPressure, gasTemperature, reyn, gasHeatCapacity, gasViscosity).heatFlux() == Approx(0));

    SUBCASE("10 degree diff")
    {
        // 10 degree difference
        const double answer = -16.2854;
        ambientTemperature = 275;
        gasTemperature = 285;
        double zeroburial = 0; // no burial stuff
        auto wall = PipeWall({PipeWall::Layer(1, 1000, -1, -1)}); // 1 meter wall with conductivity of 1000
        auto ambientFluid = AmbientFluid(1e6, AmbientFluid::seawater.viscosity(), Material::seawater); // infinite outer film coefficient (or zero??)
        double q = SteadyStateHeatTransfer(diameter, wall, zeroburial, BurialMedium::soil, ambientFluid).evaluateInternal(ambientTemperature, gasPressure, gasTemperature, gasReynoldsNumber, gasHeatCapacity, gasViscosity).heatFlux();
        q = -4*q/(1.0*gasDensity);
        CHECK(q == Approx(answer));

        // swap around temps
        ambientTemperature = 285;
        gasTemperature = 275;
        q = SteadyStateHeatTransfer(diameter, wall, zeroburial, BurialMedium::soil, ambientFluid).evaluateInternal(ambientTemperature, gasPressure, gasTemperature, gasReynoldsNumber, gasHeatCapacity, gasViscosity).heatFlux();
        q = -4*q/(1.0*gasDensity);
        CHECK(q == Approx(-answer));
    }

    double zeroburial = 0; // no burial stuff
    auto wall = PipeWall({PipeWall::Layer(1, 1000, -1, -1)}); // 1 meter wall with conductivity of 1000
    auto ambientFluid = AmbientFluid(1e6, AmbientFluid::seawater.viscosity(), Material::seawater); // infinite outer film coefficient (or zero??)
    auto heat = SteadyStateHeatTransfer(diameter, wall, zeroburial, BurialMedium::soil, ambientFluid);
    CHECK(heat.getOverallHeatTransferCoefficient() == Approx(1820.48));
}

TEST_SUITE_END();
