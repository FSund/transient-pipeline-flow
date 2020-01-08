#include "heattransfer/utils.hpp"

#include <cstdlib>
#include <cmath>
#include <iostream>

#include "utilities/stringbuilder.hpp"
#include "heattransfer/ambientfluid.hpp"

using arma::vec;
using arma::zeros;
using arma::uword;

double utils::calcGasThermalConductivity(const double pressure)
{
    // correlation from JFH code, not sure of where it comes from
    double lambda = (0.002*pressure/1.0e6 + 0.024); // [W/m K] thermal conductivity of natural gas mixture

    return lambda;
}

double utils::calcOuterWallFilmCoefficient(
        const double diameter,
        const double heatCapacityConstantPressure,
        const double viscosity,
        const double thermalConductivity,
        const double density,
        const double velocity)
{
    // calculates outer film heat transfer coefficient for heating and cooling outside tubes
    const double Pr = heatCapacityConstantPressure*viscosity/thermalConductivity; // Prandtl number
    const double Re = density*velocity*diameter/viscosity; // Reynolds number

    // eq. 7.52 from Bergman, Lavine, Incropera, DeWitt (7th Ed)
    double C;
    double m;
    if (Re >= 4e3 && Re < 4e4)
    {
        C = 0.193;
        m = 0.618;
    }
    else if (Re >= 4e4 && Re < 4e5)
    {
        C = 0.027;
        m = 0.805;
    }
    else if (Re >= 4e5) // this is outside the range defined in the source...
    {
        C = 0.027;
        m = 0.805;
    }
    else
    {
        throw std::invalid_argument(utils::stringbuilder() << "Reynolds number out of range (" << Re << ")");
    }

    const double Nu = C*pow(Re, m)*pow(Pr, 1.0/3.0);

    const double ho = Nu*thermalConductivity/diameter; // outer wall film transfer coefficient

    return ho;
}

double utils::calcOuterWallFilmCoefficient(
        const double diameter,
        const AmbientFluid& fluid)
{
    return calcOuterWallFilmCoefficient(
                diameter, fluid.heatCapacity(),fluid.viscosity(),
                fluid.conductivity(), fluid.density(), fluid.velocity());
}

double utils::calcInnerWallFilmCoefficient(
        const double diameter,
        const double pressure,
        const double ReynoldsNumber,
        const double heatCapacityConstantPressure,
        const double viscosity)
{
    double lambda = utils::calcGasThermalConductivity(pressure); // [W/m K]
    double Nu;
    if (ReynoldsNumber > 1e4)
    {
        // calculates Nusselt number using Dittus-Boelter equation
        // only valid for Re > 10 000 and 0.6 < Pr < 160
        double Pr = heatCapacityConstantPressure*viscosity/lambda; // Prandtl number of gas
        double n = 0.4;
        Nu = 0.023*std::pow(ReynoldsNumber, 0.8)*std::pow(Pr, n); // Nusselt number, using Dittus-Boelter correlation
    }
    else if (ReynoldsNumber <= 1e4 && ReynoldsNumber > 4000)
    {
        // eq. 8.55 in Bergman et al. - "Fundamentals of Heat and Mass Transfer"
        // Note that in using Equation 8.53 or 8.55 to determine h, the thermal conductivity should be
        // evaluated at T_m (the mean or bulk temperature)
        Nu = 3.66;
    }
    else // ReynoldsNumber < 4000
    {
        // disable heat transfer at low Reynolds numbers to avoid numerical instabilities
        Nu = 0;
    }

    double hi = Nu*lambda/diameter; // Inner wall film transfer coefficient [W/(m^2K)]

    return hi;
}

double utils::calcEquivalentBurialLayerWidth(
        const double innerDiameter,
        const double wallThickness,
        const double distanceFromTopOfPipeToSoil,
        const double soilConductivity)
{
    // calculates equivalent soil layer radi and widths, based on equations in OLGA documentation
    // using these widths gives the same results with Unsteady as SteadyState, when thermalized
    double equivalentSoilLayerRadius = calcEquivalentBurialLayerRadius(innerDiameter, wallThickness, distanceFromTopOfPipeToSoil, soilConductivity);
    double equivalentSoilLayerWidth = equivalentSoilLayerRadius -  innerDiameter/2.0 - wallThickness;

    return equivalentSoilLayerWidth;
}

double utils::calcEquivalentBurialLayerRadius(
        const double innerDiameter,
        const double wallThickness,
        const double distanceFromTopOfPipeToSoil,
        const double soilConductivity)
{
    double outerPipeRadius = innerDiameter/2.0 + wallThickness;
    double distanceFromCenterOfPipeToSoil = outerPipeRadius + distanceFromTopOfPipeToSoil;
    double outerPipeDiameter = 2.0*outerPipeRadius;
    double equivalentHeatTransferCoefficient =
            soilConductivity/(outerPipeRadius*std::acosh(2.0*distanceFromCenterOfPipeToSoil/outerPipeDiameter));
    double equivalentSoilLayerRadius = outerPipeRadius*exp(soilConductivity/(equivalentHeatTransferCoefficient*outerPipeRadius));

    return equivalentSoilLayerRadius;
}

vec utils::calcLogSpacedShellWidths(
        const double innerRadius,
        const double outerRadius,
        const uword nShells)
{
    vec logSpacedLayerRadii = arma::logspace(std::log10(innerRadius), std::log10(outerRadius), nShells+1);
    vec logSpacedLayerWidths = zeros<vec>(nShells);
    for (uword i = 0; i < nShells; i++)
    {
        logSpacedLayerWidths(i) = logSpacedLayerRadii(i+1) - logSpacedLayerRadii(i);
    }

    return logSpacedLayerWidths;
}

vec utils::calcEquivalentBurialLayerWidths(
        const double innerDiameter,
        const double wallThickness,
        const double distanceFromTopOfPipeToSoil,
        const double soilConductivity,
        const uword nSoilShells)
{
    double equivalentSoilLayerWidth = utils::calcEquivalentBurialLayerWidth(innerDiameter, wallThickness, distanceFromTopOfPipeToSoil, soilConductivity);
    double outerPipeRadius = innerDiameter/2.0 + wallThickness;

    return calcLogSpacedShellWidths(outerPipeRadius, outerPipeRadius + equivalentSoilLayerWidth, nSoilShells);
}
