#include "heattransfer/steadystate.hpp"

#include <cmath>

#include "heattransfer/utils.hpp"
#include "heattransfer/burialmedium.hpp"
#include "heattransfer/pipewall.hpp"
#include "heattransfer/ambientfluid.hpp"

using arma::uword;
using arma::vec;

SteadyStateHeatTransfer::SteadyStateHeatTransfer(
        const double diameter,
        const PipeWall& pipeWall,
        const double burialDepth,
        const BurialMedium& burialMedium,
        const AmbientFluid& ambientFluid):
    RadialHeatTransfer(diameter, pipeWall, burialDepth, burialMedium, ambientFluid)
{
    // calculate thermal resistance and heat transfer coefficient of the
    // pipewall and possible burial layers
    m_overallThermalResistance = 0.0;
    m_overallThermalResistance += arma::sum(log(m_ro/m_ri)/m_conductivity); // thermal resistance of pipe wall

    // outer heat transfer coefficient (outer film coefficient)
    const double outerRadius = m_ro.tail(1)(0);
    const double ho = calculateOuterFilmCoefficient(); // doesn't depend on gas
    m_overallThermalResistance += 1.0/(outerRadius*ho);
    m_overallHeatTransferCoefficient = 1.0/(m_ri(0)*m_overallThermalResistance);
}

SteadyStateHeatTransfer::SteadyStateHeatTransfer(
        const double diameter,
        const double burialDepth):
    SteadyStateHeatTransfer(diameter,
                            PipeWall::defaultPipeWall,
                            burialDepth,
                            BurialMedium::soil,
                            AmbientFluid::seawater)
{}

HeatTransferState SteadyStateHeatTransfer::evaluate(
        const HeatTransferState& /*current*/,
        const double /*timeStep*/,
        const double ambientTemperature,
        const double gasPressure,
        const double gasTemperature,
        const double gasReynoldsNumber,
        const double gasHeatCapacity,
        const double gasViscosity) const
{
    return evaluateInternal(
                ambientTemperature, gasPressure, gasTemperature,
                gasReynoldsNumber, gasHeatCapacity, gasViscosity);
}

HeatTransferState SteadyStateHeatTransfer::evaluateInternal(
        const double ambientTemperature,
        const double gasPressure,
        const double gasTemperature,
        const double gasReynoldsNumber,
        const double gasHeatCapacity,
        const double gasViscosity) const
{
    double U = calculateHeatTransferCoefficient(
                gasPressure, gasReynoldsNumber, gasHeatCapacity, gasViscosity);
    // JFH definition
//    double q = -4.0*U/(m_diameter*gasDensity)*(gasTemperature - ambientTemperature);

    // total heat transfer per meter pipe Q [W/m]
//    double q = U*(M_PI*m_diameter)*(gasTemperature - ambientTemperature);

    // heat flux [W/m2] (energy per square meter)
    const double q = U*(gasTemperature - ambientTemperature);

    return HeatTransferState(q);
}

double SteadyStateHeatTransfer::calculateHeatTransferCoefficient(
        const double gasPressure,
        const double gasReynoldsNumber,
        const double gasHeatCapacityConstantPressure,
        const double gasViscosity) const
{
    const double hi = utils::calcInnerWallFilmCoefficient(
                m_diameter, gasPressure, gasReynoldsNumber,
                gasHeatCapacityConstantPressure, gasViscosity);
    const double inverseU = 1.0/hi + m_ri(0)*m_overallThermalResistance;
    const double U = 1.0/inverseU;

    return U;
}
