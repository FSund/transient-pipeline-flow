#pragma once

#include <armadillo>

#include "heattransfer/radial.hpp"

class Pipeline;
class PipeWall;
class BurialMedium;
class AmbientFluid;

/*!
 * \brief Class that implements steady state heat transfer between gas and
 * pipeline surroundings.
 *
 * This is well documented in Jan Fredrik Helgaker's PhD thesis.
 */
class SteadyStateHeatTransfer : public RadialHeatTransfer
{
public:
    /*!
     * \brief Construct from full description of pipeline.
     * \param diameter Inner diameter [m]
     * \param pipeWall PipeWall instance
     * \param burialDepth Distance from top of pipe to top of burial medium [m]
     * \param burialMedium BurialMedium instace
     * \param ambientFluid AmbientFluid instance
     */
    SteadyStateHeatTransfer(
            const double diameter,
            const PipeWall& pipeWall,
            const double burialDepth,
            const BurialMedium& burialMedium,
            const AmbientFluid& ambientFluid);

    /*!
     * \brief Constructor with default pipe wall, burial medium and ambient medium.
     * \param diameter Inner diameter [m]
     * \param burialDepth Distance from top of pipe to top of burial medium [m]
     */
    SteadyStateHeatTransfer(
            const double diameter = 1.0,
            const double burialDepth = 1.0);

    /*!
     * \brief Evaluate 1d radial steady state heat transfer.
     *
     * Operates on a HeatTransferState and returns a new HeatTransferState,
     * but does not require discretization temperature.
     *
     * \param current Current HeatTransferState
     * \param timeStep Time step [s]
     * \param ambientTemperature Ambient temperature [K]
     * \param gasPressure Gas pressure [Pa]
     * \param gasTemperature Gas temperature [K]
     * \param gasReynoldsNumber Reynolds number of gas [-]
     * \param gasHeatCapacity Gas heat capacity (\f$c_p\f$) [J/(kg K)]
     * \param gasViscosity Gas dynamic viscosity [Pa s] = [kg/m*s]
     * \return HeatTransferState with new heat flux.
     */
    virtual HeatTransferState evaluate(
            const HeatTransferState& current,
            const double timeStep,
            const double ambientTemperature,
            const double gasPressure,
            const double gasTemperature,
            const double gasReynoldsNumber,
            const double gasHeatCapacity,
            const double gasViscosity) const override;

    /*!
     * \brief Internal method used for evaluating steady state heat transfer.
     *
     * This is exposed for testing purposes. We typically use pointers anyway,
     * so this is not accessible without casting to SteadyStateHeatTransfer.
     *
     * \param ambientTemperature Ambient temperature [K]
     * \param gasPressure Gas pressure [Pa]
     * \param gasTemperature Gas temperature [K]
     * \param gasReynoldsNumber Reynolds number of gas [-]
     * \param gasHeatCapacity Gas heat capacity (\f$c_p\f$) [J/(kg K)]
     * \param gasViscosity Gas dynamic viscosity [Pa s] = [kg/m*s]
     * \return HeatTransferState with new heat flux.
     */
    HeatTransferState evaluateInternal(
            const double ambientTemperature,
            const double gasPressure,
            const double gasTemperature,
            const double gasReynoldsNumber,
            const double gasHeatCapacity,
            const double gasViscosity) const;

    /*!
     * \brief Calculate the total heat transfer coefficient U.
     *
     * This is exposed for testing purposes. We typically use pointers anyway,
     * so this is not accessible without casting to SteadyStateHeatTransfer.
     *
     * \param gasPressure Gas pressure [Pa]
     * \param gasReynoldsNumber Gas Reynolds number [-]
     * \param gasHeatCapacityConstantPressure Gas heat capacity (\f$c_p\f$) [J/(kg K)]
     * \param gasViscosity Gas dynamic viscosity [Pa s] = [kg/m*s]
     * \return Total heat transfer coefficient U [W/(m2 K)]
     */
    double calculateHeatTransferCoefficient(
            const double gasPressure,
            const double gasReynoldsNumber,
            const double gasHeatCapacityConstantPressure,
            const double gasViscosity) const;

    //! Get the total heat transfer coefficient of all radial discretization shells.
    //! Does not include the inner film coefficient.
    double getOverallHeatTransferCoefficient() const { return m_overallHeatTransferCoefficient; }
    //! Get the total thermal resistance of all radial discretization shells.
    //! Does not include the inner film coefficient.
    double getOverallThermalResistance() const { return m_overallThermalResistance; }

private:
    //! Total heat transfer coefficient of all radial discretization shells.
    //! Does not include the inner film coefficient.
    double m_overallHeatTransferCoefficient;

    //! Total thermal resistance of all radial discretization shells.
    //! Does not include the inner film coefficient.
    double m_overallThermalResistance;
};
