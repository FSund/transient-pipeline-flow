#pragma once

#include <armadillo>

#include "heattransfer/radial.hpp"

/*!
 * \brief Implements 1d radial unsteady heat transfer.
 *
 * This is documented in Jan Fredrik Helgaker's PhD thesis and in
 *  - <a href="10.1016/j.apm.2009.07.017">Sensitivity of pipeline gas flow model to the selection of the equation of state</a>
 *  - <a href="10.1016/j.cherd.2009.06.008">Transient flow in natural gas pipeline - The effect of pipeline thermal model</a>
 */
class UnsteadyHeatTransfer : public RadialHeatTransfer
{
public:
    /*!
     * \brief Construct from complete description of pipeline and surroundings.
     * \param diameter Inner diameter [m]
     * \param pipeWall PipeWall instance
     * \param burialDepth Distance from top of pipe to top of burial medium [m]
     * \param burialMedium BurialMedium instance
     * \param ambientFluid AmbientFluid instanc.
     */
    UnsteadyHeatTransfer(
            const double diameter,
            const PipeWall& pipeWall,
            const double burialDepth,
            const BurialMedium& burialMedium,
            const AmbientFluid& ambientFluid);

    /*!
     * \brief Evaluate 1d radial unsteady heat transfer model.
     *
     * Operates on a HeatTransferState and returns a new HeatTransferState.
     * Requires the discretization temperature HeatTransferState::m_temperature.
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
     * \brief Internal method used for evaluating the unsteady heat transfer model.
     *
     * This is exposed for testing purposes. We typically use pointers anyway,
     * so this is not accessible without casting to UnsteadyHeatTransfer.
     *
     * \param shellTemperature The temperature of each discretization layer [K]
     * \param timeStep Time step [s]
     * \param ambientTemperature Ambient temperature [K]
     * \param gasPressure Gas pressure [Pa]
     * \param gasTemperature Gas temperature [K]
     * \param gasReynoldsNumber Reynolds number of gas [-]
     * \param gasHeatCapacity Gas heat capacity (\f$c_p\f$) [J/(kg K)]
     * \param gasViscosity Gas dynamic viscosity [Pa s] = [kg/m*s]
     * \return HeatTransferState with new heat flux.
     */
    HeatTransferState evaluateInternal(
            const arma::vec& shellTemperature,
            const double timeStep,
            const double ambientTemperature,
            const double gasPressure,
            const double gasTemperature,
            const double gasReynoldsNumber,
            const double gasHeatCapacity,
            const double gasViscosity) const;

    /*!
     * \brief Thermalize the unsteady heat transfer model to steady state.
     *
     * This just performs an infinite time step to find the temperature
     * distribution in the discretization layers at steady state.
     * This could also be done analytically, if we had the analytic solution.
     *
     * \param ambientTemperature Ambient temperature [K]
     * \param gasPressure Gas pressure [Pa]
     * \param gasTemperature Gas temperature [K]
     * \param gasReynoldsNumber Reynolds number of gas [-]
     * \param gasHeatCapacity Gas heat capacity (\f$c_p\f$) [J/(kg K)]
     * \param gasViscosity Gas dynamic viscosity [Pa s] = [kg/m*s]
     * \return HeatTransferState with new heat flux.
     */
    HeatTransferState thermalizeToSteadyState(
            const double ambientTemperature,
            const double gasPressure,
            const double gasTemperature,
            const double gasReynoldsNumber,
            const double gasHeatCapacity,
            const double gasViscosity) const;

private:
    //! Heat transfer coefficient for each discretization layer.
    //! This is the k_i from eq. (2.26) in JFH PhD thesis.
    arma::vec m_heatTransferCoefficient;

    /*!
     * \brief Internal (private) method used for solving the equations in the 1d
     * radial unsteady heat transfer model.
     *
     * \param shellTemperature The temperature of each discretization layer [K]
     * \param timeStep Time step [s]
     * \param ambientTemperature Ambient temperature [K]
     * \param gasPressure Gas pressure [Pa]
     * \param gasTemperature Gas temperature [K]
     * \param gasReynoldsNumber Reynolds number of gas [-]
     * \param gasHeatCapacity Gas heat capacity (\f$c_p\f$) [J/(kg K)]
     * \param gasViscosity Gas dynamic viscosity [Pa s] = [kg/m*s]
     * \return HeatTransferState with new heat flux.
     */
    arma::vec solveEquations(
            const arma::vec& shellTemperature,
            const double timeStep,
            const double gasPressure,
            const double gasTemperature,
            const double ambientTemperature,
            const double gasReynoldsNumber,
            const double gasHeatCapacity,
            const double gasViscosity) const;
};
