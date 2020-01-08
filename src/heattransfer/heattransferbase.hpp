#pragma once

#include <armadillo>

#include "heattransferstate.hpp"

/*!
 * \brief The HeatTransferBase class is an abstract class, the base class for all
 * heat transfer implementations. It mostly just defines a stencil for the
 * evaluate() function and the HeatTransferState class.
 */
class HeatTransferBase
{
public:
    //! Have to declare virtual destructor to avoid compiler warnings.
    //! Only declared here, to avoid the inline compiler-generated default destructor.
    virtual ~HeatTransferBase();

    // pure virtual method
    // returns q [W/m2] = U*(gasTemperature - ambientTemperature)
    /*!
     * \brief Evaluate heat transfer.
     * \param current Current heat transfer state.
     * \param timeStep Time step [s]
     * \param ambientTemperature Ambient temperature [K]
     * \param gasPressure Gas pressure [Pa]
     * \param gasTemperature Gas temperature [K]
     * \param gasReynoldsNumber Reynolds number of gas [-]
     * \param gasHeatCapacityConstantPressure Heat capacity (\f$c_p\f$) [J/kg K]
     * \param gasViscosity Gas dynamic viscosity [Pa s] = [kg/m*s]
     * \return HeatTransferState instance with heat flux and updated temperatures.
     */
    virtual HeatTransferState evaluate(
            const HeatTransferState& current,
            const double timeStep,
            const double ambientTemperature,
            const double gasPressure,
            const double gasTemperature,
            const double gasReynoldsNumber,
            const double gasHeatCapacityConstantPressure,
            const double gasViscosity) const = 0;

    /*!
     * \brief Make instance of HeatTransferState from heat flux.
     *
     * This is meant to be overriden in more advanced heat transfer
     * implementations that that use the optional temperature property.
     *
     * \param heatFlux Heat flux [W/m2]
     * \return HeatTransferState instance.
     */
    virtual HeatTransferState makeState(const double heatFlux) const;

    /*!
     * \brief Make instance of HeatTransferState from heat flux. Overload.
     *
     * This version just calls HeatTransferBase::makeState(const double), but more
     * advanced implementations should initialize the temperature property from
     * the two temperature arguments.
     *
     * \param heatFlux Heat flux [W/m2]
     * \param ambientTemperature Ambient temperature [K]
     * \param gasTemperature Gas temperature [K]
     * \return HeatTransferState instance with temperature.
     */
    virtual HeatTransferState makeState(
            const double heatFlux,
            const double gasTemperature,
            const double ambientTemperature) const;
};
