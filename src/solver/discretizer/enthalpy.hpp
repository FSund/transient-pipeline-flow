#pragma once

#include <armadillo>

#include "solver/discretizer/discretizer.hpp"

class Pipeline;

/*!
 * \brief Implementation of Discretizer for the enthalpy form of the energy
 * equation.
 */
class EnthalpyDiscretizer : public Discretizer
{
public:
    /*!
     * \brief Constructor. Calls Discretizer constructor with nEquations = 3.
     * \param nGridPoints Number of grid points
     */
    explicit EnthalpyDiscretizer(
            const arma::uword nGridPoints);

    /*!
     * \brief Discretize implementation. Calculates the coefficients of the
     * discretized governing equations, with the enthalpy form of the energy
     * equation.
     *
     * This is just a wrapper around discretizeFromPrimitives().
     *
     * \param dt Time step [s]
     * \param currentState Current pipeline state
     * \param newState New/guess pipeline state
     */
    virtual void discretize(
            const arma::uword dt,
            const Pipeline& currentState,
            const Pipeline& newState) override;

    /*!
     * \brief Calculates the coefficients of the discretized governing
     * equations, with the enthalpy form of the energy equation, and stores
     * the result in Discretizer::m_term_i, Discretizer::m_term_ipp and
     * Discretizer::m_boundaryTerm.
     *
     * \param dt Time step [s]
     * \param diameter Inner diameter [m]
     * \param height Height profile [m]
     * \param gridPoints Grid points [m]
     * \param currentSpecificGasConstant Current specific gas constant \f$R_{\rm specific}\f$ [J/(kg K)]
     * \param currentMassFlow Current mass flow [kg/s]
     * \param currentPressure Current pressure [Pa]
     * \param currentTemperature Temperature [K]
     * \param guessMassFlow New mass flow [kg/s]
     * \param guessPressure New pressure [Pa]
     * \param guessTemperature New temperature [K]
     * \param guessFriction New friction factor [-]
     * \param guessHeatCapacityConstantPressure New heat capacity \f$c_p\f$ [J/(kg K)]
     * \param guessHeatFlux New heat flux [W/m2]
     * \param guessDensity New gas density [kg/m3]
     * \param guessCompressibilityFactor New gas compressibility factor \f$Z\f$ [-]
     * \param guess_dZdT_p New partial derivative \f$\frac{\partial Z}{\partial T}|_p\f$ [-]
     * \param guess_dZdp New partial derivative \f$\frac{\partial Z}{\partial p}|_T\f$ [-]
     */
    void discretizeFromPrimitives(
            const arma::uword dt,

            const arma::vec& diameter,
            const arma::vec& height,
            const arma::vec& gridPoints,

            const arma::vec& currentSpecificGasConstant,
            const arma::vec& currentMassFlow,
            const arma::vec& currentPressure,
            const arma::vec& currentTemperature,

            const arma::vec& guessMassFlow,
            const arma::vec& guessPressure,
            const arma::vec& guessTemperature,
            const arma::vec& guessFriction,
            const arma::vec& guessHeatCapacityConstantPressure,
            const arma::vec& guessHeatFlux,
            const arma::vec& guessDensity,
            const arma::vec& guessCompressibilityFactor,
            const arma::vec& guess_dZdT_p,
            const arma::vec& guess_dZdp);
};
