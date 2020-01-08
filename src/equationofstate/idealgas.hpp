#pragma once

#include <armadillo>

#include "constants.hpp"
#include "equationofstate/equationofstatebase.hpp"

/*!
 * \brief The IdealGas class implements the equation of state for an ideal gas.
 */
class IdealGas : public EquationOfStateBase
{
public:
    /*!
     * \brief IdealGas constructor. Composition is only used to determine molar
     * mass (and density) of the gas.
     */
    explicit IdealGas(const Composition&) {}

    /*!
     * \brief Override. Independent of pressure and temperature.
     *
     * The compressibility is 1, derivatives are zero, and the heat capacities
     * are the heat capacities of an ideal monoatomic gas.
     *
     * \return Z = 1, derivatives = 0 and \f$c_p\f$ and \f$c_v\f$ of an ideal monoatomic gas.
     */
    virtual arma::vec evaluate(const double, const double) const override
    {
        arma::vec Z = arma::zeros<arma::vec>(6);
        Z(0) = 1.0; // Z factor
        Z(1) = 0.0; // derivatives are zero
        Z(2) = 0.0; // derivatives are zero
        Z(3) = 0.0; // derivatives are zero
        Z(4) = 5/2*constants::gasConstant; // c_p of ideal monoatomic gas
        Z(5) = 3/2*constants::gasConstant; // c_v of ideal monoatomic gas
        return Z;
    }

    /*!
     * \brief Override. Independent of pressure and temperature.
     * \return Always returns 1.0 (ideal gas).
     */
    virtual double calculateCompressibility(const double, const double) const override
    {
        return 1.0;
    }
};
