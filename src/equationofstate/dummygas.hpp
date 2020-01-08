#pragma once

#include <armadillo>

#include "constants.hpp"
#include "equationofstate/equationofstatebase.hpp"

/*!
 * \brief Dummy EOS implementation used for unit testing. This uses a special
 * molar mass for the gas components, which makes for easier testing.
 */
class DummyGas : public EquationOfStateBase
{
public:
    /*!
     * \brief Construct from Composition.
     * \param composition Composition
     */
    explicit DummyGas(const Composition& composition):
        EquationOfStateBase(composition)
    {
        m_molarMass = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        setComposition(m_composition); // update molar mass of mixture using above
    }

    /*!
     * \brief Override of evaluate, for unit testing.
     * \return A different integer for each property.
     */
    virtual arma::vec evaluate(const double, const double) const override
    {
        arma::vec Z = arma::zeros<arma::vec>(6);
        Z(0) = 1.0; // Z factor
        Z(1) = 2.0; // dZdT_p
        Z(2) = 3.0; // dZdp_T
        Z(3) = 4.0; // dZdT_rho
        Z(4) = 5.0; // c_p
        Z(5) = 6.0; // c_v

        return Z;
    }

    //! Override which just returns 1.0.
    virtual double calculateCompressibility(const double, const double) const override
    {
        return 1.0;
    }
};
