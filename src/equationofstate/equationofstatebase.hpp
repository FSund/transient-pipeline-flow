#pragma once

#include <armadillo>

#include "composition.hpp"

/*!
 * \brief The EquationOfStateBase is an abstract class, the base class for
 * different equations of state.
 *
 * The equation of state is used to determine the compressibility factor Z,
 * partial derivatives of Z (
 * \f$Z\f$, \f$\frac{\partial Z}{\partial T}|_p\f$,
 * \f$\frac{\partial Z}{\partial p}|_T\f$, and
 * \f$\frac{\partial Z}{\partial T}|_\rho\f$
 * ), as well as the heat capacity at constant pressure
 * (\f$c_p\f$) and the heat capacity at constant volume (\f$c_v\f$), from the gas
 * composition, pressure and temperature.
 */
class EquationOfStateBase
{
public:
    //! Have to declare virtual destructor to avoid compiler warnings.
    //! Only declared here, to avoid the inline compiler-generated default destructor.
    virtual ~EquationOfStateBase();

    /*!
     * \brief EquationOfStateBase constructor.
     * \param composition Gas composition (fraction).
     */
    explicit EquationOfStateBase(const arma::vec& composition = Composition::defaultComposition);

    /*!
     * \brief Pure virtual function for evaluation the EOS at constant composition.
     *
     * This function is a stencil for a method for calculating the
     * compressibility factor Z, three partial derivatives, and the heat
     * capacity at constant pressure and constant volume.
     * \param pressure Gas pressure [Pa].
     * \param temperature Gas temperature [K].
     * \return arma::vec containing \f$Z\f$, \f$\frac{\partial Z}{\partial T}|_p\f$,
     * \f$\frac{\partial Z}{\partial p}|_T\f$, \f$\frac{\partial Z}{\partial T}|_\rho\f$,
     * \f$c_p\f$ and \f$c_v\f$
     */
    virtual arma::vec evaluate(const double pressure, const double temperature) const = 0;

    /*!
     * \brief Virtual function for evaluating the EOS at a new composition.
     *
     * This function is a stencil for the method that calculates the
     * compressibility factor Z, three partial derivatives, and the heat
     * capacity at constant pressure and constant volume.
     *
     * This function only calls setComposition() before calling
     * EquationOfStateBase::evaluate(),
     * so it can be overriden if more efficient ways of evaluating the equation
     * of state at a new composition is known.
     *
     * \param pressure Gas pressure [Pa].
     * \param temperature Gas temperature [K].
     * \param composition Gas composition fractions [-].
     * \return arma::vec containing \f$Z\f$, \f$\frac{\partial Z}{\partial T}|_p\f$,
     * \f$\frac{\partial Z}{\partial p}|_T\f$, \f$\frac{\partial Z}{\partial T}|_\rho\f$,
     * \f$c_p\f$ and \f$c_v\f$
     */
    virtual arma::vec evaluate(const double pressure, const double temperature, const arma::vec& composition);

    /*!
     * \brief Method for calculating just the compressibility factor Z at a given pressure and temperature.
     * \param pressure Gas pressure [Pa].
     * \param temperature Gas temperature [K].
     * \return Gas compressibility at the given pressure and temperature.
     */
    virtual double calculateCompressibility(const double pressure, const double temperature) const = 0;

    virtual double calculateStandardDensity() const;

    /*!
     * \brief Set a new composition for the equation of state.
     * \param composition The new composition, fractions in order C1, C2, C3, iC4, nC4, iC5, nC5, C6, N2, CO2.
     * \param force If the composition should be changed even if it's within machine precision of the previous composition.
     * \return True if composition was changed, else false.
     */
    virtual bool setComposition(const arma::vec& composition, const bool force = true);

    double getMolarMassOfMixture() const; //!< Get the molar mass of the gas.
    const arma::vec& getComposition() const; //!< Get the current composition stored in the EOS instance.

protected:
    // components               C1,         C2,         C3,         iC4,        nC4,        iC5,        nC5,        C6          N2,         CO2
    //! The molar mass [g/mol] of the different gas components, in order C1, C2, C3, iC4, nC4, iC5, nC5, C6, N2, CO2.
    arma::vec m_molarMass = {   16.04,      30.07,      44.1,       58.12,      58.12,      72.15,      72.15,      86.18,      28.13,      44.01}; // molar weight of components [g/mol]
    //! The composition of the gas as fractions, in order C1, C2, C3, iC4, nC4, iC5, nC5, C6, N2, CO2.
    arma::vec::fixed<10> m_composition = arma::vec(10);

    //! The molar mass of the gas mixture [g/mol]
    double m_molarMassOfMixture;

    //! Cache density for optimization [kg/m3]
    mutable double m_density = 0;
};
