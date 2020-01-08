#pragma once

#include <string>
#include <armadillo>

#include "equationofstate/equationofstatebase.hpp"

/*!
 * \brief Implements the Benedict-Webb-Rubin-Starling (BWRS) equation of state.
 *
 * This class implements the
 * <a href="https://en.wikipedia.org/wiki/Benedict%E2%80%93Webb%E2%80%93Rubin_equation#The_BWRS_equation_of_state">Benedict-Webb-Rubin-Starling</a>
 * (BWRS) equation of state. BWRS depends on 11 mixture parameters B_0 and A_0,
 * which are calculated from the pure component parameters \f$A_{0i}\f$ and
 * \f$B_{0i}\f$, as well as the binary interaction coefficients \f$k_{ij}\f$.
 *
 * More information on BWRS can be found in
 * <i>Fluid Properties for Light Petroleum Systems</i>
 * (Kenneth E. Starling, 1973, ISBN 978-0872012936).
 *
 * Critical properties refer to the critical temperature (`m_Tc`), critical
 * density (`m_rhoc`), critical pressure (`m_pc`), and accentric factor (`m_w`),
 * which are used in the evaluation of the BWRS equation. These are found in
 * Table 5 at page 223 of <i>Fluid Properties for Light Petroleum Systems</i>.
 */
class BWRS : public EquationOfStateBase
{
public:
    /*!
     * \brief Constructor that uses a string to select which mixture parameters
     * \f$A_0\f$ and \f$B_0\f$, binary interaction coefficients \f$k_{ij}\f$,
     * and critical properties \f$T_c\f$ etc. to use. The options are "Calsep",
     * "JFH" and "Starling".
     *
     * "Calsep" refers to the parameters tuned for Statpipe, see
     * <i>Tuning af parametre i BWRS-tilstandsligningen til brug i STATPIPES termodynamiske programpakke</i>
     * (Jan Munch, 1985) for more info.
     *
     * "JFH" refer to the parameters found in the original code by Helgaker.
     * These have been shown to give bad results after the solver was rewritten
     * to use molar density (the original form used regular density), so I would
     * recommend not using this.
     *
     * "Staring" refers to the parameters given by Starling in
     * <i>Fluid Properties for Light Petroleum Systems</i>
     * (Kenneth E. Starling, 1973, ISBN 978-0872012936).
     *
     * \param composition Gas composition fraction, in order C1, C2, C3, iC4, nC4, iC5, nC5, C6, N2, CO2.
     * \param parameterSet Which parameter set to use ("Starling", "Calsep" or "JFH").
     * \return An instance of BWRS.
     */
    BWRS(
        const arma::vec& composition = Composition::defaultComposition,
        const std::string& parameterSet = "Calsep"
    );

    /*!
     * \brief Explicit BWRS constructor which allows specifying which files
     * to load mixture parameters \f$A_0\f$ and \f$B_0\f$ and binary interaction
     * coefficients \f$k_{ij}\f$ from, and which critical properties to use.
     *
     * \see loadJFHCriticalProperties()
     * \see loadCalsepCriticalProperties()
     * \see loadStarlingCriticalProperties()
     *
     * \param composition Gas composition fractions, in order C1, C2, C3, iC4, nC4, iC5, nC5, C6, N2, CO2.
     * \param ABparameterFile Path to file which contains the pure component parameters \f$A_{0i}\f$ and \f$B_{0i}\f$
     * \param binaryInteractionTableFile Path to file which contains the interaction parameters \f$k_{ij}\f$
     * \param criticalProperties Which set of critical gas properties should be used.
     */
    static BWRS fromFilePaths(
        const arma::vec& composition = Composition::defaultComposition,
        const std::string& ABparameterFile = std::string(TRANSFLOW_RESOURCE_PATH) + "/equationofstate/bwrs/calsepABparameters.csv",
        const std::string& binaryInteractionTableFile = std::string(TRANSFLOW_RESOURCE_PATH) + "/equationofstate/bwrs/calsepBinaryInteraction.csv",
        const std::string& criticalProperties = "Starling"
    );

    /*!
     * \brief Evaluate the BWRS equation of state at the given pressure and temperature.
     * \param pressure Gas pressure [Pa]
     * \param temperature Gas temperature [K]
     * \return See EquationOfStateBase::evaluate().
     */
    virtual arma::vec evaluate(
            const double pressure,
            const double temperature) const override;

    /*!
     * \brief Calculate the compressibility factor (Z) of the gas at a given pressure and temperature.
     *
     * This function is useful if we just need the compressibility and not all
     * the other properties we can find from the EOS.
     *
     * It contains a call to
     * BWRS::findMolarDensity(),
     * which contains some Newton-Raphson root finding, which are among the
     * heaviest part of this whole class, so I don't expect much time-save over
     * a call to evaluate().
     *
     * \param pressure Gas pressure [Pa].
     * \param temperature Gas temperature [K].
     * \return Compressibility factor Z [-].
     */
    virtual double calculateCompressibility(
            const double pressure,
            const double temperature) const override;

    /*!
     * \brief Set the composition of the EOS.
     *
     * This calls EquationOfStateBase::setComposition(), then updates the non-zero
     * component indices BWRS::m_indices, calculates the new critical pressure
     * and temperature of the gas mixture, before finally updating all the
     * coefficients `m_A0`, `m_B0`, through `m_GAMMA`.
     *
     * \param composition New gas composition.
     * \param force If the composition should be changed even if it's within machine precision of the previous composition.
     * \return True if composition was changed, else false.
     */
    virtual bool setComposition(const arma::vec& composition, const bool force = true) override;

    /*!
     * \brief loadCriticalProperties Load a set of critical properties. Either
     * "Calsep", "JFH" or "Starling".
     * \param name Name of critical properties. Either "Calsep", "JFH" or
     * "Starling".
     */
    void loadCriticalProperties(const std::string name);

    /*!
     * \brief loadParametersAndCriticalProperties Load parameters and critical
     * properties.
     * \param parameterSet Name of parameter set. Either "Calsep", "JFH" or
     * "Starling".
     */
    void loadParametersAndCriticalProperties(const std::string parameterSet);

    /*!
     * \brief Load the critical properties from Starling.
     *
     * This function sets the values of `m_Tc`, `m_rhoc`, `m_w`, `m_pc`,
     * `m_molarMass`, `m_expW`, and `m_R` according to the values found in
     * <i>Fluid Properties for Light Petroleum Systems</i>
     * (Kenneth E. Starling, 1973, ISBN 978-0872012936).
     *
     * The values set by this function have been found to produce the best results.
     */
    void loadStarlingCriticalProperties();

    /*!
     * \brief Load the critical properties from Helgaker.
     *
     * This function sets the values of `m_Tc`, `m_rhoc`, `m_w`, `m_pc`,
     * `m_molarMass`, `m_expW`, and `m_R` according to the values found in the
     * original Helgaker Matlab code.
     *
     * \warning The values set by this function have been found to produce bad results.
     * This is likely due to the rewrite to a molar density based solver, as
     * the equations are originally given by Starling.
     * Use loadStarlingCriticalProperties() instead.
     */
    void loadJFHCriticalProperties();

    /*!
     * \brief Load the critical properties from the Calsep report.
     *
     * This function sets the values of `m_Tc`, `m_rhoc`, `m_w`, `m_pc`,
     * `m_molarMass`, `m_expW`, and `m_R` according to the values found in the
     * original Helgaker Matlab code.
     *
     * \warning The values set by this function have been found to produce bad results.
     * This is likely due to the rewrite to a molar density based solver, as
     * the equations are originally given by Starling.
     * Use loadStarlingCriticalProperties() instead.
     */
    void loadCalsepCriticalProperties();

    /*! Load the mixture parameters \f$A_0\f$ and \f$B_0\f$, and binary
     * interaction coefficients \f$k_{ij}\f$ from Gassco.
     */
    void loadGasscoParameters();

    /*! Load the mixture parameters \f$A_0\f$ and \f$B_0\f$, and binary
     * interaction coefficients \f$k_{ij}\f$ from the Calsep report.
     */
    void loadCalsepParameters();

    /*! Load the mixture parameters \f$A_0\f$ and \f$B_0\f$, and binary
     * interaction coefficients \f$k_{ij}\f$ from the Starling BWRS book.
     */
    void loadStarlingParameters();

    /*!
     * \brief Load mixture parameters \f$A_0\f$ and \f$B_0\f$, and binary
     * interaction coefficients \f$k_{ij}\f$ coefficients from specific files.
     * \param ABparameterFile Path to AB parameter file.
     * \param binaryInteractionTableFile Path to binary interaction coefficients file.
     */
    void loadParameterFiles(
            const std::string& ABparameterFile,
            const std::string& binaryInteractionTableFile);

    //! Enable constant heat capacity \f$c_p\f$ and \f$c_v\f$
    void enableConstantHeatCapacities();

    /*!
     * \brief Find the gas density at a given pressure and temperature.
     *
     * This is just a wrapper around findMolarDensity() and a conversion from
     * molar density [mol/m3] to density [kg/m3]
     *
     * \param pressure Gas pressure [Pa].
     * \param temperature Gas temperature [K].
     * \param tolerance Newton-Raphson tolerance.
     * \return Gas density [kg/m3].
     */
    double findDensity(
            const double pressure,
            const double temperature,
            const double tolerance = 1e-4) const;

//    double getDerivativeOfZwrtDensityAtConstantTemperature(const double density, const double temperature) const;

    double getGasConstant() const { return m_R; } //!< Get the gas constant R.

    //! Get the critical pressure of the gas mixture [Pa]
    double getMixtureCriticalPressure() const { return m_criticalPressureOfMixture; }

    //! Get the critical temperature of the gas mixture [K]
    double getMixtureCriticalTemperature() const { return m_criticalTemperatureOfMixture; }

protected:

    // parameters from Starling book
    //                      C1          C2          C3          iC4         nC4         iC5         nC5         C6          N2          CO2

    //! Critical temperature [K]
    arma::vec m_Tc        = arma::vec({ 190.69,     305.39,     369.89,     408.13,     425.19,     460.37,     469.49,     507.29,     126.15,     304.15});           // critical temperature
    //! Critical density [kg/m3]
    arma::vec m_rhoc      = arma::vec({ 1.00500e+4, 6.75659e+3, 4.99936e+3, 3.80118e+3, 3.92132e+3, 3.24694e+3, 3.21491e+3, 2.71673e+3, 1.10992e+4, 1.06379e+4});       // critical molar density [mol/m3]
    //! Accentric factor [-]
    arma::vec m_w         = arma::vec({ 0.013,      0.1018,     0.157,      0.183,      0.197,      0.226,      0.252,      0.302,      0.035,      0.21});             // accentric factor
    //! Critical pressure [Pa]
    arma::vec m_pc        = arma::vec({ 45.96,      48.839,     42.5,       36.48,      37.96,      33.81,      33.69,      27.34,      33.99,      73.825})*1e5;       // critical pressure
    //! Molar mass of the different gas components [g/mol]
    arma::vec m_molarMass = arma::vec({ 16.042,     30.068,     44.094,     58.12,      58.12,      72.146,     72.146,     86.172,     28.016,     44.01});            // [g/mol]
    //! Exponential of the accentric factor.
    arma::vec m_expW = exp(-3.8*m_w);

    /*!
     * \brief The gas constant [J/(K mol)]
     *
     * The exact value of the gas constant is approx 8.3145 J/Kmol,
     * but if we are using the critical properties and other parameters from
     * Starling, these have been determined using an older definition of the
     * gas constant, which is 8.3160 J/Kmol.
     */
    double m_R = 8.3160; // gas constant from Starling [m3 Pa / K mol]

    /*!
     * \brief Binary interaction coefficients \f$k_{ij}\f$.
     *
     * The binary interaction coefficients are stored in matrix of size 10x10,
     * where the binary interactions between component i and j are stored at
     * location \f$(i, j)\f$ in the matrix. The matrix is symmetric, so
     * \f$k(i, j) == k(j, i)\f$. The components are in the usual order
     * (C1, C2, C3, iC4, nC4, iC5, nC5, C6, N2, CO2).
     *
     * These can be found in Table 1 at page 227 of <i>Fluid Properties for Light Petroleum Systems</i>.
     */
    arma::mat m_binaryInteractionParameterTable;

    /*!
     * \brief Pure component parameters Bi.
     *
     * From Table 1 at page 221 of <i>Fluid Properties for Light Petroleum Systems</i>.
     */
    arma::vec m_Bi;

    /*!
     * \brief Pure component parameters Ai.
     *
     * From Table 1 at page 221 of <i>Fluid Properties for Light Petroleum Systems</i>.
     */
    arma::vec m_Ai;

    double m_A0; //!< A coefficient used when evaluating the BWRS-equation. Independent of pressure and temperature.
    double m_B0; //!< A coefficient used when evaluating the BWRS-equation. Independent of pressure and temperature.
    double m_C0; //!< A coefficient used when evaluating the BWRS-equation. Independent of pressure and temperature.
    double m_D0; //!< A coefficient used when evaluating the BWRS-equation. Independent of pressure and temperature.
    double m_E0; //!< A coefficient used when evaluating the BWRS-equation. Independent of pressure and temperature.
    double m_a; //!< A coefficient used when evaluating the BWRS-equation. Independent of pressure and temperature.
    double m_b; //!< A coefficient used when evaluating the BWRS-equation. Independent of pressure and temperature.
    double m_c; //!< A coefficient used when evaluating the BWRS-equation. Independent of pressure and temperature.
    double m_d; //!< A coefficient used when evaluating the BWRS-equation. Independent of pressure and temperature.
    double m_ALPHA; //!< A coefficient used when evaluating the BWRS-equation. Independent of pressure and temperature.
    double m_GAMMA; //!< A coefficient used when evaluating the BWRS-equation. Independent of pressure and temperature.

    double m_criticalPressureOfMixture; //!< Critical pressure of the gas mixture.
    double m_criticalTemperatureOfMixture; //!< Critical temperature of the gas mixture.

    bool m_useConstantHeatCapacities = false; //!< Flag to set if we want to use constant heat capacity \f$c_p\f$ and \f$c_v\f$.

    /*!
     * \brief Calculates all the coefficients used for evaluating BWRS that are independent of pressure and temperature.
     *
     * A lot of different coefficients are used when evaluating the BWRS
     * equation to find the compressibility and the partial derivatives of Z.
     * Many of these are independent of pressure and temperature, so they can be
     * tabulated for the current gas composition.
     *
     * This optimization is most effective if we use constant composition, but
     * should also improve the evaluation of the formulas otherwise.
     */
    void calculateCoefficients();

    /*!
     * \brief Set up BWRS::m_indices to reflect which gas fractions are non-zero.
     *
     * The BWRS equation has some issues if it is evaluated with gas fractions
     * of zero, so to avoid this we store the indices of the non-zero fractions
     * and only loop over those when evaluating the equation.
     *
     * Call this function to find the non-zero indices and store them in BWRS::m_indices.
     */
    void findNonZeroComponents();

    arma::uvec m_indices; //!< The indices of the non-zero gas fractions (components).

    /*!
     * \brief Find the molar density of the gas at a given pressure and temperature.
     *
     * This is most likely the most computationally heavy function in this
     * class, since it contains a Newton-Raphson root finding, meaning that the
     * BWRS equation will be evaluated multiple times until it reaches
     * convergence.
     *
     * The convergence criterion is implemented as
     *
     *     if (std::abs(change/previous) < tolerance) break;
     *
     * \param pressure Gas pressure [Pa].
     * \param temperature Gas temperature [K].
     * \param tolerance Convergence criterion for the Newton-Raphson method.
     * \return Molar density [mol/m3]
     */
    double findMolarDensity(
            const double pressure,
            const double temperature,
            const double tolerance = 1e-4) const;
};
