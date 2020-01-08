#pragma once

#include <armadillo>
#include <string>

#include "equationofstate/equationofstatebase.hpp"

/*!
 * \brief The GERG04 class implements the GERG 2004 equation of state.
 *
 * This is a higly accurate, but extremely complex equation of state. It has
 * explicit forms for a lot of properties like heat capacity etc., which are
 * utilized where possible.
 *
 * \warning Note that the composition stored internally in GERG04::X is not in
 * the typical order, but in order C1, C2, C3, iC4, nC4, iC5, nC5, C6, N2, CO2.
 *
 * More info:
 *  - <a href="http://www.gerg.eu/public/uploads/files/publications/technical_monographs/tm15_04.pdf"><i>The GERG-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures</i> (O. Kunz and W. Wagner, J. Chem. Eng. Data. 2012, 57, 11, 3032-3091)</a>
 *  - <a href="https://doi.org/10.1021/je300655b"><i>The GERG-2008 Wide-Range Equation of State for Natural Gases and Other Mixtures: An Expansion of GERG-2004</i> (O. Kunz and W. Wagner, J. Chem. Eng. Data. 2012, 57, 11, 3032-3091)</a>
 */
class GERG04 : public EquationOfStateBase
{
public:
    /*!
     * \brief GERG04 constructor.
     *
     * Initializes the coefficients that only depend on composition.
     *
     * \param composition Gas composition in arma::vec, in the order C1, C2, C3, iC4, nC4, iC5, nC5, C6, N2, CO2.
     */
    explicit GERG04(const arma::vec& composition = Composition::defaultComposition);

    /*!
     * \brief Evaluate the GERG 2004 equation of state at the given pressure and temperature.
     *
     * This is a wrapper around evaluateAllProperties(), which selects the
     * wanted properties from that function. This could perhaps be optimized a
     * bit, since some of the calculations are not required. But the brunt of
     * the cpu time is spent in findDensity() anyway, so there probably isn't
     * that much to save.
     *
     * \param pressure Gas pressure [Pa]
     * \param temperature Gas temperature [K]
     * \return See EquationOfStateBase::evaluate().
     */
    virtual arma::vec evaluate(const double pressure, const double temperature) const override;

    /*!
     * \brief Evaluate all available gas properties at a given pressure and temperature.
     *
     * We have implemented 15 of the explicit equations for different properties
     * given in the GERG 2004 documentation.
     *
     * The method returns an arma::vec containing:
     *  - compressibility factor <i>Z</i>
     *  - partial derivative of Z wrt. temperature at constant pressure \f$\frac{\partial Z}{\partial T}|_p\f$
     *  - partial derivative of Z wrt. pressure at constant temperature \f$\frac{\partial Z}{\partial p}|_T\f$
     *  - partial derivative of Z wrt. temperature at constant density \f$\frac{\partial Z}{\partial T}|_\rho\f$
     *  - entropy <i>S</i>
     *  - internal energy <i>U</i>
     *  - heat capacity at constant volume \f$c_v\f$
     *  - enthalpy <i>H</i>
     *  - heat capacity at constant pressure \f$c_p\f$
     *  - gibbs free energy <i>G</i>
     *  - Joule-Thomson coefficient
     *  - speed of sound
     *  - isothermal throttling coefficient
     *  - density <i>ρ</i>
     *  - isentropic exponent \f$\gamma = \frac{c_p}{c_v}\f$
     *
     * \warning Other than the compressibility factor, the partial derivatives,
     * and the heat capacity, the output of this function has not been tested
     * properly. Feel free to implement tests for the other results.
     *
     * \param pressure Gas pressure [Pa].
     * \param temperature Gas temperature [K].
     * \return An arma::vec containing all the properties. See the detailed description for more details.
     */
    arma::vec evaluateAllProperties(const double pressure, const double temperature) const;

    /*!
     * \brief Calculate compressibility Z at given pressure and temperature.
     * \param pressure Gas pressure [Pa].
     * \param temperature Gas temperature [K].
     * \return Compressibility factor Z [-].
     */
    virtual double calculateCompressibility(const double pressure, const double temperature) const override;

    /*!
     * \brief Set the composition.
     *
     * This calls EquationOfStateBase::setComposition(), updates the non-zero
     * component indices GERG04::m_indices, and calculates the inverse reducing function
     * for mixture density GERG04::rhored and the reducing function for mixture
     * temperature GERG04::tred.
     *
     * \param composition New gas composition.
     * \param force If the composition should be changed even if it's within machine precision of the previous composition.
     * \return True if composition was changed, else false.
     */
    virtual bool setComposition(const arma::vec& composition, const bool force = true) override;

    /*!
     * \brief Find the density of the gas at a given pressure and temperature.
     * \param pressure Gas pressure [Pa].
     * \param temperature Gas temperature [K].
     * \return Gas density [kg/m3]
     */
    double findDensity(const double pressure, const double temperature) const;

    /*!
     * \brief Find the soundspeed in the gas at a given temperature and density.
     * \param temperature Gas temperature [K].
     * \param density Gas density [kg/m3]
     * \return Speed of sound in gas [m/s].
     */
    double findSpeedOfSound(const double temperature, const double density) const;

//    void useBadHeatCapacities(bool useBadHeatCapacities);

//    arma::vec getOutput() const;

    /*!
     * \brief Get the indices of the non-zero components in the composition.
     *
     * \warning Note that the components are stored internally in GERG04 in a
     * different order than the input order and the usual order, so it's not
     * straight-forward to relate this to usual composition.
     *
     * \return The indices of the non-zero components.
     */
    const arma::uvec& indicesOfNonZeroComponents() const { return m_indices; }

private:
    /*!
     * \brief Internal (private) function used in the process of evaluating the GERG 2004 equations.
     *
     * Uses Newton-Raphson to find the density, so will be pretty cpu heavy
     * (depending on the starting point).
     *
     * \param pressure Gas pressure [P].
     * \param temperature Gas Temperature [K].
     * \param density Gas density [kg/m3].
     * \param aroidelta Output argument - left sum in eq. (7.21b) in TM15 (also appears in other equations).
     * \param arijdelta Output argument - the right (double) sum in eq. (7.21b) in TM15 (also appears in other equations).
     * \param aroideltadelta Output argument - the left sum in eq. (7.21c) in TM15 (also appears in other equations).
     * \param arijdeltadelta Output argument - the right (double) sum in eq. (7.21c) in TM15 (also appears in other equations).
     * \return The gas density [kg/m3] as well as many other factors via the output arguments.
     */
    double findDensity(
            const double pressure,
            const double temperature,
            double& density, // output
            double& aroidelta, // output
            double& arijdelta, // output
            double& aroideltadelta, // output
            double& arijdeltadelta // output
            ) const;

    /*!
     * \brief Internal (private) function used in the process of evaluating the GERG 2004 equations.
     * \param tred_temperature Inverse reduced mixture temperature.
     * \param start_rhored Reduced mixture density.
     * \param start_rhored_pow_minusOne pow(reduced mixture density, -1)
     * \param start_rhored_pow_minusTwo pow(reduced mixture density, -2)
     * \param aroidelta Output argument -- left sum in eq. (7.21b) in TM15 (also appears in other equations).
     * \param arijdelta Output argument -- the right (double) sum in eq. (7.21b) in TM15 (also appears in other equations).
     * \param aroideltadelta Output argument -- the left sum in eq. (7.21c) in TM15 (also appears in other equations).
     * \param arijdeltadelta Output argument -- the right (double) sum in eq. (7.21c) in TM15 (also appears in other equations).
     */
    void evaluateAlpha_roi_deltas(
            const double tred_temperature,
            const double start_rhored,
            const double start_rhored_pow_minusOne,
            const double start_rhored_pow_minusTwo,
            double& aroidelta, // output
            double& arijdelta, // output
            double& aroideltadelta, // output
            double& arijdeltadelta // output
            ) const;

    /*!
     * \brief Internal (private) function that updates the non-zero components.
     *
     * This function updates GERG04::m_indices, GERG04::m_firstIndices, and
     * GERG04::m_lastIndices according to the current composition stored in
     * GERG04::X.
     */
    void setNonZeroComponents();

    /*!
     * \brief Set the pressure and temperature independent coefficients.
     *
     * This calculates and sets the pressure and temperature independent factors
     * GERG04::rhored and GERG04::tred, which are stored for efficient'
     * computation.
     */
    void calculateCoefficients();

    const arma::uword N = 10; //!< Number of components
    arma::vec X; //!< Composition in order CH4, N2, CO2, C2H6, C3H8, nC4H10, iC4H10, nC5H12, iC5H12, nC6H14

    //! Critical temperature [K] for components CH4, N2, CO2, C2H6, C3H8, nC4H10, iC4H10, nC5H12, iC5H12, nC6H14
    arma::vec Tc =   {190.56,         126.19,         304.13,         305.32,         369.83,         425.13,         407.82,         469.7,          460.35,         507.82};
    //! Critical density [kg/m3] for components CH4, N2, CO2, C2H6, C3H8, nC4H10, iC4H10, nC5H12, iC5H12, nC6H14
    arma::vec rhoc = {10.14*16.04,    11.18*28.01,    10.62*44.01,    6.87*30.07,     5*44.1,         3.92*58.12,     3.86*58.12,     3.21*72.15,     3.27*72.15,     2.71*86.18};
    //                ch4,            n2,             c02,            c2h6,           c3h8,           nc4h10,         ic4h10,         nc5h12,         ic5h12          nc6h14

    /*!
     * \brief The indices of the non-zero gas fractions (components).
     *
     * \warning Note that the components are stored internally in GERG04 in a
     * different order than the input order and the usual order, so it's not
     * straight-forward to relate this to usual composition.
     */
    arma::uvec m_indices;

    /*!
     * \brief The indices of the two first non-zero gas fractions (components).
     *
     * \warning Note that the components are stored internally in GERG04 in a
     * different order than the input order and the usual order, so it's not
     * straight-forward to relate this to usual composition.
     */
    arma::uvec m_firstIndices;

    /*!
     * \brief The indices of the two last non-zero gas fractions (components).
     *
     * \warning Note that the components are stored internally in GERG04 in a
     * different order than the input order and the usual order, so it's not
     * straight-forward to relate this to usual composition.
     */
    arma::uvec m_lastIndices;

    static arma::mat betav; //!< \f$\beta_{v, ij}\f$ density interaction coefficient, from table A3.8
    static arma::mat betat; //!< \f$\beta_{T, ij}\f$ temperature coefficient, from table A3.8
    static arma::mat gammav; //!< \f$\gamma_{v, ij}\f$ density interaction coefficient, from table A3.8
    static arma::mat gammat; //!< \f$\gamma_{T, ij}\f$ temperature interaction coefficient, from table A3.8

    // Kpol coefficients - used in the first sum in alpha^r_{oi} and derivatives
    static arma::mat noipol; //!< \f$n_{oi, k}\f$, from Table A3.2 in TM15, used in the first sum of alpha^r_{oi}
    static arma::mat doipol; //!< \f$d_{oi, k}\f$, from Table A3.2 in TM15, used in the first sum of alpha^r_{oi}
    static arma::mat toipol; //!< \f$t_{oi, k}\f$, from Table A3.2 in TM15, used in the first sum of alpha^r_{oi}

    // Kexp coefficient - used in the second sum in alpha^r_{oi} and derivatives
    static arma::mat noiexp; //!< \f$n_{oi, k}\f$, from Table A3.3 in TM15, used in the second sum of alpha^r_{oi}
    static arma::mat doiexp; //!< \f$d_{oi, k}\f$, from Table A3.3 in TM15, used in the second sum of alpha^r_{oi}
    static arma::mat coiexp; //!< \f$c_{oi, k}\f$, from Table A3.3 in TM15, used in the second sum of alpha^r_{oi}
    static arma::mat toiexp; //!< \f$t_{oi, k}\f$, from Table A3.3 in TM15, used in the second sum of alpha^r_{oi}

    static arma::mat Fij; //!< \f$F_{ij}\f$ interaction coefficient

    static arma::cube nijpol; //!< \f$n_{ij, k}\f$, from Table A3.7 in TM15, used in the first sum of alpha^r_{ij}
    static arma::cube dijpol; //!< \f$d_{ij, k}\f$, from Table A3.7 in TM15, used in the first sum of alpha^r_{ij}
    static arma::cube tijpol; //!< \f$t_{ij, k}\f$, from Table A3.7 in TM15, used in the first sum of alpha^r_{ij}

    static arma::cube nijexp; //!< \f$n_{ij, k}\f$, from Table A3.7 in TM15, used in the second sum of alpha^r_{ij}
    static arma::cube dijexp; //!< \f$d_{ij, k}\f$, from Table A3.7 in TM15, used in the second sum of alpha^r_{ij}
    static arma::cube tijexp; //!< \f$t_{ij, k}\f$, from Table A3.7 in TM15, used in the second sum of alpha^r_{ij}

    static arma::cube nuijexp; //!< \f$\nu_{ij, k}\f$, from Table A3.7 in TM15, used in the second sum of alpha^r_{ij}
    static arma::cube epijexp; //!< \f$\epsilon_{ij, k}\f$, from Table A3.7 in TM15, used in the second sum of alpha^r_{ij}
    static arma::cube beijexp; //!< \f$\beta_{ij, k}\f$, from Table A3.7 in TM15, used in the second sum of alpha^r_{ij}
    static arma::cube gaijexp; //!< \f$\gamma_{ij, k}\f$, from Table A3.7 in TM15, used in the second sum of alpha^r_{ij}

    // other coefficients
    static arma::mat noik; //!< \f$n^o_{oi, k}\f$, from Table A3.1 in TM15, used in alpha^o_{oi}
    static arma::mat voik; //!< \f$n^o_{oi, k}\f$, from Table A3.1 in TM15, used in alpha^o_{oi}

    // constant terms
    static arma::cube nijpol_times_tijpol; //!< \f$n_{ij, k} \times t_{ij, k}\f$ (GERG04::nijpol x GERG04::tijpol)
    static arma::cube nijpol_times_tijpol_times_tijpol_minus_one; //!<  \f$(n_{ij, k} \times t_{ij, k}) \times (t_{ij, k} - 1)\f$ ((GERG04::nijpol x GERG04::tijpol) x (GERG04::tijpol - 1.0))
    static arma::cube nijexp_times_tijexp; //!< \f$n_{ij, k} \times t_{ij, k}\f$ (GERG04::nijexp x GERG04::tijexp)
    static arma::cube nijexp_times_tijexp_times_tijexp_minus_one; //!< \f$(n_{ij, k} \times t_{ij, k}) \times (t_{ij, k} - 1)\f$ ((GERG04::nijpol x GERG04::tijpol) x (GERG04::tijpol - 1.0))
    static arma::cube nijpol_times_dijpol; //!< \f$n_{ij, k} \times d_{ij, k}\f$ (GERG04::nijpol x GERG04::dijpol)
    static arma::cube nijpol_times_dijpol_times_tijpol; //!< \f$n_{ij, k} \times d_{ij, k} \times t_{ij, k}\f$ (GERG04::nijpol x GERG04::dijpol x GERG04::tijpol)
    static arma::cube nijpol_times_dijpol_times_dijpol_minus_one; //!< \f$n_{ij, k} \times d_{ij, k} \times(d_{ij, k} - 1) \f$ (GERG04::nijpol x GERG04::dijpol x (GERG04::dijpol - 1.0))

    double rhored; //!< Inverse reducing function for mixture density \f$1/\rho_r(\bar x)\f$.
    double tred; //!< Reducing function for mixture temperature \f$T_r(\bar x)\f$.

    double Ra; //!< Gas constant of mixture [J/mol K]
};
