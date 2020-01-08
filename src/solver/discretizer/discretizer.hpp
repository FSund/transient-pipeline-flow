#pragma once

#include <armadillo>

class Pipeline;

/*!
 * \brief Discretizer is an abstract class, the base class the implementation of
 * the discretization of the two different sets of governing equations; internal
 * energy and enthalpy.
 *
 * The governing equations are discretized using a implicit backward difference
 * method, as detailed in
 *  - <a href="https://doi.org/10.1016/j.jngse.2018.03.014">Gas composition tracking in transient pipeline flow</a>, Chaczykowski et. al., Journal of Natural Gas Science and Engineering 55 (2018)
 *
 * This finite difference scheme discretizes the governing properties flow,
 * pressure and temperature onto a grid. The derivatives and values for each
 * grid section \f$I\f$ are approximated by different combination of the flow,
 * pressure and temperatures at the grid points, \f$y_i\f$ and \f$y_{i+1}\f$,
 * as follows
 *
 * \f[
 *     \frac{\partial y(I, t_{n+1})}{\partial t} = \frac{y_{i+1}^{n+1} + y_{i}^{n+1} - y_{i+1}^{n} - y_{i}^{n}}{2\Delta t}
 * \f]
 *
 * \f[
 *     \frac{\partial y(I, t_{n+1})}{\partial x} = \frac{y_{i+1}^{n+1} - y_{i}^{n+1}}{\Delta x}
 * \f]
 *
 * \f[
 *     y(I, t_{n+1}) = \frac{y_{i+1}^{n+1} + y_{i}^{n+1}}{2}
 * \f]
 *
 * This sets up a matrix equation \f$Ax = b\f$, where the matrix \f$A\f$
 * contains the coeffcients of each term \f$y_i\f$, the vector \f$x\f$ contain
 * the unknowns \f$\dot m_i\f$ (flow), \f$p_i\f$ (pressure) and \f$T_i\f$
 * (temperature), and the vector \f$b\f$ contain the constant/known terms
 * (from boundary conditions and other knowns).
 *
 * \see MatrixEquation
 */
class Discretizer
{
public:
    //! Have to declare virtual destructor to avoid compiler warnings.
    //! Only declared here, to avoid the inline compiler-generated default destructor.
    virtual ~Discretizer();

    /*!
     * \brief Construct from number of grid points and number of equations and
     * variables.
     *
     * Allocates the matrices Discretizer::m_term_i, Discretizer::m_term_ipp,
     * and the vector Discretizer::m_boundaryTerm.
     *
     * \param nGridPoints Number of grid points
     * \param nEquationsAndVariables Number of equations and variables
     */
    Discretizer(
            const arma::uword nGridPoints,
            const arma::uword nEquationsAndVariables);

    /*!
     * \brief Pure virtual function, stencil for subclasses. This method
     * calculates in the coefficients of \f$y_i\f$ and \f$y_{i+1}\f$ and stores
     * them in Discretizer::m_term_i and Discretizer::m_term_ipp.
     * \param dt Time step
     * \param currentState Current pipeline state
     * \param newState New/guess pipeline state
     */
    virtual void discretize(
            const arma::uword dt,
            const Pipeline& currentState,
            const Pipeline& newState) = 0;

    const arma::cube& term_i() const { return m_term_i; } //!< Get coefficients of \f$y_i\f$
    const arma::cube& term_ipp() const { return m_term_ipp; } //!< Get coefficients of \f$y_{i+1}\f$
    const arma::mat& boundaryTerms() const { return m_boundaryTerm; } //!< Get constant terms

protected:
    // TODO: better description of how m_term_i and m_term_ipp are organized
    /*!
     * \brief The coefficients of \f$y_i\f$ in the discretized governing equations.
     *
     * This is an `arma::cube`, which is organized as follows.
     * `m_term_i(grid point, equation number, variable number)`,
     * where grid points are the usual grid points, equations are the
     * continuity (0), momentum (1) and energy equations (2), and the variables
     * are in order flow (0), pressure (1) and temperature(1).
     * So for example if you want the coefficient at grid point number 5, in the
     * energy equation, for flow, you want
     * `m_term_i(4, 2, 0)` (zero-indexed).
     */
    arma::cube m_term_i;

    /*!
     * \brief The coefficients of \f$y_i\f$ in the discretized governing equations.
     *
     * See Discretizer::m_term_i for a description of the organization of this.
     */
    arma::cube m_term_ipp;

    /*!
     * \brief The constant/known terms in the discretized governing equations.
     *
     * This is a matrix, organized as
     * `m_boundaryTerm(grid point, equation number)`.
     */
    arma::mat m_boundaryTerm;

    double m_gravity = 9.81; //!< Gravity
};
