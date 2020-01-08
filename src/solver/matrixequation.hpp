#pragma once

#include <armadillo>

class BoundaryConditions;

/*!
 * \brief The MatrixEquation class sets up the matrix equation from the system
 * of equations found from the governing equations.
 *
 * The governing equations are discretized using a implicit backward difference
 * method, as detailed in
 *  - <a href="https://doi.org/10.1016/j.jngse.2018.03.014">Gas composition tracking in transient pipeline flow</a>, Chaczykowski et. al., Journal of Natural Gas Science and Engineering 55 (2018)
 *
 * This class supports all combinations of inlet and outlet boundary conditions,
 * so we can for example use flow at the inlet, pressure at the outlet, and
 * temperature at the inlet, or any other combination. Also supports
 * over-determined systems.
 *
 * This class is designed for use together with a Discretizer instance,
 * which sets up the individual elements of the matrix. MatrixEquation just
 * fills in the matrix using results from Discretizer, and solves the matrix
 * equation.
 *
 * The matrix equation is solved using the sparse matrix solver arma::spsolve,
 * which uses superlu internally. This solver only works for critically
 * determined systems, not for over-determined systems. In that case the regular
 * solver is used, and the sparse matrix has to be converted to a regular
 * matrix.
 *
 * Most of the elements in the matrix are zero, so we use sparse matrices. The
 * non-zero elements are all located near the diagonal.
 *
 * \see Discretizer
 */
class MatrixEquation
{
public:
    //! Declared to avoid the inline compiler-generated default destructor.
    ~MatrixEquation();

    /*!
     * \brief Fill in the coefficient matrix A (MatrixEquation::m_coefficients)
     * and constants vector b (MatrixEquation::m_constants) in the matrix
     * equation Ax = b.
     *
     * Here `term_i` and `term_ipp` refer to the coefficients of the variables
     * in the discretized governing equations (flow, pressure and temperature --
     * \f$y_i\f$ in general). See for example eq. (14)-(16) in
     * <a href="https://doi.org/10.1016/j.jngse.2018.03.014">Gas composition tracking in transient pipeline flow</a>.
     * These coefficients are calculated by the Discretizer class.
     *
     * The boundary terms are the constant terms, the elements of the vector
     * \f$b\f$ in the matrix equation \f$Ax = b\f$.
     *
     * MatrixEquation supports all combinations of inlet and outlet boundary
     * conditions, but giving more than 3 boundary conditions lead to an
     * over-determined system, and giving less than 3 leads to a
     * under-determined system.
     *
     * \see Discretizer::discretize()
     *
     * \param nGridPoints Number of grid points
     * \param nEquationsAndVariables Number of equations and variables
     * \param boundaryConditions The boundary conditions
     * \param term_i Matrix coefficients at point i (see above)
     * \param term_ipp Matrix coefficients at point i+1 (see above)
     * \param boundaryTerms The boundary terms (see above)
     */
    void fillCoefficientMatrixAndConstantsVector(
            const arma::uword nGridPoints,
            const arma::uword nEquationsAndVariables,
            const BoundaryConditions& boundaryConditions,
            const arma::cube& term_i,
            const arma::cube& term_ipp,
            const arma::mat& boundaryTerms); // fills self matrix and vector

    /*!
     * \brief Solve the matrix equation Ax = b.
     *
     * This solves the matrix equation, after the coefficient matrix
     * MatrixEquation::m_coefficients and constant vector
     * MatrixEquation::m_coefficients has been filled in by
     * fillCoefficientMatrixAndConstantsVector().
     *
     * \param nGridPoints Number of grid points
     * \param nEquationsAndVariables Number of equations and variables
     * \param boundaryConditions The boundary conditions
     * \return The solution of the matrix equation (x), reshaped to contain
     * flow, pressure and temperature in separate columns.
     */
    arma::mat solve(
            const arma::uword nGridPoints,
            const arma::uword nEquationsAndVariables,
            const BoundaryConditions& boundaryConditions) const;

    //! Get coefficient matrix A. For testing purposes.
    const arma::sp_mat& coefficients() const { return m_coefficients; }
    //! Get constants vector b. For testing purposes.
    const arma::vec& constants() const { return m_constants; }

private:
    arma::sp_mat m_coefficients; //!< Coefficient matrix A
    arma::vec m_constants; //!< Constants vector b

    /*!
     * \brief Reshape output from solving the matrix equation into a matrix
     * containing flow, pressure and temperature as columns.
     *
     * Reshapes the vector x, the solution of the matrix equation Ax = b,
     * into a matrix with three columns. Also inserts the boundary conditions
     * where at the correct locations.
     *
     * Supports all combinations of inlet/outlet boundary settings.
     *
     * \param x Solution of matrix equation
     * \param boundaryConditions Boundary conditions
     * \param nGridPoints Number of grid points
     * \param nVariables Number of variables
     * \return Matrix containing flow, pressure and temperature columns.
     */
    arma::mat reshapeSolverOutput(
            const arma::vec& x,
            const BoundaryConditions& boundaryConditions,
            const arma::uword nGridPoints,
            const arma::uword nVariables) const;

    /*!
     * \brief Internal method that solves the matrix equation.
     *
     * \see solve()
     * \return Vector x, solution of the matrix equation Ax = b.
     */
    arma::vec solveMatrixEquation() const;
};
