#pragma once

#include <armadillo>
#include <memory>

class Pipeline;
class BoundaryConditions;
class InternalEnergyDiscretizer;
class MatrixEquation;
class Discretizer;

/*!
 * \brief Simple Base class to avoid having to specify template argument for
 * GoverningEquationSolver. This means that we have to use pointers to
 * GoverningEquationSolverBase instead of simple GoverningEquationSolverBase
 * instances.
 */
class GoverningEquationSolverBase
{
public:
    //! Declared to avoid the inline compiler-generated default destructor.
    virtual ~GoverningEquationSolverBase();

    //! See GoverningEquationSolver::solve().
    virtual arma::mat solve(
            const arma::uword dt,
            const Pipeline& currentState,
            const Pipeline& newState,
            const BoundaryConditions& boundaryConditions) = 0;

    //! Returns true of the equation system is over-determined given the
    //! input boundary conditions. Just returns true if there are more than
    //! three active boundary conditions, and false else.
    bool isOverDetermined(const BoundaryConditions& boundaryConditions);
};

/*!
 * \brief The GoverningEquationSolver class is a composition of Discretizer
 * and MatrixEquation with a little bit of logic, that solves the governing
 * equations for a 1d gas pipeline.
 */
template<typename T = InternalEnergyDiscretizer>
class GoverningEquationSolver : public GoverningEquationSolverBase
{
public:
    //! Declared to avoid the inline compiler-generated default destructor.
    virtual ~GoverningEquationSolver();

    /*!
     * \brief Construct with string determining the type of energy equation.
     * \param nGridPoints Number of grid points
     */
    explicit GoverningEquationSolver(const arma::uword nGridPoints);

    /*!
     * \brief Construct from Discretizer instance.
     * \param discretizer Discretizer instance
     */
    explicit GoverningEquationSolver(std::unique_ptr<Discretizer> discretizer);

    /*!
     * \brief Solve the governing equations for a given time step and boundary
     * conditions.
     * \param dt Time step [s]
     * \param currentState Current Pipeline state
     * \param newState New/guess Pipeline state
     * \param boundaryConditions Boundary conditions
     * \return Matrix containing flow, pressure and temperature columns
     */
    virtual arma::mat solve(
            const arma::uword dt,
            const Pipeline& currentState,
            const Pipeline& newState,
            const BoundaryConditions& boundaryConditions) override;

private:
    const arma::uword m_nVariables = 3; //!< Number of flow variables (flow, pressure and temperature)

    //! Discretizer used when solving the governing equations. The main
    //! difference between the two options InternalEnergyDiscretizer and
    //! EnthalpyDiscretizer is the choice of energy equation.
    //! \see InternalEnergyDiscretizer
    //! \see EnthalpyDiscretizer
    std::unique_ptr<Discretizer> m_discretizer;

    //! MatrixEquation instance used to solve the matrix equation set up when
    //! discretizing the governing equations.
    std::unique_ptr<MatrixEquation> m_matrixEquation;
};
