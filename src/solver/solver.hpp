#pragma once

#include <string>
#include <memory>
#include <string>
#include <armadillo>

class Config;
class Physics;
class BatchTracking;
class Pipeline;
class GoverningEquationSolverBase;
class BoundaryConditions;
class BoundaryConditionsStamped;

/*!
 * \brief The Solver class combines GoverningEquationSolver and BatchTracking
 * to advance the governing equations forward in time, and advect the
 * gas composition.
 *
 * It operates on Pipeline instances, and returns a copy of the input, with
 * updated flow, pressure, temperature, and composition (if batch tracking is
 * enabled).
 */
class Solver
{
public:
    //! Declared to avoid the inline compiler-generated default destructor.
    virtual ~Solver();

    /*!
     * \brief Construct from Config.
     * \param nGridPoints Number of grid points
     * \param config Config instance
     */
    Solver(
            const arma::uword nGridPoints,
            const Config& config);

    /*!
     * \brief Construct from number of grid points and string to select
     * discretizer.
     * \param nGridPoints Number of grid points.
     * \param energyEquation Type of energy equation ("InternalEnergy" or "Enthalpy")
     * \param relaxationFactors Relaxation factors, typically around 1.
     * \param toleranceType Tolerance type ("absolute" or "relative")
     * \param tolerances Convergence criteria
     * \param bruteForce Will always do maxIterations iterations and not check
     * for convergence
     * \param maxIterations The maximum number of iterations to perform
     */
    explicit Solver(
            const arma::uword nGridPoints,
            const std::string& energyEquation = "InternalEnergy",
            const arma::vec& relaxationFactors = {1, 1, 2/3.0},
            const std::string& toleranceType = "relative",
            const arma::vec& tolerances = {0.001, 0.001, 0.001},
            const bool bruteForce = false,
            const arma::uword maxIterations = 200);

    /*!
     * \brief Solve the governing equations.
     *
     * \see solveWithIterations()
     *
     * \param dt Time step [s]
     * \param current Current pipeline state
     * \param boundaryConditions Boundary conditions with time stamp
     * \param physics Physics instance to update derived properties and heat transfer
     * \return Copy of current with updated properties.
     */
    virtual Pipeline solve(
        const arma::uword dt,
        const Pipeline &current,
        const BoundaryConditionsStamped &boundaryConditions,
        const Physics &physics) const;

    /*!
     * \brief Solve the governing equations.
     *
     * \see solveWithIterations()
     *
     * \param dt Time step [s]
     * \param current Current pipeline state
     * \param boundaryConditions Boundary conditions
     * \param physics Physics instance to update derived properties and heat transfer
     * \return Copy of current with updated properties.
     */
    virtual Pipeline solve(
            const arma::uword dt,
            const Pipeline& current,
            const BoundaryConditions& boundaryConditions,
            const Physics& physics) const;

    /*!
     * \brief Check if differences between an old and a new state are within
     * tolerances. This will either compare relative or absolute differences.
     *
     * This is used to see if the iteration over solution of the governing
     * equations have converged, by comparing differences between the two last
     * iterations.
     *
     * In the case of relative differences, we can run into issues if we have
     * zero flow, so this is handled explicitly (temperature and pressure is
     * never expected to reach zero in well-behaved simulations).
     *
     * \param guess
     * \param previous
     * \param tolerances
     * \param toleranceType
     * \param relaxationFactors
     * \return
     */
    static arma::uword differencesWithinTolerance(
            const Pipeline& guess,
            const Pipeline& previous,
            const arma::vec& tolerances,
            const std::string& toleranceType,
            const arma::vec& relaxationFactors);

    /*!
     * \brief Solve the governing equations.
     *
     * TODO: Some more documentation here
     *
     * \param dt Time step [s]
     * \param current Pipeline state
     * \param boundaryConditions Boundary conditions
     * \param physics Physics instance to update derived properties and heat transfer
     * \return Copy of current with updated properties.
     */
    Pipeline solveWithIterations(
            const arma::uword dt,
            const Pipeline& current,
            const BoundaryConditions& boundaryConditions,
            const Physics& physics) const;

    //! Enable brute force solver, which always does m_maxIterations iterations.
    void enableBruteForce();

    //! Set max iterations.
    void setMaxIterations(const arma::uword maxIterations);

    //! Return the number of iterations performed during previous solution attempt
    arma::uword nIterations() const { return m_nIterations; }

    //! Get (const ref) relaxation factors
    const arma::vec& relaxationFactors() const { return m_relaxationFactor;}
    //! Get (const ref) tolerance type
    const std::string& toleranceType() const { return m_toleranceType; }
    //! Get (const ref) tolerances
    const arma::vec& tolerances() const { return m_tolerances; }
    //! Get (const ref) GoverningEquationSolver
    const GoverningEquationSolverBase& governingEquationSolver() const { return *m_governingEquationSolver; }

private:
    //! Relaxation factors for each property (flow, pressure and temperature).
    arma::vec m_relaxationFactor;

    //! What kind of tolerance to use when checking convergence ("relative" or
    //! "absolute").
    std::string m_toleranceType;

    //! Convergence tolerances for each property (flow, pressure and temperature).
    arma::vec m_tolerances;

//    arma::uword m_maxNumberOfIterations = 20; // not used
//    bool m_hardLimitOnNumberOfIterations = true; // not used
//    bool m_bruteForce = false; // not used

    //! If we should use brute-force solving, which means not checking for
    //! convergence, but always performing m_maxIterations iterations.
    bool m_bruteForce;
    //! The maximum number of iterations to perform, or (if m_bruteForce is
    //! true) the exact number of iterations to perform.
    arma::uword m_maxIterations;
    //! The number of iterations performed during the previous solution attempt.
    //! This is mutable, and is updated in Solver::solve().
    mutable arma::uword m_nIterations = 0;

    //! GoverningEquationSolver instance used for solving the governing equations.
    std::unique_ptr<GoverningEquationSolverBase> m_governingEquationSolver;

    //! BatchTracking instance used to advect composition in the pipeline.
    std::unique_ptr<BatchTracking> m_compositionSolver;

    /*!
     * \brief Private method for making GoverningEquationSolverBase instance
     * from Config and nGridPoints.
     * \param config Config instance
     * \param nGridPoints Number of grid points
     * \return unique_ptr to GoverningEquationSolverBase
     */
    std::unique_ptr<GoverningEquationSolverBase> makeGoverningEquationSolver(
            const arma::uword nGridPoints,
            const Config& config);

    /*!
     * \brief Private method for making GoverningEquationSolverBase from string
     * and nGridPoints.
     * \param discretizer Discretizer type ("InternalEnergy" or "Enthalpy")
     * \param nGridPoints Number of grid points
     * \return unique_ptr to GoverningEquationSolverBase
     */
    std::unique_ptr<GoverningEquationSolverBase> makeGoverningEquationSolver(
            const arma::uword nGridPoints,
            const std::string& discretizer);
};
