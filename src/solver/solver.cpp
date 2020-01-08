#include "solver/solver.hpp"

#include <cmath>

#include "advection/batchtracking.hpp"
#include "advection/batchtrackingstate.hpp"
#include "solver/discretizer/enthalpy.hpp"
#include "solver/discretizer/internalenergy.hpp"
#include "solver/governingequationsolver.hpp"
#include "utilities/errors.hpp"
#include "utilities/stringbuilder.hpp"
#include "boundaryconditions.hpp"
#include "config.hpp"
#include "physics.hpp"
#include "pipeline.hpp"

using arma::uword;
using arma::vec;
using arma::uvec;
using arma::mat;
using arma::zeros;
using std::cout;
using std::endl;

Solver::~Solver()
{}

Solver::Solver(
        const arma::uword nGridPoints,
        const std::string& energyEquation,
        const arma::vec& relaxationFactors,
        const std::string& toleranceType,
        const arma::vec& tolerances,
        const bool bruteForce,
        const arma::uword maxIterations):
    m_relaxationFactor(relaxationFactors),
    m_toleranceType(toleranceType),
    m_tolerances(tolerances),
    m_bruteForce(bruteForce),
    m_maxIterations(maxIterations),
    m_governingEquationSolver(makeGoverningEquationSolver(nGridPoints, energyEquation)),
    m_compositionSolver(std::make_unique<BatchTracking>())
{}

Solver::Solver(const arma::uword nGridPoints, const Config& config):
    Solver(nGridPoints, config.discretizer, config.relaxationFactors,
           config.toleranceType, config.tolerances, config.bruteForce,
           config.maxIterations)
{}

Pipeline Solver::solve(
        const arma::uword dt,
        const Pipeline& current,
        const TimeStep& boundaryConditions,
        const Physics& physics) const
{
    return solve(dt, current, BoundaryConditions(boundaryConditions), physics);
}

Pipeline Solver::solve(
        const arma::uword dt,
        const Pipeline& current,
        const BoundaryConditions& boundaryConditions,
        const Physics& physics) const
{
    try
    {
        return solveWithIterations(dt, current, boundaryConditions, physics);
    }
    catch (const utils::no_convergence_error& e)
    {
        cout << e.what() << endl;
    }

    // if all else fails, return old state
    return current;
}

Pipeline Solver::solveWithIterations(
        const arma::uword dt,
        const Pipeline& current,
        const BoundaryConditions& boundaryConditions,
        const Physics& physics) const
{
    Pipeline guess = current; // make copy
    Pipeline previous = current;

    arma::vec relaxationFactor = m_relaxationFactor;

    const bool lowFlowState =
            (boundaryConditions.inletFlow().isActive() && boundaryConditions.inletFlow() < 10)
            || (boundaryConditions.outletFlow().isActive() && boundaryConditions.outletFlow() < 10);

    m_nIterations = 0; // mutable
    while (true)
    {
        m_nIterations++;

        // decrease relaxation factor after some iterations
        if (m_nIterations >= 50)
            relaxationFactor *= 0.95;

        // calculate new flow, pressure, temperature
        const mat output = m_governingEquationSolver->solve(dt, current, guess, boundaryConditions);
        guess.flow()        = guess.flow()        + (output.col(0) - guess.flow())       *relaxationFactor(0);
        guess.pressure()    = guess.pressure()    + (output.col(1) - guess.pressure())   *relaxationFactor(1);
        guess.temperature() = guess.temperature() + (output.col(2) - guess.temperature())*relaxationFactor(2);

        // TODO: do some validation of flow, pressure and temperature here?

        if (!guess.constantComposition())
        {
            if (!guess.batchTrackingIsInitialized())
                throw std::runtime_error("batch tracking not initialized");

            // calculate new composition
            // here it's important to use the original ("current") batch
            // tracking state and the new ("guess") velocity, or else we risk
            // advecting the batches each iteration
            guess.batchTrackingState() = m_compositionSolver->advect(current.batchTrackingState(), dt, guess, boundaryConditions);
            guess.setCompositionUnsafe(guess.batchTrackingState().sample());
        }

        // update derived properties and heat transfer
        physics.updateDerivedProperties(guess);
        physics.heatTransfer().evaluate(current.heatTransferState(), dt, guess); // updates heat transfer state and heat flow

        if (!m_bruteForce)
        {
            // stop if converged
            const bool converged = differencesWithinTolerance(guess, previous, m_tolerances, m_toleranceType, relaxationFactor);
            if (converged)
                break;

            // this solver has issues converging at low flows, especially if using
            // unsteady heat transfer
            // hotfix is to just limit the number of iterations in these cases,
            // and return even if we haven't converged
            if (lowFlowState && m_nIterations >= 5)
                break;
        }

        if (m_nIterations >= m_maxIterations)
            break;

        // for comparing differences
        previous = guess;
    }

    if (!m_bruteForce)
    {
        if (m_nIterations >= m_maxIterations)
        {
            throw utils::no_convergence_error(utils::stringbuilder() << "no convergence after " << m_nIterations << " iterations", m_nIterations);
        }
    }

    return guess;
}

void Solver::enableBruteForce()
{
    m_bruteForce = true;
}

void Solver::setMaxIterations(const arma::uword maxIterations)
{
    m_maxIterations = maxIterations;
}

std::unique_ptr<GoverningEquationSolverBase> Solver::makeGoverningEquationSolver(
        const arma::uword nGridPoints,
        const Config& config)
{
    return makeGoverningEquationSolver(nGridPoints, config.discretizer);
}

std::unique_ptr<GoverningEquationSolverBase> Solver::makeGoverningEquationSolver(
        const arma::uword nGridPoints,
        const std::string& discretizer)
{
    if (discretizer == "InternalEnergy")
    {
        return std::make_unique<GoverningEquationSolver<InternalEnergyDiscretizer>>(nGridPoints);
    }
    else if (discretizer == "Enthalpy")
    {
        return std::make_unique<GoverningEquationSolver<EnthalpyDiscretizer>>(nGridPoints);
    }

    throw std::invalid_argument("discretizer");
}

uword Solver::differencesWithinTolerance(
        const Pipeline& guess,
        const Pipeline& previous,
        const vec& tolerances,
        const std::string& toleranceType,
        const vec& relaxationFactors)
{
    vec flowDiff;
    vec pressureDiff;
    vec temperatureDiff;
    if (toleranceType == "absolute")
    {
        flowDiff = arma::abs(guess.flow() - previous.flow());
        pressureDiff = arma::abs(guess.pressure() - previous.pressure());
        temperatureDiff = arma::abs(guess.temperature() - previous.temperature());
    }
    else if (toleranceType == "relative")
    {
        flowDiff = zeros<vec>(guess.size());
        for (uword i = 0; i < guess.size(); i++)
        {
            const double currentAbsFlow = std::abs(guess.flow()(i));
            if (currentAbsFlow > 10*tolerances(0))
            {
                // use relative difference
                flowDiff(i) = std::abs(guess.flow()(i) - previous.flow()(i))/std::abs(previous.flow()(i));
            }
            else
            {
                // if actual flow is lower than 10*tolerances(0), don't check error
                flowDiff(i) = 0.0;
            }
        }
        pressureDiff = arma::abs(guess.pressure() - previous.pressure())/previous.pressure();
        temperatureDiff = arma::abs(guess.temperature() - previous.temperature())/previous.temperature();
    }
    else
    {
        throw std::invalid_argument(utils::stringbuilder() << "unknown tolerance type \"" << toleranceType << "\"");
    }

    flowDiff = flowDiff/relaxationFactors(0);
    pressureDiff = pressureDiff/relaxationFactors(1);
    temperatureDiff = temperatureDiff/relaxationFactors(2);

    uvec withinTolerance = zeros<uvec>(3);
    withinTolerance(0) = arma::sum(flowDiff > tolerances(0)) == 0 ? true : false;
    withinTolerance(1) = arma::sum(pressureDiff > tolerances(1)) == 0 ? true : false;
    withinTolerance(2) = arma::sum(temperatureDiff > tolerances(2)) == 0 ? true : false;

    return arma::all(withinTolerance);
}
