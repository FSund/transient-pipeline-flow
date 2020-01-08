#pragma once

#include <string>
#include <vector>
#include <armadillo>

/*!
 * \brief The Config struct stores all settings for the different parts of the
 * pipeline simulation.
 */
struct Config
{
    //! Equation of state, either "BWRS", "GERG04" or "IdealGas"
    std::string equationOfState = "BWRS";
    //! Type of heat transfer, either "SteadyState", "Unsteady", "FixedUValue" or "FixedQValue"
    std::string heatTransfer = "SteadyState";
    //! Type of energy equation, either "InternalEnergy" or "Enthalpy"
    std::string discretizer = "InternalEnergy";

    //! The tolerance of Solver, either "relative" or "absolute"
    std::string toleranceType = "relative";
    //! The relaxation factors used by Solver.
    arma::vec relaxationFactors = {1, 1, 2/3.0};
    //! The convergence criteria used by Solver.
    arma::vec tolerances = {0.001, 0.001, 0.001};
    //! If we want to force the number of iterations in Solver::solve
    bool bruteForce = false;
    //! Max number of iterations to use in Solver::solve
    arma::uword maxIterations = 200;

    //! Where to put results. Will not output any results if equal to empty
    //! string.
    std::string outputPath = ""; // no output by default

    //! How often to sample results [s]
    arma::uword samplingInterval = 60;

    //! If we should append to (true) or overwrite (false) existing output files
    bool appendResults = false;
};
