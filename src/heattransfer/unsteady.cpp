#include "heattransfer/unsteady.hpp"

#include <cmath>

#include "utilities/errors.hpp"
#include "utilities/numerics.hpp"
#include "heattransfer/utils.hpp"
#include "constants.hpp"

using arma::uword;
using arma::vec;
using arma::uvec;
using arma::zeros;
using std::pow;
using std::log;

UnsteadyHeatTransfer::UnsteadyHeatTransfer(
        const double diameter,
        const PipeWall& pipeWall,
        const double burialDepth,
        const BurialMedium& burialMedium,
        const AmbientFluid& ambientFluid):
    RadialHeatTransfer(diameter, pipeWall, burialDepth, burialMedium, ambientFluid)
{
    // Transfer coefficients between each shell in each pipe section
    //    m_heatTransferCoefficient; // == h
    m_heatTransferCoefficient = zeros<vec>(size());

    // Inner n-1 shells
    for (uword j = 0; j < size() - 1; j++)
    {
        m_heatTransferCoefficient(j) =
            pow( // Sum of ln(ro/ri)/lambda for half shells
                log(                      m_ro(j)/(m_ri(j) + m_width(j)/2.0) )/( 2.0*constants::pi*m_conductivity(j) )
                + log( (m_ro(j) + m_width(j+1)/2.0)/m_ro(j)                  )/( 2.0*constants::pi*m_conductivity(j+1) )
            , -1);
    }
    // Outermost shell
    uword j = size() - 1;
    m_heatTransferCoefficient(j) =
            pow(log( m_ro(j)/(m_ri(j) + m_width(j)/2.0) )/( 2.0*constants::pi*m_conductivity(j) ), -1);
}

// override
HeatTransferState UnsteadyHeatTransfer::evaluate(
        const HeatTransferState& current,
        const double timeStep,
        const double ambientTemperature,
        const double gasPressure,
        const double gasTemperature,
        const double gasReynoldsNumber,
        const double gasHeatCapacity,
        const double gasViscosity) const
{
    HeatTransferState newState = evaluateInternal(
                current.temperature(), timeStep, ambientTemperature,
                gasPressure, gasTemperature, gasReynoldsNumber, gasHeatCapacity,
                gasViscosity);

    return newState;
}

// private
HeatTransferState UnsteadyHeatTransfer::evaluateInternal(
        const vec& shellTemperature,
        const double timeStep,
        const double ambientTemperature,
        const double gasPressure,
        const double gasTemperature,
        const double gasReynoldsNumber,
        const double gasHeatCapacity,
        const double gasViscosity) const
{
    if (shellTemperature.n_elem != m_density.n_elem)
        throw std::runtime_error("incompatible size");

    vec x = solveEquations(
                shellTemperature, timeStep, gasPressure, gasTemperature,
                ambientTemperature, gasReynoldsNumber, gasHeatCapacity,
                gasViscosity);
    double heatFlux = x(0);
    vec newShellTemperature = x(arma::span(1, size()));

    if (arma::any(newShellTemperature < 0))
    {
        throw utils::temperature_range_error("wall layer temperature less than 0 K");
    }

    return HeatTransferState(heatFlux, newShellTemperature);
}

// private
vec UnsteadyHeatTransfer::solveEquations(
        const vec& shellTemperature,
        const double timeStep,
        const double pressure,
        const double temperature,
        const double ambientTemperature,
        const double reynoldsNumber,
        const double heatCapacity,
        const double viscosity
        ) const
{
    /* First do some global calculations ----------------------------------- */

    const double hi = utils::calcInnerWallFilmCoefficient(m_diameter, pressure, reynoldsNumber, heatCapacity, viscosity);
    const double hw = pow( // Heat transfer coefficient for inner film + half the first shell
        1.0/(2.0*constants::pi*m_ri(0)*hi)
        + log( (m_ri(0) + m_width(0)/2.0)/m_ri(0) )/( 2.0*constants::pi*m_conductivity(0) )
    , -1);

    vec mass = m_density%m_crossSection;
    vec factor = mass%m_heatCapacity/timeStep; // Factor used a lot in A and b

    /*
    Equation 2.4.1 in "Numerical Recipies" (2nd edition)
    Store tridiagonal matrix as 3 vectors instead of matrix, and solve with custom function tridag() from Numerical Recipies.

    | b1  c1   0  ...                 |   |  x1  |   |  r1  |
    | a2  b2  c2  ...                 |   |  x2  |   |  r2  |
    |             ...                 | * | ...  | = | ···  |
    |             ...  aN−1 bN−1 cN−1 |   | xN−1 |   | rN−1 |
    |             ...   0    aN   bN  |   |  xN  |   |  rN  |
    */
    vec at = zeros<vec>(size()+1);
    vec bt = zeros<vec>(size()+1);
    vec ct = zeros<vec>(size()+1);
    vec rt = zeros<vec>(size()+1);

    // note that a(0) and c(n) are undefined and are not referenced by the
    // solution routiune

    /* Setting up equations ------------------------------------------------ */
    // Heat equation, first part of eq. (2.26) in JFH thesis
    bt(0) = 1;
    ct(0) = hw/(constants::pi*m_diameter);
    rt(0) = temperature*hw/(constants::pi*m_diameter);

    // Pipe layer 1
    // this is the second line of eq. (2.26), multiplied by A_1
    bt(1) = factor(0) + hw + m_heatTransferCoefficient(0); // Could get this into loops below if we incoorporated hw in h
    ct(1) = -m_heatTransferCoefficient(0);
    rt(1) = factor(0)*shellTemperature(0) + hw*temperature;

    // Filling row 2:(n-1) of A and b using loops
    // this is the third line of eq. (2.26), multiplied by A_i
    for (uword i = 2; i < size(); i++)
    {
        bt(i) = m_heatTransferCoefficient(i-1) + m_heatTransferCoefficient(i-2) + factor(i-1); // Diagonal
        ct(i) = -m_heatTransferCoefficient(i-1); // Right of diagonal
        at(i) = -m_heatTransferCoefficient(i-2); // Left of diagonal
        rt(i) = factor(i-1)*shellTemperature(i-1);
    }

    // Calculate final heat transfer coefficient
    const uword end = size() - 1;
    const double k0N = calculateOuterFilmCoefficient();
    const double hN = pow(
        log( m_ro(end) / (m_ri(end) + m_width(end)/2.0))/(2*constants::pi*m_conductivity(end) )
        + 1.0/(m_ro(end)*2.0*constants::pi*k0N)
    , -1); // heat transfer coefficient for last, outer shell

    // Fill last row of A and b
    // this is the last line of eq. (2.26)
    uword i = size();
    bt(i) = hN + m_heatTransferCoefficient(i-2) + factor(i-1); // Diagonal
    at(i) = -m_heatTransferCoefficient(i-2); // Left of diagonal
    rt(i) = factor(i-1)*shellTemperature(i-1) + hN*ambientTemperature;

    // Solve Ax = b
    vec x = utils::tridag(at, bt, ct, rt, size()+1);
    return x;
}

HeatTransferState UnsteadyHeatTransfer::thermalizeToSteadyState(
        const double ambientTemperature,
        const double gasPressure,
        const double gasTemperature,
        const double gasReynoldsNumber,
        const double gasHeatCapacity,
        const double gasViscosity) const
{
    // thermalize until unsteady solution heat flux is equal to steady state heat flux
    // just update with infinite timestep
    // proper way would be to solve unsteady equation without time dependency, but no point in doing that when this works just fine

    // initial temp isn't important
    const arma::vec shellTemperature = arma::zeros<vec>(size()) + 273.15;
    HeatTransferState state = evaluateInternal(
                shellTemperature,
                std::numeric_limits<double>::infinity(), ambientTemperature,
                gasPressure, gasTemperature, gasReynoldsNumber, gasHeatCapacity,
                gasViscosity);

    return state;
}
