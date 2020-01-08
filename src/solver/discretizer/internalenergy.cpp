#include "solver/discretizer/internalenergy.hpp"

#include <cmath>

#include "pipeline.hpp"
#include "utilities/utilities.hpp"
#include "constants.hpp"

using arma::vec;
using arma::uword;
using arma::span;
using arma::zeros;

InternalEnergyDiscretizer::InternalEnergyDiscretizer(
        const uword nGridPoints):
    Discretizer(nGridPoints, 3)
{}

void InternalEnergyDiscretizer::discretize(
        const arma::uword dt,
        const Pipeline& currentState,
        const Pipeline& newState)
{
    return discretizeFromPrimitives(
                dt,
                currentState.diameter(),
                currentState.height(),
                currentState.gridPoints(),

                currentState.specificGasConstant(),
                currentState.flow(),
                currentState.pressure(),
                currentState.temperature(),

                newState.flow(),
                newState.pressure(),
                newState.temperature(),
                newState.frictionFactor(),
                newState.heatCapacityConstantVolume(),
                newState.heatFlow(),
                newState.density(),
                newState.compressibilityFactor(),
                newState.dZdtAtConstantPressure(),
                newState.dZdpAtConstantTemperature(),
                newState.dZdtAtConstantDensity());
}

void InternalEnergyDiscretizer::discretizeFromPrimitives(
        const arma::uword dt,

        const vec& diameter_,
        const vec& height,
        const vec& gridPoints,

        const vec& currentSpecificGasConstant,
        const vec& currentMassFlow,
        const vec& currentPressure,
        const vec& currentTemperature,

        const vec& guessMassFlow,
        const vec& guessPressure,
        const vec& guessTemperature,
        const vec& guessFriction,
        const vec& guessHeatCapacityConstantVolume,
        const vec& guessHeatFlux,
        const vec& guessDensity,
        const vec& guessCompressibilityFactor,
        const vec& guess_dZdT_p,
        const vec& guess_dZdp,
        const vec& guess_dZdT_rho)
{
    const uword nGridPoints = gridPoints.n_elem;
    const uword n = nGridPoints;

    const vec gasConstant           = utils::centerAverage(currentSpecificGasConstant);

    const vec diameter              = utils::centerAverage(diameter_);
    const vec crossSection          = constants::pi*arma::pow(diameter/2.0, 2.0);
    const vec dh                    =  utils::centerDifference(height); // difference
    const vec dx                    =  utils::centerDifference(gridPoints); // difference

    const vec friction_             = utils::centerAverage(guessFriction);
    const vec heatCapacity_         = utils::centerAverage(guessHeatCapacityConstantVolume);
    const vec rho_                  = utils::centerAverage(guessDensity);

    // convert from Q/A_h = U*(T-T_a) [W/m2] to to -4*U/(D*rho)*(T-T_a)
    const vec q                     = utils::centerAverage(guessHeatFlux);
    const vec heatTransfer_         = -4.0*q/(diameter%rho_);

    const vec massFlow_             = utils::centerAverage(guessMassFlow);
    const vec pressure_             = utils::centerAverage(guessPressure);
    const vec temperature_          = utils::centerAverage(guessTemperature);

    const vec Z_                    = utils::centerAverage(guessCompressibilityFactor);
    const vec dZdT_p_               = utils::centerAverage(guess_dZdT_p);
    const vec dZdp_                 = utils::centerAverage(guess_dZdp);
    const vec dZdT_rho_             = utils::centerAverage(guess_dZdT_rho);

    // indexing in term_i and term_ipp is (grid point, equation #, variable)
    // order of variables are m, p, T
    // order of equations are (continuity, momentum, energy)

    uword col; // For indexing

    // Continuity equation
    col = 0;
    const vec c1c = 1.0/(1.0/pressure_ - (1.0/Z_)%dZdp_);
    const vec c2c = 1.0/temperature_ + (1.0/Z_)%dZdT_p_;
    const vec c3c = Z_%gasConstant%temperature_/(pressure_%crossSection);
    m_term_i  .slice(0).col(col) = - c1c%c3c/dx;
    m_term_ipp.slice(0).col(col) = + c1c%c3c/dx;
    m_term_i  .slice(1).col(col) = zeros<vec>(nGridPoints-1) + 1.0/(2.0*dt);
    m_term_ipp.slice(1).col(col) = zeros<vec>(nGridPoints-1) + 1.0/(2.0*dt);
    m_term_i  .slice(2).col(col) = - c1c%c2c/(2.0*dt);
    m_term_ipp.slice(2).col(col) = - c1c%c2c/(2.0*dt);
    m_boundaryTerm.col(col) =
            - c1c%c2c%(currentTemperature(span(1, n-1)) + currentTemperature(span(0, n-2)))/(2.0*dt)
            + (currentPressure(span(1, n-1)) + currentPressure(span(0, n-2)))/(2.0*dt);

    // Momentum equation
    col = 1;
    const vec c1m = massFlow_%Z_%gasConstant%temperature_/(pressure_%crossSection);
    const vec c2m = massFlow_%(1.0/pressure_ - (1.0/Z_)%dZdp_);
    const vec c3m = massFlow_%(1.0/temperature_ + (1.0/Z_)%dZdT_p_);
    const vec c4m = friction_%Z_%gasConstant%temperature_%arma::abs(massFlow_)/(2.0*diameter%crossSection%pressure_);
    const vec sinTheta              = dh/dx;
    const vec c5m = crossSection/(Z_%gasConstant%temperature_)*m_gravity%sinTheta;
    m_term_i  .slice(0).col(col) = + ( 1.0/(2.0*dt) + c4m/2.0 ) - 2.0*c1m/dx; // default
    m_term_ipp.slice(0).col(col) = + ( 1.0/(2.0*dt) + c4m/2.0 ) + 2.0*c1m/dx; // default
    m_term_i  .slice(1).col(col) = - ( crossSection/dx - c1m%c2m/dx ) + c5m/2.0; // default
    m_term_ipp.slice(1).col(col) = + ( crossSection/dx - c1m%c2m/dx ) + c5m/2.0; // default
    m_term_i  .slice(2).col(col) = - c1m%c3m/dx; // default - ok
    m_term_ipp.slice(2).col(col) = + c1m%c3m/dx; // default - ok
    m_boundaryTerm.col(col) =
            (currentMassFlow(span(1, n-1)) + currentMassFlow(span(0, n-2)))/(2.0*dt);

    // Energy equation
    col = 2;
    const vec c1e = massFlow_%Z_%gasConstant%temperature_/(pressure_%crossSection);
    const vec c2e = c1e%Z_%gasConstant%temperature_/heatCapacity_%temperature_%(1.0/temperature_ + (1.0/Z_)%dZdT_rho_);
    const vec c2eoverm = Z_%gasConstant%temperature_/(pressure_%crossSection)  %  Z_%gasConstant%temperature_/heatCapacity_%temperature_%(1.0/temperature_ + (1.0/Z_)%dZdT_rho_);
    const vec c3e = 1.0/pressure_ - (1.0/Z_)%dZdp_;
    const vec c4e = 1.0/temperature_ + (1.0/Z_)%dZdT_p_;
    const vec c5e =
            + friction_/(2.0*heatCapacity_%diameter)
                %arma::pow(temperature_, 2.0)
                %arma::pow(Z_%gasConstant%massFlow_/(pressure_%crossSection), 3.0)
            + 1.0/(temperature_%heatCapacity_)%heatTransfer_; // heatTransfer = -4*U/(D*rho)*(T-T_a)
    m_term_i  .slice(0).col(col) = - c2eoverm/dx; // default
    m_term_ipp.slice(0).col(col) = + c2eoverm/dx; // default
    m_term_i  .slice(1).col(col) = + c2e%c3e/dx; // default
    m_term_ipp.slice(1).col(col) = - c2e%c3e/dx; // default
    m_term_i  .slice(2).col(col) = + ( 1.0/(2.0*dt) - c5e/2.0 ) - ( c1e/dx + c2e%c4e/dx ); // default
    m_term_ipp.slice(2).col(col) = + ( 1.0/(2.0*dt) - c5e/2.0 ) + ( c1e/dx + c2e%c4e/dx ); // default
    m_boundaryTerm.col(col) =
            + (currentTemperature(span(1, n-1)) + currentTemperature(span(0, n-2)))/(2.0*dt);
}
