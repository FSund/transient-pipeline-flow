#include "solver/discretizer/enthalpy.hpp"

#include "pipeline.hpp"
#include "utilities/utilities.hpp"
#include "constants.hpp"

using arma::vec;
using arma::zeros;
using arma::span;
using arma::uword;

EnthalpyDiscretizer::EnthalpyDiscretizer(
        const uword nGridPoints):
    Discretizer(nGridPoints, 3)
{}

void EnthalpyDiscretizer::discretize(
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
                newState.heatCapacityConstantPressure(),
                newState.heatFlow(),
                newState.density(),
                newState.compressibilityFactor(),
                newState.dZdtAtConstantPressure(),
                newState.dZdpAtConstantTemperature());
}

void EnthalpyDiscretizer::discretizeFromPrimitives(
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
        const vec& guessHeatCapacityConstantPressure,
        const vec& guessHeatFlux,
        const vec& guessDensity,
        const vec& guessCompressibilityFactor,
        const vec& guess_dZdT_p,
        const vec& guess_dZdp)
{
    const uword nGridPoints = gridPoints.n_elem;
    const uword n = nGridPoints;

    const vec gasConstant           = utils::centerAverage(currentSpecificGasConstant);

    const vec diameter              = utils::centerAverage(diameter_);
    const vec crossSection          = constants::pi*pow(diameter/2.0, 2.0);
    const vec dh                    = utils::centerDifference(height); // difference
    const vec dx                    = utils::centerDifference(gridPoints); // difference

    const vec friction_             = utils::centerAverage(guessFriction);
    const vec heatCapacityAtConstantPressure_ = utils::centerAverage(guessHeatCapacityConstantPressure);
    const vec heatTransfer_         = utils::centerAverage(guessHeatFlux);
    const vec rho_                  = utils::centerAverage(guessDensity);

    const vec massFlow_             = utils::centerAverage(guessMassFlow);
    const vec pressure_             = utils::centerAverage(guessPressure);
    const vec temperature_          = utils::centerAverage(guessTemperature);

    const vec Z_                    = utils::centerAverage(guessCompressibilityFactor);
    const vec dZdT_p_               = utils::centerAverage(guess_dZdT_p);
    const vec dZdp_                 = utils::centerAverage(guess_dZdp);

    // indexing in term_i and term_ipp is (grid point, equation #, variable)
    // order of variables are m, p, T
    // order of equations are (continuity, momentum, energy)

    uword col; // For indexing

    // Common terms
    const vec ZRToverpA             = Z_%gasConstant%temperature_/(pressure_%crossSection);

    // Continuity equation
    col = 0;
    // Helgaker form
    const vec c1c                   = 1.0/(1.0/pressure_ - (1.0/Z_)%dZdp_);
    const vec c2c                   = 1.0/temperature_ + (1.0/Z_)%dZdT_p_;
    const vec c3c                   = ZRToverpA;
    m_term_i  .slice(0).col(col)    = - c1c%c3c/dx;
    m_term_ipp.slice(0).col(col)    = + c1c%c3c/dx;
    m_term_i  .slice(1).col(col)    = zeros<vec>(nGridPoints-1) + 1.0/(2.0*dt);
    m_term_ipp.slice(1).col(col)    = zeros<vec>(nGridPoints-1) + 1.0/(2.0*dt);
    m_term_i  .slice(2).col(col)    = - c1c%c2c/(2.0*dt);
    m_term_ipp.slice(2).col(col)    = - c1c%c2c/(2.0*dt);
    m_boundaryTerm.col(col) =
            - c1c%c2c%(currentTemperature(span(1, n-1)) + currentTemperature(span(0, n-2)))/(2.0*dt)
            + (currentPressure(span(1, n-1)) + currentPressure(span(0, n-2)))/(2.0*dt);

    // Momentum equation
    col = 1;
    // Helgaker form
    const vec c1m                   = massFlow_%ZRToverpA;
    const vec c2m                   = massFlow_%(1.0/pressure_ - (1.0/Z_)%dZdp_);
    const vec c3m                   = massFlow_%(1.0/temperature_ + (1.0/Z_)%dZdT_p_);
    const vec c4m                   = friction_%arma::abs(massFlow_)/(2.0*diameter)%ZRToverpA;
    const vec sinTheta              = dh/dx;
    const vec c5m                   = crossSection/(Z_%gasConstant%temperature_)*m_gravity%sinTheta;
    m_term_i  .slice(0).col(col)    = + ( 1.0/(2.0*dt) + c4m/2.0 ) - 2.0*c1m/dx;
    m_term_ipp.slice(0).col(col)    = + ( 1.0/(2.0*dt) + c4m/2.0 ) + 2.0*c1m/dx;
    m_term_i  .slice(1).col(col)    = - ( crossSection/dx - c1m%c2m/dx ) + c5m/2.0;
    m_term_ipp.slice(1).col(col)    = + ( crossSection/dx - c1m%c2m/dx ) + c5m/2.0;
    m_term_i  .slice(2).col(col)    = - c1m%c3m/dx;
    m_term_ipp.slice(2).col(col)    = + c1m%c3m/dx;
    m_boundaryTerm.col(col) =
            (currentMassFlow(span(1, n-1)) + currentMassFlow(span(0, n-2)))/(2.0*dt);

    // Energy equation
    // Common terms
    const vec oneMinusDZDp_T        = 1.0 - (pressure_/Z_)%dZdp_;
    const vec onePlusDZDT_p         = 1.0 + (temperature_/Z_)%dZdT_p_;
    const vec Vw2overZRT            = 1.0/( oneMinusDZDp_T - (Z_%gasConstant/heatCapacityAtConstantPressure_)%onePlusDZDT_p%onePlusDZDT_p );
    const vec Vw2overT              = Vw2overZRT%Z_%gasConstant;
    const vec Vw2                   = Vw2overZRT%Z_%gasConstant%temperature_;

    // convert from Q/A_h = U*(T-T_a) [W/m2] to to 4*U/D*(T-T_a)
    const vec actualHeatTransfer_   = 4*heatTransfer_/(diameter);

    col = 2;
    {
        const vec c1                    = onePlusDZDT_p;
        const vec c2                    = oneMinusDZDp_T;
        const vec c3                    = ZRToverpA;
        const vec c4                    = massFlow_%(1.0 + (Vw2overT/heatCapacityAtConstantPressure_)%c1%c1);
        const vec c5                    = 1.0/(heatCapacityAtConstantPressure_%pressure_)%c2;
        const vec c6                    = Vw2overT%actualHeatTransfer_;
        const vec c7                    = Vw2%massFlow_%arma::abs(massFlow_)%(friction_/(2.0*diameter%crossSection))%c3%c3;
        m_term_i  .slice(0).col(col)    = - Vw2/heatCapacityAtConstantPressure_%c1%c3/dx - c5%c7/2.0;
        m_term_ipp.slice(0).col(col)    = + Vw2/heatCapacityAtConstantPressure_%c1%c3/dx - c5%c7/2.0;
        m_term_i  .slice(1).col(col)    = + Vw2/heatCapacityAtConstantPressure_%c1%c3/pressure_%massFlow_%c2/dx;
        m_term_ipp.slice(1).col(col)    = - Vw2/heatCapacityAtConstantPressure_%c1%c3/pressure_%massFlow_%c2/dx;
        m_term_i  .slice(2).col(col)    = + 1.0/(2.0*dt) + c5%c6/2.0 - c3%c4/dx;
        m_term_ipp.slice(2).col(col)    = + 1.0/(2.0*dt) + c5%c6/2.0 + c3%c4/dx;
        m_boundaryTerm.col(col) =
                + (currentTemperature(span(1, n-1)) + currentTemperature(span(0, n-2)))/(2.0*dt);
    }
}
