#include "pipeline.hpp"

#include "constants.hpp"
#include "heattransfer/heattransferbase.hpp"
#include "solver/boundaryconditions.hpp"
#include "advection/batchtracking.hpp"
#include "advection/batchtrackingstate.hpp"

using arma::zeros;
using arma::vec;
using std::vector;

Pipeline::Pipeline(const arma::uword size, const double length):
    m_length(length),
    m_gridPoints(arma::linspace(0, length, size)),
    m_diameter(zeros<vec>(size) + 1),
    m_height(zeros<vec>(size)),
    m_roughness(zeros<vec>(size)),

    m_burialDepth(zeros<vec>(size)),
    m_pipeWall(size, PipeWall::defaultPipeWall),
    m_burialMedium(size, BurialMedium::soil),
    m_ambientFluid(size, AmbientFluid::seawater),

    m_constantComposition(true),
    m_state(
        m_gridPoints,
        zeros<vec>(size) + constants::standardPressure,
        zeros<vec>(size) + constants::standardTemperature,
        zeros<vec>(size), // flow
        Composition::defaultComposition)
{
    ambientTemperature() = temperature();
}

BoundaryConditions Pipeline::getBoundaryConditions() const
{
    return BoundaryConditions(*this); // use BoundaryConditions copy constructor to avoid code duplication
}

void Pipeline::setLength(const double length)
{
    m_length = length;
    m_gridPoints = arma::linspace(0, length, size());
    prop().m_batchTrackingState = BatchTrackingState(m_gridPoints, prop().m_composition);
}

const arma::vec& Pipeline::inletComposition() const
{
    return prop().composition().front().vec();
}

const arma::vec& Pipeline::outletComposition() const
{
    return prop().composition().back().vec();
}

void Pipeline::updateComposition(const std::vector<Composition>& composition)
{
    prop().m_composition = composition;
    initializeBatchTracking();
}

void Pipeline::updateComposition(const Composition& composition)
{
    updateComposition(std::vector<Composition>(size(), composition));
}

void Pipeline::setCompositionUnsafe(const std::vector<Composition>& composition)
{
    prop().m_composition = composition;
}

void Pipeline::enableBatchTracking()
{
    m_constantComposition = false;
    initializeBatchTracking();
}

void Pipeline::initializeBatchTracking()
{
    prop().m_batchTrackingState = BatchTrackingState(m_gridPoints, prop().m_composition);
    prop().m_batchTrackingIsInitialized = true;
}


////////////////////////////////////////////////////////////////////////////////
Pipeline::State::State(
        const arma::vec& gridPoints,
        const arma::vec& pressure,
        const arma::vec& temperature,
        const arma::vec& flow,
        const std::vector<Composition>& composition):
    Pipeline::State(gridPoints)
{
    if (pressure.n_elem != temperature.n_elem
            || temperature.n_elem != flow.n_elem
            || flow.n_elem != composition.size())
        throw std::invalid_argument("incompatible size");

    m_pressure = pressure;
    m_temperature = temperature;
    m_flow = flow;
    m_composition = composition;

    m_batchTrackingState = BatchTrackingState(gridPoints, composition); // fill with correct composition
    m_batchTrackingIsInitialized = true;
}

Pipeline::State::State(
        const arma::vec& gridPoints,
        const arma::vec& pressure,
        const arma::vec& temperature,
        const arma::vec& flow,
        const Composition& composition):
    Pipeline::State(
        gridPoints, pressure, temperature, flow,
        vector<Composition>(gridPoints.n_elem, composition))
{}

std::ostream& operator <<(std::ostream& out, const Pipeline::State& state)
{
    const arma::uword end = state.flow().size() - 1;
    out << "Flow:        " << state.flow()(0) << ", " << state.flow()(end) << std::endl;
    out << "Pressure:    " << state.pressure()(0) << ", " << state.pressure()(end) << std::endl;
    out << "Temperature: " << state.temperature()(0) << ", " << state.temperature()(end) << std::endl;
    out << "Inlet comp:  " << state.m_composition.front().vec().t();
    out << "Outlet comp: " << state.m_composition.back().vec().t();

    return out;
}

// protected
Pipeline::State::State(const arma::vec& gridPoints):
    m_flow(zeros<vec>(gridPoints.n_elem)),
    m_pressure(zeros<vec>(gridPoints.n_elem)),
    m_temperature(zeros<vec>(gridPoints.n_elem)),
    m_composition(gridPoints.n_elem),

    m_heatCapacityConstantVolume(zeros<vec>(gridPoints.n_elem)),
    m_heatCapacityConstantPressure(zeros<vec>(gridPoints.n_elem)),
    m_density(zeros<vec>(gridPoints.n_elem)),
    m_viscosity(zeros<vec>(gridPoints.n_elem)),
    m_specificGasConstant(zeros<vec>(gridPoints.n_elem)),
    m_molarMass(zeros<vec>(gridPoints.n_elem)),

    m_compressibilityFactor(zeros<vec>(gridPoints.n_elem)),
    m_temperatureDerivativeConstantPressure(zeros<vec>(gridPoints.n_elem)),
    m_pressureDerivativeConstantTemperature(zeros<vec>(gridPoints.n_elem)),
    m_temperatureDerivativeConstantDensity(zeros<vec>(gridPoints.n_elem)),

    m_velocity(zeros<vec>(gridPoints.n_elem)),
    m_frictionFactor(zeros<vec>(gridPoints.n_elem)),
    m_reynoldsNumber(zeros<vec>(gridPoints.n_elem)),

    m_ambientTemperature(zeros<vec>(gridPoints.n_elem)),
    m_heatFlow(zeros<vec>(gridPoints.n_elem)),

    m_heatTransferState(gridPoints.n_elem, HeatTransferState()),
    m_heatTransferIsInitialized(false),

    m_batchTrackingState(gridPoints, m_composition),
    m_batchTrackingIsInitialized(false)
{}
////////////////////////////////////////////////////////////////////////////////
