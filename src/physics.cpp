#include "simulator.hpp"

#include <cmath>
#include <vector>

#include "pipeline.hpp"
#include "solver/boundaryconditions.hpp"
#include "utilities/physics.hpp"
#include "equationofstate/equationofstate.hpp"
#include "equationofstate/equationofstatebase.hpp"
#include "equationofstate/bwrs.hpp"
#include "equationofstate/gerg04.hpp"
#include "equationofstate/idealgas.hpp"
#include "equationofstate/dummygas.hpp"
#include "heattransfer/heattransferbase.hpp"
#include "heattransfer/radial.hpp"
#include "heattransfer/unsteady.hpp"

using std::cout;
using std::endl;
using std::vector;
using std::make_unique;
using std::unique_ptr;
using arma::uword;
using arma::mat;
using arma::vec;

Physics::~Physics()
{}

Physics::Physics(
        const Pipeline& state,
        const Config& config):
    Physics(state, config.equationOfState, config.heatTransfer)
{}

Physics::Physics(
        const Pipeline& state,
        const std::string& eos,
        const std::string& heat):
    m_eos(make_unique<EquationOfState>(state, eos)),
    m_heat(make_unique<HeatTransfer>(state, heat))
{}

void Physics::updateDerivedProperties(Pipeline& state) const
{
    // calculate derived properties based on current pressure, temperature and composition

    // evaluate equation of state
    mat out = m_eos->evaluate(state);

    state.compressibilityFactor() = out.col(0);

    state.dZdtAtConstantPressure() = out.col(1);
    state.dZdpAtConstantTemperature() = out.col(2);
    state.dZdtAtConstantDensity() = out.col(3);

    state.heatCapacityConstantPressure() = out.col(4);
    state.heatCapacityConstantVolume() = out.col(5);

    state.molarMass() = out.col(6);

    // other stuff
    state.specificGasConstant() = constants::gasConstant/(state.molarMass()/1000.0 /*g->kg*/); // J/(kg*K)
    state.density() = state.pressure()/(state.compressibilityFactor() % state.specificGasConstant() % state.temperature());
    state.viscosity() = utils::calculateViscosity(state.molarMass(), state.temperature(), state.density());
    state.reynoldsNumber() = utils::calculateReynoldsNumber(state.flow(), state.diameter(), state.viscosity());
    state.velocity() = state.flow()/(state.density() % state.diameter());

    // this isn't technically a derived property, but...
    state.frictionFactor() =
            utils::calculateColebrookWhiteFrictionFactor(
                state.roughness(), state.diameter(), state.reynoldsNumber());
}

void Physics::initializeHeatTransferState(Pipeline& state) const
{
    vector<HeatTransferState> heatTransferState;
    for (uword i = 0; i < m_heat->size(); i++)
    {
        heatTransferState.push_back(
                    m_heat->at(i).makeState(
                        state.heatFlow()(i),
                        state.ambientTemperature()(i),
                        state.temperature()(i)));
    }

    state.heatTransferState() = heatTransferState;
    state.heatTransferIsInitialized() = true;
}

void Physics::thermalizeHeatTransfer(Pipeline& pipeline) const
{
    if (!pipeline.heatTransferIsInitialized())
        throw std::runtime_error("heat transfer not initialized");

    vector<HeatTransferState> heatTransferState;
    arma::vec heatFlow(m_heat->size());
    for (uword i = 0; i < m_heat->size(); i++)
    {
        // only possible for Unsteady
        const UnsteadyHeatTransfer* heat = dynamic_cast<const UnsteadyHeatTransfer*>(&m_heat->at(i));
        if (heat)
        {
            HeatTransferState state = heat->thermalizeToSteadyState(
                pipeline.ambientTemperature()(i),
                pipeline.pressure()(i),
                pipeline.temperature()(i),
                pipeline.reynoldsNumber()(i),
                pipeline.heatCapacityConstantPressure()(i),
                pipeline.viscosity()(i));
            heatTransferState.push_back(state);
            heatFlow(i) = state.heatFlux();
        }
        else
        {
            // just copy existing state if not Unsteady
            heatTransferState.push_back(pipeline.heatTransferState().at(i));
            heatFlow(i) = pipeline.heatTransferState().at(i).heatFlux();
        }
    }

    pipeline.heatTransferState() = heatTransferState;
    pipeline.heatFlow() = heatFlow;
}
