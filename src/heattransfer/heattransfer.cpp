#include "heattransfer.hpp"

#include "heattransferbase.hpp"
#include "heattransferstate.hpp"
#include "pipeline.hpp"
#include "unsteady.hpp"
#include "steadystate.hpp"
#include "fixedqvalue.hpp"
#include "fixeduvalue.hpp"

using std::make_unique;
using std::vector;
using std::unique_ptr;
using arma::uword;
using arma::vec;

// in this file only
unique_ptr<HeatTransferBase> makeSingle(
        const Pipeline& pipeline,
        const arma::uword i,
        const std::string& type);

HeatTransfer::~HeatTransfer()
{}

HeatTransfer::HeatTransfer(
        const Pipeline& pipeline,
        const std::string& type):
    m_heat(make_unique<vector<unique_ptr<HeatTransferBase>>>())
{
    // fill vector
    for (arma::uword i = 0; i < pipeline.size(); i++)
    {
        m_heat->push_back(makeSingle(pipeline, i, type));
    }
}

void HeatTransfer::evaluate(
        const std::vector<HeatTransferState>& state,
        const double timeStep,
        Pipeline& pipeline) const
{
    if (!pipeline.heatTransferIsInitialized())
    {
        throw std::runtime_error("heat transfer state is not initialized");
    }
    if (m_heat->size() != pipeline.size())
    {
        throw std::runtime_error("incompatible size)");
    }

    for (uword i = 0; i < m_heat->size(); i++)
    {
        HeatTransferState heatTransferState = m_heat->at(i)->evaluate(
                    state.at(i),
                    timeStep,
                    pipeline.ambientTemperature()(i),
                    pipeline.pressure()(i),
                    pipeline.temperature()(i),
                    pipeline.reynoldsNumber()(i),
                    pipeline.heatCapacityConstantPressure()(i),
                    pipeline.viscosity()(i));

        pipeline.heatFlow()(i) = heatTransferState.heatFlux();
        pipeline.heatTransferState().at(i) = heatTransferState;
    }
}

unique_ptr<HeatTransferBase> makeSingle(
        const Pipeline& pipeline,
        const arma::uword i,
        const std::string& type)
{
    if (type == "Unsteady")
    {
        return make_unique<UnsteadyHeatTransfer>(
                    pipeline.diameter()(i),
                    pipeline.pipeWall().at(i),
                    pipeline.burialDepth()(i),
                    pipeline.burialMedium().at(i),
                    pipeline.ambientFluid().at(i)
                    );
    }
    else if (type == "SteadyState")
    {
        return make_unique<SteadyStateHeatTransfer>(
                    pipeline.diameter()(i),
                    pipeline.pipeWall().at(i),
                    pipeline.burialDepth()(i),
                    pipeline.burialMedium().at(i),
                    pipeline.ambientFluid().at(i)
                    );
    }
    else if (type == "FixedQValue")
    {
        const double q = 0;
        return make_unique<FixedQValue>(q);
    }
    else if (type == "FixedUValue")
    {
        const double U = 0;
        return make_unique<FixedUValue>(U);
    }
    else
    {
        const std::string what = "invalid type transfer type \"" + type + "\"";
        throw std::invalid_argument(what);
    }
}
