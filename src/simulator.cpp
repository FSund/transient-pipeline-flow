#include "simulator.hpp"

#include <cmath>
#include <stdexcept>

#include "solver/discretizer/enthalpy.hpp"
#include "solver/discretizer/internalenergy.hpp"
#include "solver/boundaryconditions.hpp"
#include "solver/solver.hpp"
#include "timeseries.hpp"

Simulator::Simulator(const Pipeline& pipeline, const Config& config):
    m_state(std::make_unique<Pipeline>(pipeline)),
    m_physics(std::make_unique<Physics>(*m_state, config)),
    m_solver(std::make_unique<Solver>(pipeline.size(), config)),
    m_sampler(makeSampler(config)) // optional
{
    m_physics->updateDerivedProperties(*m_state);

    m_physics->initializeHeatTransferState(*m_state);
    m_physics->thermalizeHeatTransfer(*m_state);

    m_state->initializeBatchTracking();
}

Simulator::Simulator(
        const Pipeline& pipeline,
        std::unique_ptr<Physics>& physics,
        std::unique_ptr<Solver>& solver):
    m_state(std::make_unique<Pipeline>(pipeline)),
    m_physics(std::move(physics)),
    m_solver(std::move(solver))
{
    m_physics->updateDerivedProperties(*m_state);

    m_physics->initializeHeatTransferState(*m_state);
    m_physics->thermalizeHeatTransfer(*m_state);

    m_state->initializeBatchTracking();
}

void Simulator::enableBatchTracking()
{
    m_state->enableBatchTracking();
}

std::optional<Sampler> Simulator::makeSampler(const Config& config)
{
    if (config.outputPath == "")
    {
        return {};
    }
    else
    {
        return Sampler(config.outputPath, config.samplingInterval, config.appendResults);
    }
}

arma::vec Simulator::simulate(const TimeSeries& ts)
{
    const std::vector<TimeStep>& timeSeries = ts; // user-defined cast
    arma::vec nIterations(ts.size());

    if (m_sampler && m_state->timestamp() == 0)
    {
        m_sampler->sample(*m_state);
    }

    arma::uword dt = 0;
    for (std::size_t i = 0; i < timeSeries.size(); i++)
    {
        const TimeStep& bc = timeSeries.at(i);

        // timestamps of m_state and boundary conditions should now be synced
        if (bc.timestamp() < m_state->timestamp())
            throw std::runtime_error("negative time step, likely error with timestamps");

        if (bc.timestamp() - m_state->timestamp() > 24*60*60)
            throw std::runtime_error("time step larger than 24 hours, likely error with timestamps");

        // calculate time step
        dt = bc.timestamp() - m_state->timestamp();

        if (dt == 0)
        {
            nIterations(i) = 0;
            continue; // skip
        }

        *m_state = m_solver->solve(dt, *m_state, bc, *m_physics);
        m_state->timestamp() = bc.timestamp(); // sync with boundary conditions

        nIterations(i) = m_solver->nIterations();

        if (m_sampler)
        {
            m_sampler->sample(*m_state);
        }
    }

    return nIterations;
}
