#pragma once

#include <string>
#include <vector>
#include <memory>
#include <optional>
#include <armadillo>

// can replace some of these with forward declarations, but keep them for now
// so the program is easier to use (no need to include so much stuff)
// TODO: should probably make a "transflow.hpp" file with the most useful stuff
// included
#include "pipeline.hpp"
#include "solver/solver.hpp"
#include "physics.hpp"
#include "timeseries.hpp"
#include "config.hpp"
#include "sampler.hpp"

/*!
 * \brief The Simulator class combines Physics and Solver to advance the state
 * of the pipeline in time. It contains a Pipeline instance which has thes state
 * of the pipeline at all times.
 */
class Simulator
{
public:
    /*!
     * \brief Construct using Pipeline and Config. Will construct Physics and
     * Solver using the settings in Config, and will make local copy of
     * pipeline.
     * \param pipeline Pipeline description
     * \param config Simulator configuration
     */
    Simulator(
            const Pipeline& pipeline = Pipeline(),
            const Config& config = Config());

    /*!
     * \brief Construct from finished instances of Pipeline, Physics and Solver.
     * \param pipeline Pipeline description
     * \param physics Physics description
     * \param solver Solver instance
     */
    Simulator(
            const Pipeline& pipeline,
            std::unique_ptr<Physics>& physics,
            std::unique_ptr<Solver>& solver);

    /*!
     * \brief Advance the pipeline in time by any number of time steps. If
     * initialize is true the pipeline pressure, temperature and flow will be
     * set to reasonable values according to the boundary conditions before
     * the simulation is started.
     * \param timeSeries Boundary conditions with timestamps
     * \return Number of iterations at each time step
     */
    arma::vec simulate(const TimeSeries& timeSeries);

    //! Enables batch tracking. Wrapper around Pipeline::enableBatchTracking().
    void enableBatchTracking();

    //! Get size (number of grid points)
    arma::uword size() const { return m_state->size(); }

    //! Get (const ref) Pipeline
    const Pipeline& pipeline() const { return *m_state; }

    //! Get (const ref) Pipeline::State
    const Pipeline::State& state() const { return m_state->state(); }

    //! Get (const ref) Solver
    const Solver& solver() const { return *m_solver; }

    //! Get (const ref) Physics
    const Physics& physics() const { return *m_physics; }

    //! Get (ref) sampler
    Sampler& sampler() { return m_sampler.value(); }

private:
    std::unique_ptr<Pipeline> m_state; //!< Pipeline state
    std::unique_ptr<Physics> m_physics; //!< Physics instance, contains HeatTransfer and EquationOfState
    std::unique_ptr<Solver> m_solver; //!< Solver instance, contains GoverningEquationSolver and BatchTracking
    std::optional<Sampler> m_sampler; //!< (Optional) sampler instance, for writing results to file during simulation

    /*!
     * \brief Make optional Sampler instance. Returns empty optional if
     * config.outputPath is empty.
     * \return config Config instance
     */
    std::optional<Sampler> makeSampler(const Config& config);
};
