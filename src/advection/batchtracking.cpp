#include "batchtracking.hpp"

#include <cmath>

#include "pipeline.hpp"
#include "solver/boundaryconditions.hpp"
#include "utilities/utilities.hpp"
#include "batchtrackingstate.hpp"

using arma::mat;
using arma::vec;
using arma::sword;
using arma::uvec;
using arma::uword;
using arma::zeros;
using std::vector;
using std::cout;
using std::endl;

/*!
 * \ingroup typedefs
 * Shorthand for BatchTrackingState::Batch.
 */
typedef BatchTrackingState::Batch Batch;

BatchTracking::~BatchTracking()
{}

BatchTrackingState BatchTracking::advect(
        const BatchTrackingState& state,
        const arma::uword dt,
        const Pipeline& pipeline,
        const BoundaryConditions& boundaryConditions)
{
    // check that pipeline and batch tracking state are consistent
    // (check that all batches are within the grid points)
    if (pipeline.batchTrackingState().batches().front().position() < pipeline.gridPoints()(0))
        throw std::runtime_error("first batch outside grid");
    if (pipeline.batchTrackingState().batches().back().position() > pipeline.gridPoints().tail(1)(0))
        throw std::runtime_error("last batch outside grid");

    if (!pipeline.batchTrackingIsInitialized())
        throw std::runtime_error("batch tracking not initialized");

    const vec velocity = utils::centerAverage(pipeline.velocity());
    mat inletAndOutletComposition = arma::join_horiz(
                vec(boundaryConditions.inletComposition()),
                vec(boundaryConditions.outletComposition()));

    return advect(state, dt, inletAndOutletComposition, velocity);
}

// keep this version using inletAndOutletComposition so we can test with simpler
// compositions (with only one component for example)
BatchTrackingState BatchTracking::advect(
        const BatchTrackingState& state,
        const arma::uword dt,
        const arma::mat& inletAndOutletComposition,
        const arma::vec& velocity)
{
    // TODO: Implement fix for negative velocity. Have in mind that velocity might have varying sign along the pipeline
    if (arma::any(velocity < 0))
    {
        cout << "WARNING: BatchTracker::advect(): Negative velocity (not implemented), so skipping advection." << endl;
        return state; // makes copy
    }

    // skip advection if all velocities are zero
    if (arma::all(velocity == 0))
    {
        return state;
    }

    // we use centered velocities as input, so velocity should be one larger
    // than the number of grid points
    if (velocity.n_elem != state.m_gridPoints.n_elem - 1)
    {
        throw std::invalid_argument("inconsistent sizes (velocity.n_elem != state.m_gridPoints.n_elem - 1)");
    }

    BatchTrackingState newState(state);

    // shorthand
    const arma::vec& gridPoints = newState.m_gridPoints;
    std::vector<Batch>& batches = newState.m_batches;

    // do a backwards loop, so we can do m_batches.pop_back() if the last batch
    // moves outside domain - this also works if several batches move outside
    // the domain
    for (sword i = sword(batches.size()) - 1; i >= 0; i--)
    {
        auto& batch = batches.at(uword(i));

        // find where on grid this batch is
        // find last grid point smaller than batch position
        uvec tmp = arma::find(gridPoints <= batch.m_position, 1, "last");
        if (tmp.n_elem == 0)
        {
            throw std::runtime_error("batch outside grid (tmp.n_elem > 0)");
        }

        uword j = tmp(0);

        if (j+1 >= gridPoints.n_elem || j >= velocity.n_elem)
        {
            throw std::runtime_error("invalid index");
        }

        // move this batch according to velocity
        double timeTravelled = 0;
        double distanceTravelled = 0;
        while (timeTravelled < dt)
        {
            double distanceToEndOfCell = gridPoints(j+1) - batch.m_position;
            double maxTimeInThisCell = distanceToEndOfCell/velocity(j);
            if (maxTimeInThisCell >= (dt - timeTravelled)) // if we don't reach the end of this grid cell within the remaining time
            {
                double dx = velocity(j)*(dt - timeTravelled);
                distanceTravelled += dx;
                batch.m_position += dx;

                timeTravelled += (dt - timeTravelled); // equivalent to setting timeTravelled = dt, but += is more clear
            }
            else // if we move outside this grid cell within the remaining time
            {
                double dx = velocity(j)*maxTimeInThisCell;
                batch.m_position += dx;
                distanceTravelled += dx;

                timeTravelled += maxTimeInThisCell;
                j++; // go to next cell

                if (j >= gridPoints.n_elem - 1)
                {
                    // batch have moved outside the grid
                    break; // break loop, this batch will be removed from the array later
                }
            }
        }

        if (i == sword(batches.size()) - 1 && batch.m_position >= gridPoints.tail(1)(0))
        {
            // if at the last batch in the vector, and the position of the batch is outside the grid
            if (batches.size() > 1)
            {
                batches.pop_back();
            }
            else
            {
                // do nothing if we only have one batch left
            }
        }
    }

    if (batches.size() == 0)
    {
        throw std::runtime_error("no batches left, something terrible has happened");
    }

    // extract inlet composition
    const vec inletComposition = inletAndOutletComposition.col(0);

    // check if inlet composition has changed
    const vec diffElem = arma::abs(batches.at(0).concentration() - inletComposition);
    const double diff = arma::sum(arma::abs(diffElem));
    if (diff < 1e-10)
    {
        // if inlet composition hasn't changed, we don't insert a new batch, but
        // reset the position of the first batch
        batches.at(0).m_position = gridPoints(0);
    }
    else
    {
        // if the inlet composition has changed we insert a new batch at the
        // front, using composition from boundary conditions
        // size of batch will be determined by velocity of the (previous)
        // frontmost batch -- meaning that the new batch will take up the space
        // between the first gridpoint and the next batch
        if (batches.at(0).position() > gridPoints(0))
        {
            // only do this if the first batch has actually moved
            // if we have zero inlet velocity this batch might not move during
            // advection, so then we just keep it as is (no gas actually enters
            // the pipeline with zero flow, so this seems reasonable)
            batches.insert(batches.begin(), Batch(0, inletComposition));
        }
    }

    return newState;
}
