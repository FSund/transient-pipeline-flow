#pragma once

#include <vector>
#include <armadillo>

class Pipeline;
class BoundaryConditions;
class BatchTrackingState;

/*!
 * \brief A class for calculating the time development of the gas composition from the gas velocity.
 *
 * This class implements "batch tracking", which is a method for calculating the
 * time development of the gas composition. This is done by setting up "batches"
 * of gas with a given composition, and advancing each batch according to the
 * gas velocity around each batch.
 *
 * The batch tracking procedure is documented in <a href="https://doi.org/10.1016/j.jngse.2018.03.014"><i>Gas composition tracking in transient pipeline flow</i> (Chaczykowski et. al., Journal of Natural Gas Science and Engineering 2018)</a>.
*/

class BatchTracking
{
public:
    //! Declared to avoid the inline compiler-generated default destructor.
    ~BatchTracking();

    /*!
     * \brief A wrapper around BatchTracking::advect(const BatchTrackingState&, const double, const arma::mat&, const arma::vec&)
     * @see BatchTracking::advect(const BatchTrackingState&, const double, const arma::mat&, const arma::vec&)
     * \param dt Time step [s].
     * \param state BatchTrackingState to advect
     * \param pipeline Instance of Pipeline containing the current Pipeline state. Will get velocity from this.
     * \param boundaryConditions Boundary conditions.
     * \return A copy of state with updated batch positions.
     */
    static BatchTrackingState advect(
            const BatchTrackingState& state,
            const arma::uword dt,
            const Pipeline& pipeline,
            const BoundaryConditions& boundaryConditions);

    /*!
     * \brief Calculate new Batch positions from gas velocity.
     *
     * This is the main advection function, which performs the calculations (all
     * other advect functions are just wrappers around this one).
     *
     * Advection is performed by translating each Batch according to the velocity
     * of the gas around each Batch.
     *
     * The batch tracking procedure is documented in <a href="https://doi.org/10.1016/j.jngse.2018.03.014"><i>Gas composition tracking in transient pipeline flow</i> (Chaczykowski et. al., Journal of Natural Gas Science and Engineering 2018)</a>.
     *
     * \param state The current state, the positions and concentrations of all batches.
     * \param dt Time step [s].
     * \param inletAndOutletConcentration Inlet and outlet concentration as column vectors.
     * \param velocity The velocity in the pipeline.
     * \return A copy of state with updated batch positions.
     */
    static BatchTrackingState advect(
            const BatchTrackingState& state,
            const arma::uword dt,
            const arma::mat& inletAndOutletConcentration,
            const arma::vec& velocity);
};
