#pragma once

#include <vector>
#include <armadillo>

#include "composition.hpp"

/*! \brief Contains the state which BatchTracking operates on.
 *
 * This contains the position and concentration of all batches in a pipeline.
 * This class also implements methods for sampling the composition at arbitrary
 * locations.
 */
class BatchTrackingState
{
public:
    friend class BatchTracking;

    /*!
     * \brief Contains the information for a single batch
     *
     * This contains the position and concentration of a single batch.
    */
    struct Batch
    {
        /*!
         * \brief Construct Batch from initial position and concentration.
         * \param position The position of the batch.
         * \param concentration The concentration of the batch.
         */
        Batch(const double position, const arma::vec& concentration):
            m_position(position),
            m_concentration(concentration)
        {}

        double m_position; //!< Position of the Batch.
        arma::vec m_concentration; //!< Concentration of the Batch.

    public:
        //! Get the position of the Batch.
        double position() const { return m_position; }
        //! Get the concentration of the Batch
        const arma::vec& concentration() const { return m_concentration; }
    };

    /*!
     * \brief State constructor from grid points and the same concentration for all batches.
     * \param gridPointsIncludingEndPoint Grid points (including endpoint).
     * \param concentration Concentration of all batches.
     * \param nBatches Number of Batches.
     */
    BatchTrackingState(
            const arma::vec& gridPointsIncludingEndPoint,
            const arma::vec& concentration = Composition::defaultComposition,
            const arma::uword nBatches = 1);

    /*!
     * \brief State constructor from grid points and specified Composition at each grid point.
     * \param gridPointsIncludingEndPoint Grid points (including endpoint).
     * \param composition Composition at each grid point.
     */
    BatchTrackingState(
            const arma::vec& gridPointsIncludingEndPoint,
            const std::vector<Composition>& composition);

    /*!
     * \brief sample Samples the composition at locations in m_gridPoints.
     * \return Returns a std::vector of the Composition at each grid point.
     */
    std::vector<Composition> sample() const;

    /*!
     * \brief sample Samples the composition at arbitrary locations.
     * \param locations Where to sample the composition.
     * \return Returns a std::vector of the Composition at each grid point.
     */
    std::vector<Composition> sample(const arma::vec& locations) const;

    /*!
     * \brief sampleToVec Samples the composition at locations in m_gridPoints.
     * \return Returns a std::vector of the composition at each grid point in form of an arma::vec.
     */
    std::vector<arma::vec> sampleToVec() const;

    /*!
     * \brief sampleToVec Samples the composition at arbitrary locations.
     * \param locations Where to sample the composition.
     * \return Returns a std::vector of the composition at each grid point in form of an arma::vec.
     */
    std::vector<arma::vec> sampleToVec(const arma::vec& locations) const;

    /*!
     * \brief sampleInternal Samples the composition at arbitrary locations.
     * This is the main function that does the composition sampling. It simply
     * takes the Composition in the batch each location is contained in.
     * \param locations Where to sample the composition.
     * \return Returns a std::vector of the Composition at each grid point.
     */
    std::vector<arma::vec> sampleInternal(const arma::vec& locations) const;

    /*!
     * \brief Get a const reference to the batches.
     * \return A const reference to the internal std::vector of Batch, m_batches.
     */
    const std::vector<Batch>& batches() const { return m_batches; }

protected:
    std::vector<Batch> m_batches; //!< Vector containing all batches.
    arma::vec m_gridPoints; //!< Grid points of pipeline.
};
