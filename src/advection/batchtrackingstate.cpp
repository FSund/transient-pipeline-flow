#include "advection/batchtrackingstate.hpp"

#include "utilities/stringbuilder.hpp"

using arma::mat;
using arma::vec;
using arma::sword;
using arma::uvec;
using arma::uword;
using arma::zeros;
using std::vector;
using std::cout;
using std::endl;

BatchTrackingState::BatchTrackingState(
        const vec& gridPointsIncludingEndPoint,
        const vec& concentration,
        const uword nBatches):
    m_gridPoints(gridPointsIncludingEndPoint)
{
    if (nBatches == 0)
    {
        const uword nBatches = gridPointsIncludingEndPoint.n_elem - 1;
        m_batches.clear();
        for (uword i = 0; i < nBatches; i++)
        {
            auto position = gridPointsIncludingEndPoint(i);
            m_batches.push_back(Batch(position, concentration));
        }
    }
    else
    {
        m_batches.clear();
        double dx = (gridPointsIncludingEndPoint.tail(1)(0) - gridPointsIncludingEndPoint(0))/double(nBatches); // nBatches+1??
        double position = gridPointsIncludingEndPoint(0);
        for (uword i = 0; i < nBatches; i++)
        {
            m_batches.push_back(Batch(position, concentration));
            position += dx;
        }
    }
}

BatchTrackingState::BatchTrackingState(
        const arma::vec& gridPointsIncludingEndPoint,
        const vector<Composition>& composition):
    m_gridPoints(gridPointsIncludingEndPoint)
{
    if (gridPointsIncludingEndPoint.n_elem != composition.size())
        throw std::runtime_error("incompatible size");

    m_batches.clear();
    for (uword i = 0; i < composition.size() - 1; i++) // skip outlet point
    {
        Batch batch(gridPointsIncludingEndPoint(i), (composition.at(i).vec() + composition.at(i+1).vec())/2.0);
        m_batches.push_back(batch);
    }
}

vector<Composition> BatchTrackingState::sample() const
{
    return sample(m_gridPoints);
}

vector<arma::vec> BatchTrackingState::sampleToVec() const
{
    return sampleInternal(m_gridPoints);
}

vector<arma::vec> BatchTrackingState::sampleToVec(const arma::vec& gridPoints) const
{
    return sampleInternal(gridPoints);
}

vector<vec> BatchTrackingState::sampleInternal(const vec& gridPoints) const
{
    // this samples the concentration at the positions gridPoints, by just
    // taking the composition in the batch the grid point is located

    if (arma::any(gridPoints < m_gridPoints(0)) || arma::any(gridPoints > m_gridPoints.tail(1)(0)))
    {
        throw std::out_of_range("requested sample points not within defined range");
    }

    vec batchPositions = zeros<vec>(m_batches.size());
    for (uword i = 0; i < m_batches.size(); i++)
    {
        batchPositions(i) = m_batches.at(i).m_position;
    }

    // sample composition at grid points
    vector<vec> composition;
    composition.reserve(gridPoints.n_elem);
    for (uword i = 0; i < gridPoints.n_elem; i++)
    {
        // find index
        uvec tmp = arma::find(batchPositions <= gridPoints(i), 1, "last"); // last element of batchPositions that is <= gridPoints(i)
        if (tmp.n_elem == 0)
            throw std::runtime_error(
                    utils::stringbuilder()
                    << "no elements found (batchPositions(end) = " << batchPositions.tail(1)(0)
                    << ", gridPoints(i) = " << gridPoints(i) << ")"
                );
        uword j = tmp(0);
        composition.push_back(m_batches.at(j).m_concentration);
    }

    return composition;
}

vector<Composition> BatchTrackingState::sample(const arma::vec& gridPoints) const
{
    vector<vec> compositionVecs = sampleInternal(gridPoints);

    vector<Composition> composition;
    composition.reserve(gridPoints.n_elem);
    for (uword i = 0; i < gridPoints.n_elem; i++)
    {
        composition.push_back(Composition(compositionVecs.at(i)));
    }

    return composition;

//    // this samples the concentration at the positions gridPoints, by just
//    // taking the composition in the batch the grid point is located

//    vec batchPositions = zeros<vec>(m_batches.size());
//    for (uword i = 0; i < m_batches.size(); i++)
//    {
//        batchPositions(i) = m_batches.at(i).m_position;
//    }

//    // sample composition at grid points
//    vector<Composition> composition;
//    composition.reserve(gridPoints.n_elem);
//    for (uword i = 0; i < gridPoints.n_elem; i++)
//    {
//        // find index
//        uvec tmp = arma::find(batchPositions <= gridPoints(i), 1, "last"); // last element of batchPositions that is <= gridPoints(i)
//        assert(tmp.n_elem > 0);
//        uword j = tmp(0);
//        composition.push_back(Composition(m_batches.at(j).m_concentration));
//    }

//    return composition;
}
