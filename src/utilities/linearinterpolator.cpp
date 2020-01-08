#include "linearinterpolator.hpp"

#include <iostream>
#include <stdexcept>
#include <cassert>

using arma::uvec;
using arma::vec;

utils::LinearInterpolator::LinearInterpolator(
        const vec& positions,
        const vec& values,
        const arma::uword interpolationMethod):
    m_N(positions.n_elem),
    m_positions(positions),
    m_values(values),
    m_interpolationMethod(interpolationMethod)
{
    // do some checking of size of inputs
    if (values.n_elem != positions.n_elem && values.n_elem != (positions.n_elem - 1))
    {
        // values should either be defined in between the points (N-1 values), or at the points (N values)
        throw std::runtime_error("wrong size of values vector.");
    }

    switch (interpolationMethod)
    {
        case 0:
        {
            // check that we have N-1 values
            if (values.n_elem != (positions.n_elem - 1))
            {
                if (values.n_elem == positions.n_elem)
                {
                    std::cout << "WARNING: LinearInterpolator::LinearInterpolator(): values vector 1 too large for \"no interpolation\", discarding final point" << std::endl;
                    m_values = values.head(values.n_elem - 1);

                    assert(m_values.n_elem == (positions.n_elem - 1));
                }
                else
                {
                    throw std::runtime_error("wrong size of values vector.");
                }
            }
            break;
        }
        case 1:
        {
            // check that we have N values
            if (values.n_elem != positions.n_elem)
            {
                throw std::runtime_error("wrong size of values vector");
            }
            break;
        }
        default:
        {
            throw std::runtime_error("unknown interpolation method (or not implemented yet)");
        }
    }

    m_gradients = arma::zeros<vec>(m_N - 1); // gradients are defined between positions, so we get n-1 gradients

    switch (m_interpolationMethod)
    {
        case 0:
        {
            // no interpolation, will just use value from previous point, until next point
            m_gradients = arma::zeros<vec>(m_N - 1);
            break;
        }
        case 1:
        {
            // linear interpolation
            for (arma::uword i = 0; i < m_N-1; i++)
            {
                double dx = m_positions(i+1) - m_positions(i);
                double dy = m_values(i+1) - m_values(i);
                m_gradients(i) = dy/dx;
            }
            break;
        }
    }
}

vec utils::LinearInterpolator::getValuesAtPoints(
        const vec& referencePoints,
        const vec& referenceValues,
        const vec& points,
        const arma::uword interpolationMethod)
{
    LinearInterpolator interpolator(referencePoints, referenceValues, interpolationMethod);
    return interpolator.getValuesAtPoints(points);
}

vec utils::LinearInterpolator::getValuesAtPoints(const vec& points)
{
    vec values = arma::zeros<vec>(points.n_elem);
    for (arma::uword i = 0; i < points.n_elem; i++)
    {
        values(i) = getValueAtPoint(points(i));
    }

    return values;
}

double utils::LinearInterpolator::getValueAtPoint(const double position)
{
    if (position > m_positions(m_N-1)) // if position is beyond last defined position
    {
        if (m_printOutsideRangeWarning)
        {
            std::cout << "WARNING: LinearInterpolator::getValueAtPoint(): Requested position " << position << " outside defined range. Returning value at LAST point (x = " << m_positions(m_N-1) << ", y = " << std::endl;
        }

        switch (m_interpolationMethod)
        {
            case 0:
            {
                if (m_printOutsideRangeWarning) std::cout << m_values(m_N-2) << ")." << std::endl;
                return m_values(m_N-2);
            }
            case 1:
            {
                if (m_printOutsideRangeWarning) std::cout << m_values(m_N-1) << ")." << std::endl;
                return m_values(m_N-1);
            }
        }
    }
    else if (position < m_positions(0)) // if position is before first defined position
    {
        if (m_printOutsideRangeWarning)
        {
            std::cout << "WARNING: LinearInterpolator::getValueAtPoint(): Requested position " << position << " outside defined range. Returning value at FIRST point (x = " << m_positions(0) << ", y = " << m_values(0) << ")." << std::endl;
        }
        return m_values(0);
    }
    else if (position == m_positions(m_N-1)) // if position is identical to the last defined position
    {
        switch (m_interpolationMethod)
        {
            case 0:
            {
                return m_values(m_N-2);
            }
            case 1:
            {
                return m_values(m_N-1);
            }
        }
    }

    // index of first value before or _at_ position
    // this find() should not give any errors, since we have checked that we are within range (above)
    uvec indexVec = arma::find(m_positions > position, 1, "first"); // subtract one to get point _before_ position
    arma::uword index = indexVec(0) - 1;

    double value;
    if (position == m_positions(index)) // if we are _at_ the point
    {
        // return the value in the point instead of interpolating
        // (have already checked that we are not at the last point, above)
        value = m_values(index);
    }
    else
    {
        double dx = position - m_positions(index); // distance from first point before "position", to "position"
        value = m_values(index) + m_gradients(index)*dx;
        // gradient is 0 if using no interpolation, dy/dx if using linear interpolation
    }

    return value;
}
