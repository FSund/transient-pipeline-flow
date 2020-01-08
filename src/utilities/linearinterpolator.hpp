#pragma once

#include <armadillo>

namespace utils
{
    class LinearInterpolator
    {
        // contains methods for finding interpolated points between given points
        // used for finding height and other properties when setting up the grid

    public:
        LinearInterpolator(
                const arma::vec& positions,
                const arma::vec& values,
                const arma::uword interpolationMethod = 0 // default is no interpolation, only height uses linear interpolation
                );

        double getValueAtPoint(const double position);
        arma::vec getValuesAtPoints(const arma::vec& points);
        static arma::vec getValuesAtPoints(
                const arma::vec& referencePoints,
                const arma::vec& referenceValues,
                const arma::vec& points,
                const arma::uword interpolationMethod = 0);

    private:
        const size_t m_N;
        arma::vec m_positions;
        arma::vec m_values;
        arma::vec m_gradients;
        arma::uword m_interpolationMethod; // 0 means constant, 1 means linear

        bool m_printOutsideRangeWarning = false;
    };
}
