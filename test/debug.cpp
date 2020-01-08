#include "debug.hpp"

bool equal(const arma::vec& v1, const arma::vec& v2, const double tolerance)
{
    for (std::size_t i = 0; i < v1.size(); i++)
    {
        if (v1.at(i) != Approx(v2.at(i)).epsilon(tolerance))
        {
            return false;
        }
    }

    return true;
}

bool equal(const arma::mat& v1, const arma::mat& v2, const double tolerance)
{
    for (std::size_t i = 0; i < v1.n_rows; i++)
    {
        for (std::size_t j = 0; j < v1.n_cols; j++)
        {
            if (v1.at(i, j) != Approx(v2.at(i, j)).epsilon(tolerance))
            {
                return false;
            }
        }
    }

    return true;
}
