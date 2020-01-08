#pragma once

#include <armadillo>

namespace utils
{
    arma::vec tridag(const arma::vec& a, const arma::vec& b, const arma::vec& c,
                     const arma::vec& r, const arma::uword n);
}
