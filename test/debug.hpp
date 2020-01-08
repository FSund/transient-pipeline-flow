#pragma once
#pragma clang diagnostic ignored "-Woverloaded-shift-op-parentheses"
#include "doctest.h"
using doctest::Approx;
#include <armadillo>

bool equal(const arma::vec& v1, const arma::vec& v2, const double tolerance = 1e-9);

bool equal(const arma::mat& v1, const arma::mat& v2, const double tolerance = 1e-9);
