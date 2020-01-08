#include "solver/discretizer/discretizer.hpp"

using arma::uword;
using arma::cube;
using arma::mat;
using arma::zeros;

Discretizer::~Discretizer()
{}

Discretizer::Discretizer(
        const uword nGridPoints,
        const uword nEquationsAndVariables):
    m_term_i(zeros<cube>(nGridPoints - 1, nEquationsAndVariables, nEquationsAndVariables)),
    m_term_ipp(zeros<cube>(nGridPoints - 1, nEquationsAndVariables, nEquationsAndVariables)),
    m_boundaryTerm(zeros<mat>(nGridPoints - 1, nEquationsAndVariables))
{}
