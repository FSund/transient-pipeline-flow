#include "solver/governingequationsolver.hpp"

#include <memory>
#include <stdexcept>
#include <armadillo>

#include "solver/discretizer/discretizer.hpp"
#include "solver/discretizer/enthalpy.hpp"
#include "solver/discretizer/internalenergy.hpp"
#include "solver/boundaryconditions.hpp"
#include "solver/matrixequation.hpp"
#include "utilities/utilities.hpp"
#include "pipeline.hpp"

using arma::uword;
using arma::vec;
using arma::uvec;
using arma::zeros;
using arma::mat;
using std::unique_ptr;

GoverningEquationSolverBase::~GoverningEquationSolverBase()
{}

bool GoverningEquationSolverBase::isOverDetermined(const BoundaryConditions& boundaryConditions)
{
    if (boundaryConditions.nActiveBoundaryConditions() > 3)
    {
        return true;
    }
    else
    {
        return false;
    }
}

template<typename T>
GoverningEquationSolver<T>::~GoverningEquationSolver()
{}

template<typename T>
GoverningEquationSolver<T>::GoverningEquationSolver(
        const uword nGridPoints):
    m_matrixEquation(std::make_unique<MatrixEquation>())
{
    m_discretizer = std::make_unique<T>(nGridPoints);
}

template<typename T>
GoverningEquationSolver<T>::GoverningEquationSolver(
        std::unique_ptr<Discretizer> discretizer):
    m_discretizer(std::move(discretizer)),
    m_matrixEquation(std::make_unique<MatrixEquation>())
{}

template<typename T>
mat GoverningEquationSolver<T>::solve(
        const arma::uword dt,
        const Pipeline& currentState,
        const Pipeline& newState,
        const BoundaryConditions& boundaryConditions)
{
    m_discretizer->discretize(dt, currentState, newState);
    m_matrixEquation->fillCoefficientMatrixAndConstantsVector(
                newState.gridPoints().n_elem,
                m_nVariables,
                boundaryConditions,
                m_discretizer->term_i(),
                m_discretizer->term_ipp(),
                m_discretizer->boundaryTerms());

    return m_matrixEquation->solve(
                newState.gridPoints().n_elem,
                m_nVariables,
                boundaryConditions);
}

// explicit instantiantion
template class GoverningEquationSolver<InternalEnergyDiscretizer>;
template class GoverningEquationSolver<EnthalpyDiscretizer>;
