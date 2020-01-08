#include "debug.hpp"

#include "solver/matrixequation.hpp"
#include "solver/boundaryconditions.hpp"

using arma::uword;
using arma::zeros;
using arma::cube;
using arma::endr;
using arma::mat;

TEST_CASE("MatrixEquation::fillMatrixAndVector")
{
    MatrixEquation s;

    uword nGridPoints = 11;
    uword nEquations = 3;
    uword nVariables = 3;
    uword nEquationsAndVariables = 3;
    cube term_i = zeros<cube>(nGridPoints, nEquations, nVariables);
    cube term_ipp = zeros<cube>(nGridPoints, nEquations, nVariables);
    mat boundaryConditions;
    boundaryConditions
            << 1.0 << 1.0 << endr
            << 1.0 << 1.0 << endr
            << 1.0 << 1.0;
    BoundaryConditions bc(boundaryConditions);
    bc.setBoundarySettings({"inlet", "outlet", "inlet"});

    mat boundaryTerms = zeros<mat>(nGridPoints, nEquations);
    s.fillCoefficientMatrixAndConstantsVector(nGridPoints, nEquationsAndVariables, bc, term_i, term_ipp, boundaryTerms);

    CHECK(all(all(mat(s.coefficients()) == 0)));
    CHECK(arma::all(s.constants() == 0));

    CHECK(s.coefficients().n_cols == 3*nGridPoints - 3);
    CHECK(s.coefficients().n_rows == 3*nGridPoints - 3);
    CHECK(s.constants().n_rows == 3*nGridPoints - 3);
}

// don't think it's my place to test the armadillo implementation
//TEST_CASE("MatrixEquation test Ax = b solver")
//{
//    MatrixEquation s;

//    mat A;
//    vec x, b;

//    A << 1 << 0 << 0 << endr
//      << 0 << 1 << 0 << endr
//      << 0 << 0 << 1;
//    s.setCoefficients(A);
//    b << 1 << 1 << 1;
//    s.setConstants(b);

//    x = arma::spsolve(s.coefficients(), s.constants(), "superlu"); // sparse solver, uses superlu by default
//    CHECK(arma::all(x == vec({1, 1, 1})));

//    A << 3 << 2 << -1 << endr
//      << 2 << -2 << 4 << endr
//      << -1 << 0.5 << -1;
//    s.setCoefficients(A);
//    b << 1 << -2 << 0;
//    s.setConstants(b);
//    x = arma::spsolve(s.coefficients(), s.constants(), "superlu"); // sparse solver, uses superlu by default
//    CHECK(x(0) == doctest::Approx(1));
//    CHECK(x(1) == doctest::Approx(-2));
//    CHECK(x(2) == doctest::Approx(-2));
//}
