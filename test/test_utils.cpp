#include "debug.hpp"
#include "utilities/utilities.hpp"
#include "utilities/numerics.hpp"

using arma::zeros;
using arma::uword;
using arma::vec;

TEST_SUITE_BEGIN("utils");

TEST_CASE("testing cubeRoot")
{
    vec a = zeros<vec>(1);
    SUBCASE("value")
    {
        CHECK(utils::cubeRoot(a)(0) == 0.0);
        a(0) = 1;
        CHECK(utils::cubeRoot(a)(0) == 1.0);
        a(0) = 8;
        CHECK(utils::cubeRoot(a)(0) == 2.0);
        a(0) = 10;
        CHECK(utils::cubeRoot(a)(0) == doctest::Approx(2.154434690031884));
        a(0) = -0.125;
        CHECK(utils::cubeRoot(a)(0) == doctest::Approx(-0.5));
    }

    SUBCASE("size")
    {
        CHECK(utils::cubeRoot(a).n_elem == 1);
        CHECK(utils::cubeRoot(vec({0,0,0})).n_elem == 3);
    }
}

TEST_CASE("testing utils::smoothTransition")
{
    SUBCASE("0 to 1")
    {
        vec a = zeros<vec>(10);
        vec b = zeros<vec>(10) + 1;
        vec out = utils::smoothTransition(a, b);
        CHECK(out.n_elem == 20);
        CHECK(arma::all(out <= 1.0));
        CHECK(arma::all(out >= 0.0));
        for (uword i = 0; i < out.n_elem - 1; i++)
        {
            CHECK(out(i+1) > out(i));
        }
    }
    SUBCASE("-1 to 1")
    {
        vec a = zeros<vec>(10) - 1;
        vec b = zeros<vec>(10) + 1;
        vec out = utils::smoothTransition(a, b);
        CHECK(out.n_elem == 20);
        CHECK(arma::all(out <= 1.0));
        CHECK(arma::all(out >= -1.0));
        for (uword i = 0; i < out.n_elem - 1; i++)
        {
            CHECK(out(i+1) > out(i));
        }
    }
}

TEST_CASE("numerics")
{
    SUBCASE("arma solve")
    {
        arma::mat A;
        A << 1 << 3 << -2 << arma::endr
          << 3 << 5 << 6 << arma::endr
          << 2 << 4 << 3;

        arma::vec b{5, 7, 8};

        arma::vec x = arma::solve(A, b);
        CHECK(x(0) == doctest::Approx(-15));
        CHECK(x(1) == doctest::Approx(8));
        CHECK(x(2) == doctest::Approx(2));
    }

    SUBCASE("tridag")
    {
        arma::vec a{0, 3, 4};
        arma::vec b{1, 5, 3};
        arma::vec c{3, 6, 0};
        arma::vec r{5, 7, 8};

        arma::vec x1 = utils::tridag(a, b, c, r, a.n_elem);

        arma::mat A;
        A << b(0) << c(0) << 0    << arma::endr
          << a(1) << b(1) << c(1) << arma::endr
          << 0    << a(2) << b(2);

        arma::vec x2 = arma::solve(A, r);

        for (uword i = 0; i < x1.n_elem; i++)
        {
            CHECK(x1(i) == doctest::Approx(x2(i)));
        }
    }
}

TEST_SUITE_END();
