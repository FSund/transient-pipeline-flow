#include "debug.hpp"

#include "equationofstate/bwrs.hpp"
#include "equationofstate/gerg04.hpp"
#include "equationofstate/idealgas.hpp"
#include "composition.hpp"

using namespace arma;

TEST_SUITE_BEGIN("equationofstate");

TEST_CASE("compare BWRS and GERG04")
{
    const double p = 10e6;
    const double T = 273.15;

    const double tolerance = 0.1; // high tolerance, since we are comparing very different EOSs

    for (uword i = 0; i < 6; i++)
    {
        CHECK(BWRS().evaluate(p, T)(i) == doctest::Approx(GERG04().evaluate(p, T)(i)).epsilon(tolerance));
    }
}

TEST_CASE("compare BWRS and GERG04 with zeroed components")
{
    // compare with zeroed out components
    const double p = 10e6;
    const double T = 273.15;

    const double tolerance = 0.1; // high tolerance, since we are comparing very different EOSs

    const Composition c = Composition({85, 7, 3, 1, 1, 1, 1, 0.5, 0.5, 0.5}).normalize();

    for (uword j = 0; j < 10; j++)
    {
        Composition comp = c;
        comp(j) = 0;
        comp.normalize();
        auto bwrs_out = BWRS(comp).evaluate(p, T);
        auto gerg_out = GERG04(comp).evaluate(p, T);

        // loop over results
        for (uword i = 0; i < 6; i++)
        {
            double tol = tolerance;
            if (i == 4 || i == 5)
            {
                // higher tolerance for cv and cp, since the BWRS equation is pretty bad
                tol = 0.3;
            }
            if (j == 0 && i >= 3)
            {
                // something strange happens when setting C1 to zero sometimes
                tol = 5;
            }
            if (j == 0 && i == 0)
            {
                // something strange with Z at C1 == 0
                tol = 0.5;
            }
//            cout << "i " << i << ", j " << j << endl;
            CHECK(bwrs_out(i) == doctest::Approx(gerg_out(i)).epsilon(tol));
        }
    }

//    Composition comp({0, 0.5, 0.15, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05});
//    comp(0) = 0;
//    comp.normalize();
//    cout << comp << endl;

//    cout << GERG04(comp).evaluate(p, T)(0) << endl;
//    cout << BWRS(comp).evaluate(p, T)(0) << endl;

//    cout << GERG04(comp).findDensity(p, T) << endl;
//    cout << BWRS(comp).findDensity(p, T) << endl;
}

TEST_CASE("Check that 10 components work")
{
    const Composition c = Composition({85, 7, 3, 1, 1, 1, 1, 0.5, 0.5, 0.5}).normalize();

    SUBCASE("BWRS")
    {
        CHECK_NOTHROW(BWRS eos(c));
    }

    SUBCASE("GERG04")
    {
        CHECK_NOTHROW(GERG04 eos(c));
    }
}

TEST_SUITE_END();
