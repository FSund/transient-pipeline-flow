#include "debug.hpp"
#include "composition.hpp"

using arma::zeros;
using arma::uword;
using arma::vec;

TEST_CASE("Composition constructors")
{
    CHECK_THROWS(Composition(zeros<vec>(9)));
    CHECK_THROWS(Composition(zeros<vec>(11)));
}

TEST_CASE("Composition normalize")
{
    // fixed
    CHECK(arma::sum(Composition({1, 2, 3, 4, 5, 6, 7, 8, 9, 10}).normalized().vec()) == doctest::Approx(1.0));
    CHECK(Composition({1, 2, 3, 4, 5, 6, 7, 8, 9, 10}).normalized().isNormalized());
    CHECK(Composition({1, 2, 3, 4, 5, 6, 7, 8, 9, 10}).normalize().isNormalized());

    // random composition
    arma::arma_rng::set_seed_random(); // set seed to random value
    for (arma::uword i = 0; i < 10; i ++)
    {
        vec comp = arma::randu<vec>(10);

        // the random function in armadillo doesn't always work, so have to check here
        if (arma::sum(comp) > 0)
        {
            CHECK(arma::sum(Composition(comp).normalized().vec()) == doctest::Approx(1.0));
            CHECK(arma::sum(Composition(comp).normalize().vec()) == doctest::Approx(1.0));
            CHECK(Composition(comp).normalized().isNormalized());
            CHECK(Composition(comp).normalize().isNormalized());
        }
    }
}

TEST_CASE("comparison operator")
{
    CHECK(Composition({1,2,3,4,5,6,7,8,9,10}) == Composition({1,2,3,4,5,6,7,8,9,10}));
    CHECK(Composition::defaultComposition == Composition::defaultComposition);

    CHECK_FALSE(Composition({1.1,2,3,4,5,6,7,8,9,10}) == Composition({1,2,3,4,5,6,7,8,9,10}));
    CHECK_FALSE(Composition({1,2.1,3,4,5,6,7,8,9,10}) == Composition({1,2,3,4,5,6,7,8,9,10}));
    CHECK_FALSE(Composition({1,2,3.1,4,5,6,7,8,9,10}) == Composition({1,2,3,4,5,6,7,8,9,10}));
    CHECK_FALSE(Composition({1,2,3,4.1,5,6,7,8,9,10}) == Composition({1,2,3,4,5,6,7,8,9,10}));
    CHECK_FALSE(Composition({1,2,3,4,5.1,6,7,8,9,10}) == Composition({1,2,3,4,5,6,7,8,9,10}));
    CHECK_FALSE(Composition({1,2,3,4,5,6.1,7,8,9,10}) == Composition({1,2,3,4,5,6,7,8,9,10}));
    CHECK_FALSE(Composition({1,2,3,4,5,6,7.1,8,9,10}) == Composition({1,2,3,4,5,6,7,8,9,10}));
    CHECK_FALSE(Composition({1,2,3,4,5,6,7,8.1,9,10}) == Composition({1,2,3,4,5,6,7,8,9,10}));
    CHECK_FALSE(Composition({1,2,3,4,5,6,7,8,9.1,10}) == Composition({1,2,3,4,5,6,7,8,9,10}));
    CHECK_FALSE(Composition({1,2,3,4,5,6,7,8,9,10.1}) == Composition({1,2,3,4,5,6,7,8,9,10}));
}

TEST_CASE("get by name")
{
    Composition comp({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    CHECK(comp.C1() == 1);
    CHECK(comp.C2() == 2);
    CHECK(comp.C3() == 3);
    CHECK(comp.iC4() == 4);
    CHECK(comp.nC4() == 5);
    CHECK(comp.iC5() == 6);
    CHECK(comp.nC5() == 7);
    CHECK(comp.C6() == 8);
    CHECK(comp.N2() == 9);
    CHECK(comp.CO2() == 10);
}
