#include "debug.hpp"
#include "equationofstate/gerg04.hpp"
#include "equationofstate/bwrs.hpp"
#include "composition.hpp"
#include <iostream>
#include <iomanip>
#include <string>

using namespace arma;

TEST_CASE("GERG 2004")
{
    SUBCASE("Compare to JFH implementation")
    {
        vec composition = vec({89.16, 7.3513, 0.5104, 0.0311, 0.0251, 0.0024, 0.0009, 1.0, 0.6980, 1.2208})/100;
        GERG04 eos(composition);
        double pressure = 1e6;
        double temperature = 273.15;
        vec Z = eos.evaluateAllProperties(pressure, temperature);

        vec Z_JFH = {
            9.682420267050561e-01,
            3.948866369387216e-04,
           -3.201857378379282e-08,
            2.724539991965982e-04,
           -9.780111514502604e+02,
           -1.652062131164945e+05,
            1.604097510989244e+03,
           -4.548440810264379e+04,
            2.129437090778683e+03,
            2.147539283150379e+05,
            6.338277405377913e-06,
            3.922279907682157e+02,
           -1.349696299865620e-02,
            8.354016295842076e+00,
            1.285002316196897e+00
        };

        for (uword i = 0; i < Z.n_elem; i++)
        {
            // allow 0.1 % error, since we have fixed some bugs in JFH implementation
            CHECK(Z(i) == doctest::Approx(Z_JFH(i)).epsilon(0.001));
        }
    }
}

TEST_CASE("GERG04 and zeroed components")
{
    const double p = 10e6;
    const double T = 273.15;

    const Composition comp = Composition({89.16, 7.3513, 0.5104, 0.0311, 0.0251, 0.0024, 0.0009, 0.05, 0.6980, 2.2208}).normalize();

    for (uword i = 0; i < 10; i++)
    {
        Composition c = comp;
        c(i) = 0;
        c.normalize();
        GERG04 eos(c);
        CHECK_NOTHROW(eos.findDensity(p, T));

        CHECK(eos.indicesOfNonZeroComponents().n_elem == 9);

        // check that the correct index has been detected as having zero fraction
        // this doesn't work the expected way, since GERG04 stores the components
        // in a different order than BWRS and the default way we use
//        cout << i << endl;
//        cout << c.vec().t();
//        cout << eos.indicesOfNonZeroComponents().t() << endl;
//        CHECK_FALSE(arma::any(eos.indicesOfNonZeroComponents() == i));
    }
}


TEST_CASE("GERG vs. SIM with no C1")
{
    const double p = 10e6;
    const double T = 273.15;
    Composition c({0, 0.5, 0.15, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05});
    GERG04 eos(c);
    CHECK(eos.findDensity(p, T) == doctest::Approx(333.9).epsilon(0.03));
    const double Z = eos.evaluate(p, T)(0);
    CHECK(Z == doctest::Approx(0.5616).epsilon(0.07));

//BWRS: D 555,8 kg/m³, Z 0,3373
//GERG: D 333,9 kg/m³, Z 0,5616

//    0.525029
//    0.339636
//    344.174
//    552.089
}
