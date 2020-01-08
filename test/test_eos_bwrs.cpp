#include "debug.hpp"
#include "equationofstate/bwrs.hpp"
#include "composition.hpp"

using namespace arma;

TEST_SUITE_BEGIN("BWRS");

TEST_CASE("Compare to JFH implementation")
{
    Composition composition(vec({89.16, 7.3513, 0.5104, 0.0311, 0.0251, 0.0024, 0.0009, 0.1, 0.6980, 1.2208}));
    composition.normalize();

//    BWRS eos = BWRS::stringToSelectParameterSet(composition, "JFH");
    BWRS eos(composition);

    const double pressure = 1e6;
    const double temperature = 273.15;
    const vec Z_cpp = eos.evaluate(pressure, temperature);

    const vec Z_JFH = {
        0.970179340372898, // Z-factor
        0.000352691578677283, // dzdT_p
        -3.00464363541499e-08, // dzdp_T
        0.000235401433608496, // dzdT_rho
        8.31647035050504, // density
        515.994756448608 // cp-cv
    };

    // allow some differences, since we are using different AB parameters,
    // binary interaction coefficients and critical parameters than the JFH code
    CHECK(Z_cpp(0) == doctest::Approx(Z_JFH(0)).epsilon(0.001)); // Z-factor
    CHECK(Z_cpp(1) == doctest::Approx(Z_JFH(1)).epsilon(0.0001)); // dzdT_p
    CHECK(Z_cpp(2) == doctest::Approx(Z_JFH(2))); // dzdp_T
    CHECK(Z_cpp(3) == doctest::Approx(Z_JFH(3))); // dzdT_rho

    const double Z = Z_cpp(0);
    const double R = eos.getGasConstant(); // gas constant
    const double rho_m = pressure/(Z*R*temperature); // molar density [mol / m3] (?)
    const double M = eos.getMolarMassOfMixture();
    const double rho = rho_m*M/1000.0;

    // check that we have actually calculated the density correctly from the molar density
    CHECK(rho == doctest::Approx(eos.findDensity(pressure, temperature)));

    CHECK(rho == doctest::Approx(Z_JFH(4)).epsilon(0.03));
}

TEST_CASE("BWRS and zeroed components")
{
    const double p = 10e6;
    const double T = 273.15;

    for (uword i = 0; i < 10; i++)
    {
        Composition c = Composition::defaultComposition;
        c(i) = 0;
        c.normalize();
        BWRS eos(c);
        CHECK_NOTHROW(eos.findDensity(p, T));
    }
}

TEST_CASE("BWRS vs. SIM with no C1")
{
    const double p = 10e6;
    const double T = 273.15;
    Composition c({0, 0.5, 0.15, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05});
    BWRS eos(c);
    CHECK(eos.findDensity(p, T) == doctest::Approx(555.8).epsilon(0.01));
    const double Z = eos.evaluate(p, T)(0);
    CHECK(Z == doctest::Approx(0.3373).epsilon(0.01));

//BWRS: D 555,8 kg/m³, Z 0,3373
//GERG: D 333,9 kg/m³, Z 0,5616
}

TEST_SUITE_END();
