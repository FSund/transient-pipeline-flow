#include "debug.hpp"

#include "utilities/physics.hpp"
#include "constants.hpp"

using namespace arma;
using namespace std;

TEST_CASE("viscosity")
{
    const double temperature = 280;
    const double molarMass = 17.94565218;
    const double density = 102.0119264377514;
    const double viscosity = 1.426135632811405e-05; // from Matlab code
    CHECK(utils::calculateViscosity({molarMass}, {temperature}, {density})(0) == Approx(viscosity));
}

TEST_CASE("heat capacity cv TGNet")
{
    const double temperature = 280;
    const double molarMass = 17.94565218;
    const double pressure = 10e6;
    CHECK(utils::calculateHeatCapacityConstantVolumeTGNet(molarMass, pressure, temperature) == Approx(1677.83));
}

TEST_CASE("heat capacity cv JFH")
{
    const double p = 10e6;
    const double cv = 1.797600000000000e+03;
    CHECK(utils::calculateHeatCapacityConstantVolumeJFH(p) == Approx(cv));
}

TEST_CASE("heat capacity cp")
{
    const double temperature = 280;
    const double molarMass = 17.94565218;
    const double pressure = 10e6;

    // all these forms use more accurate conversions to/from imperial units than
    // were in the original code, so there are some deviations, so we allow up
    // to 0.1 % difference
    SUBCASE("Langelandsvik")
    {
        const double cp = 3.532583095670939e+03;
        CHECK(utils::calculateHeatCapacityConstantPressureLangelandsvik(molarMass, pressure, temperature) == Approx(cp).epsilon(0.001));
    }
    SUBCASE("JFH/TGNet")
    {
        const double cp = 2.141823171762757e+03;

        CHECK(utils::calculateHeatCapacityConstantPressureJFH(molarMass, pressure, temperature) == Approx(cp).epsilon(0.001));
        CHECK(utils::calculateHeatCapacityConstantPressureTGNet(molarMass, pressure, temperature) == Approx(cp).epsilon(0.001));
    }
}

TEST_CASE("heat capacity cp KIO")
{
    SUBCASE("ideal gas part from Example 1, with imperial units")
    {
        // these are the imperial unit constants
        const double a_1 =  8.0211;
        const double a_2 =  3.3359;
        const double b_1 =  2.0744e-2;
        const double b_2 = -4.2441e-3;
        const double c_1 = -8.1528e-6;
        const double c_2 =  4.8536e-6;
        const double d_1 =  1.2887e-9;
        const double d_2 = -1.1626e-8;

        const double specificGravity = 0.6;
        const double T = 400; // F
        const double cp_ideal =
                + (a_1*specificGravity + a_2)
                + (b_1*specificGravity + b_2)*T
                + (c_1*specificGravity + c_2)*std::pow(T, 2.0)
                + (d_1*specificGravity + d_2)*std::pow(T, 3.0);

        CHECK(cp_ideal == Approx(10.7288));
    }

    SUBCASE("residual part from Example 1")
    {
        const double P = 5000; // [psi]
        const double T = 400 + 460; // [Rankine]
        const double Ppc = 676.862; // critical pressure [psi]
        const double Tpc = 352.26; // critical temperature [Rankine]
        const double reducedPressure = P/Ppc;
        const double reducedTemperature = T/Tpc;
        const double cp_residual = utils::details::KIOdimensionlessResidualCP(reducedPressure, reducedTemperature);
        CHECK(cp_residual == Approx(0.81115)); // from article
        CHECK(cp_residual*constants::gasConstantBTU == Approx(1.6109).epsilon(1e-4)); // from article
    }

    SUBCASE("ideal gas part")
    {
        // extracted approximate values from fig. 6 in article
        CHECK(utils::details::KIOidealGasCP(0.90, 200) == Approx(40).epsilon(0.01));
        CHECK(utils::details::KIOidealGasCP(0.65, 300) == Approx(40).epsilon(0.01));
        CHECK(utils::details::KIOidealGasCP(0.53, 400) == Approx(40).epsilon(0.001));

        CHECK(utils::details::KIOidealGasCP(0.82, 1000) == Approx(100).epsilon(0.01));
        CHECK(utils::details::KIOidealGasCP(0.92, 600) == Approx(80).epsilon(0.01));
        CHECK(utils::details::KIOidealGasCP(0.65, 600) == Approx(60).epsilon(0.01));
    }

//    SUBCASE("residual part")
//    {
//        // from fig. 10
//        // should be around 1.0
//        cout << utils::details::KIOdimensionlessResidualCP(0.5, 1.1) << " close to 1.0?" << endl;

//        // should be around 0.7
//        cout << utils::details::KIOdimensionlessResidualCP(1.0, 1.6) << " close to 0.7?" << endl;

//        // should be around 4
//        cout << utils::details::KIOdimensionlessResidualCP(1.0, 1.1) << " close to 4.0?" << endl;

//        // should be around 0.1
//        cout << utils::details::KIOdimensionlessResidualCP(0.5, 2.4) << " close to 0.1?" << endl;

//        // should be around 8
//        cout << utils::details::KIOdimensionlessResidualCP(2.0, 1.2) << " close to 8.0?" << endl;
//    }
}

#include "composition.hpp"
#include "equationofstate/gerg04.hpp"
TEST_CASE("J-H-EOS cp" * doctest::skip(true))
{
//    const Composition comp = Composition::defaultComposition;

    // from Table 1
    const Composition comp = Composition({71.5, 11.0, 6.5, 0.9, 1.9, 0.4, 0.4, 0, 1.8, 4.0}).normalized();
//    const Composition comp = Composition({80, 11.0, 6.5, 0.9, 1.9, 0.4, 0.4, 0, 1.8, 4.0}).normalized();

    arma::mat results(10, 3);
    const double temperature = 295.9; // from article Fig. 9 a)
    for (uword i = 0; i < 10; i++)
    {
        const double pressure = (i+1)*2e6; // 2 to 20 MPa
        GERG04 gerg(comp);
        const double Z = gerg.evaluate(pressure, temperature)(0);
        std::cout << "JKH: " << utils::calculateIsobaricHeatCapacityJKH(comp, pressure, temperature)
                  << ", JKH + GERG: " << utils::calculateIsobaricHeatCapacityJKH(comp, pressure, temperature, Z)
                  << ", Kareem et al: " << utils::calculateHeatCapacityConstantPressureKIO(comp, pressure, temperature)*(gerg.getMolarMassOfMixture()/1000.0) // convert from [J/kg K] to [J/mol K]
                  << " @ " << pressure/1e5 << " bar"
                  << std::endl;

        results(i, 0) = pressure;
        results(i, 1) = utils::calculateIsobaricHeatCapacityJKH(comp,pressure, temperature, Z);
        results(i, 2) = utils::calculateHeatCapacityConstantPressureKIO(comp, pressure, temperature)*(gerg.getMolarMassOfMixture()/1000.0); // convert from [J/kg K] to [J/mol K]
    }
    cout << results.t() << endl;
}

TEST_CASE("Colebrook-White" * doctest::skip(true))
{}

TEST_CASE("Haaland friction")
{
    const double diameter = 0.9;
    const double Re = 6275904;
    const double ep = 1.7e-6;
    const double friction = 0.008807638512811;
    CHECK(utils::calculateHaalandFrictionFactor(ep, diameter, Re) == Approx(friction));
}
