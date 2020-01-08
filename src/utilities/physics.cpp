#include "utilities/physics.hpp"

#include <cmath>

#include "constants.hpp"
#include "composition.hpp"
#include "utilities/utilities.hpp"

using arma::vec;
using arma::uword;
using arma::zeros;

double utils::calculateHeatCapacityConstantVolumeJFH(const double pressure)
{
    // specific heat capacity [J/kg K]
    const double heatCapacityConstantVolume = (-3.1*1e-16*std::pow(pressure, 2) + 1.46*1e-8*pressure + 1.6826)*1000.0;

    return heatCapacityConstantVolume;
}

double utils::calculateHeatCapacityConstantVolumeTGNet(
        const double molarMass,
        const double pressure,
        const double temperature)
{
    // have to use this heat capacity to get the best results
    const double cp = utils::calculateHeatCapacityConstantPressureTGNet(molarMass, pressure, temperature);

    // specific heat capacity [J/kg K]
    // from TGNet manual page 11
    // valid for specific gravity 0.55 to 0.80, temperature 0 F to 150 F, pressure 15 psia to 1400 psia
    const double specificGravity = molarMass/constants::molarMassOfAir;
    const double cv = cp - (1720.0/specificGravity)*constants::footPoundForcePerSlugRankine;

    return cv;
}

double utils::calculateHeatCapacityConstantPressureJFH(
        const double molarMass,
        const double pressure,
        const double temperature)
{
    return utils::calculateHeatCapacityConstantPressureTGNet(molarMass, pressure, temperature);
}

double utils::calculateHeatCapacityConstantPressureLangelandsvik(
        const double molarMass,
        const double pressure,
        const double temperature)
{
    // specific heat capacity [J/kg K]

    // heat capacity constant pressure - isobaric heat capacity
    // relation from Langelandsvik, eq. 2.26 and 2.27
    // not accurate for high pressures
    // (valid for 0-100 barg)
    const double specificGravity = molarMass/constants::molarMassOfAir;
    const double EXP = (
                15.69*1e-2
                *std::pow(pressure*constants::pascalToPoundForcePerSquareFoot, 1.106)
                *std::exp(-6.203*temperature*constants::kelvinToRankine*1e-3)
            )/specificGravity; // Langelandsvik eq. 2.27
    const double heatCapacityConstantPressure = (
                + 1.432*1e4
                - 1.045*1e4*specificGravity
                + 3.255*temperature*constants::kelvinToRankine
                + 10.01*specificGravity*temperature*constants::kelvinToRankine
                + EXP)*constants::footPoundForcePerSlugRankine;

    return heatCapacityConstantPressure;
}

double utils::calculateHeatCapacityConstantPressureTGNet(
        const double molarMassOfMixture,
        const double pressure,
        const double temperature)
{
    // specific heat capacity [J/kg K]
    // from TGNet manual page 11
    // valid for specific gravity 0.55 to 0.80, temperature 0 F to 150 F, pressure 15 psia to 1400 psia
    const double specificGravity = molarMassOfMixture/constants::molarMassOfAir;

    // expression from TGNet (note exp difference from Langelandsvik expression)
    const double EXP = (
                15.69*1e-2
                *std::pow(pressure*constants::pascalToPoundForcePerSquareFoot, 1.106)
                *std::exp(-6.203)*1e-3*temperature*constants::kelvinToRankine
            )/specificGravity;
    const double heatCapacityConstantPressure = (
                + 1.432*1e4
                - 1.045*1e4*specificGravity
                + 3.255*temperature*constants::kelvinToRankine
                + 10.01*specificGravity*temperature*constants::kelvinToRankine
                + EXP)*constants::footPoundForcePerSlugRankine;

    return heatCapacityConstantPressure;
}

double utils::calculateHeatCapacityConstantPressureKIO(
        const Composition& comp,
        const double pressure,
        const double temperature)
{
    // critical pressure [Pa]
    arma::vec Pc = arma::vec({ 45.96,      48.839,     42.5,       36.48,      37.96,      33.81,      33.69,      27.34,      33.99,      73.825})*1e5;
    // critical temperature [K]
    arma::vec Tc = arma::vec({ 190.69,     305.39,     369.89,     408.13,     425.19,     460.37,     469.49,     507.29,     126.15,     304.15});
    // Molar mass of the different gas components [g/mol]
    arma::vec molarMassOfComponents = arma::vec({ 16.042,     30.068,     44.094,     58.12,      58.12,      72.146,     72.146,     86.172,     28.016,     44.01});

    const double molarMass = sum(comp.vec() % molarMassOfComponents);
    const double criticalPressureOfMixture = sum(comp.vec() % Pc);
    const double criticalTemperatureOfMixture = sum(comp.vec() % Tc);

    // correlation from Kareem, Iwalewa, Omeke (2014)
    // valid for pseudo-reduced temperatures in the range [1.1, 3],
    // pseudo-reduced pressures in the range [0.01, 15] and specific gravity
    // in the range of [0.55, 1.0]

    const double specificGravity = molarMass/constants::molarMassOfAir;
    const double cp_ideal = details::KIOidealGasCP(specificGravity, temperature); // [J/mol K]

    const double reducedPressure = pressure/criticalPressureOfMixture;
    const double reducedTemperature = temperature/criticalTemperatureOfMixture;
    const double cp_residual = details::KIOdimensionlessResidualCP(reducedPressure, reducedTemperature); // dimensionless

    const double heatCapacityConstantPressure =
            (
                cp_ideal
                + constants::gasConstant*cp_residual // multiply with R to get dimension [J/mol K]
             )
            /(molarMass/1000.0); // convert from [J/mol K] (molar heat capacity) to [J/kg K]

    return heatCapacityConstantPressure;
}

vec utils::calculateViscosity(
        const vec& molarMass,
        const vec& temperature,
        const vec& density)
{
    const uword nGridPoints = uword(temperature.n_elem);
    const vec temperatureRankine = temperature*constants::kelvinToRankine;

    // dynamic viscosity, using Lee-Gonzales-Eakin correlation, eq. 2.27 in Helgaker
    const vec K = (9.4 + 0.02*molarMass) % arma::pow(temperatureRankine, 1.5)/(209.0 + 19.0*molarMass + temperatureRankine);
    const vec X = 3.5 + 986.0/(temperatureRankine) + 0.01*molarMass;
    const vec Y = 2.4 - 0.2*X;
    vec YY = zeros<vec>(nGridPoints);
    for (uword i = 0; i < nGridPoints; i++)
    {
        YY(i) = std::pow(density(i)/1000.0, Y(i));
    }
    vec viscosity = K % arma::exp(X % YY)/1.0e7; // divide by 1e7 in the end to convert from micropoise [1e-7 kg/(m*s)] to kg/(m*s)

    return viscosity;
}

vec utils::calculateReynoldsNumber(
        const vec& massFlow,
        const vec& diameter,
        const vec& viscosity)
{
    return arma::abs(massFlow)*4.0/(constants::pi*diameter % viscosity);
}

vec utils::calculateColebrookWhiteFrictionFactor(
        const vec& sandGrainEquivalentRoughness,
        const vec& diameter,
        const vec& reynoldsNumber)
{
    vec f = zeros<vec>(diameter.n_elem);
    for (uword i = 0; i < f.n_elem; i++)
    {
        f(i) = calculateColebrookWhiteFrictionFactor(sandGrainEquivalentRoughness(i), diameter(i), reynoldsNumber(i));
    }

    return f;
}

double utils::calculateColebrookWhiteFrictionFactor(
        const double sandGrainEquivalentRoughness,
        const double diameter,
        const double reynoldsNumber)
{
    return details::colebrookWhiteFrictionFactor(sandGrainEquivalentRoughness, diameter, reynoldsNumber);
}

double utils::details::colebrookWhiteFrictionFactor(
        const double sandGrainEquivalentRoughness,
        const double diameter,
        const double reynoldsNumber)
{
    // Friction factor with Coolebrook-White correlation

    if (reynoldsNumber < 0)
    {
        throw std::runtime_error("negative Reynolds number");
    }
    else if (reynoldsNumber < 1e-10)
    {
        return 0; // or inf..?
    }
    else if (reynoldsNumber < 4000)
    {
        return double(64.0)/reynoldsNumber; // valid for Re < 2320
    }

    if (sandGrainEquivalentRoughness < 1e-12)
    {
        return 0;
    }

    // Use Newtons method, guess a starting value of 0.01
    double fp = 1.0;
    double f = 0.01;
    double delta = 1;
    uword counter = 0;
    uword counterLimit = 100;
    double eps = 0.001;
    while (std::abs(delta/fp) > eps) // until relative change converges
    {
        fp = f;
        double eval = -1.0/std::sqrt(fp) + colebrookWhite(fp, sandGrainEquivalentRoughness, diameter, reynoldsNumber);
        double diff = 1.0/(2.0*fp*std::sqrt(fp)) + colebrookWhiteDerivative(fp, sandGrainEquivalentRoughness, diameter, reynoldsNumber);
        if (diff == 0.0)
        {
            // do something (will divide by zero below)
            // TODO: consider throwing no_convergence_error instead
            break; // break out of while()
        }
        delta = eval/diff;
        f = fp - delta;
        if (f < 0)
        {
            f = fp/2.0; // avoid negative friction factor
        }
        if (f <= 0)
            throw std::runtime_error("f <= 0");
        if (counter >= counterLimit)
        {
            // TODO: consider throwing no_convergence_error instead
            std::cout << "Colebrook-White using more than " << counterLimit << " iterations. Breaking loop!" << std::endl;
            break; // break out of while()
        }
        counter++;
    }

    return f;
}

double utils::details::colebrookWhite(
        const double f,
        const double sandGrainEquivalentRoughness,
        const double diameter,
        const double reynoldsNumber)
{
    double c1 = 3.7;
    double c2 = 2.51;
    double c3 = 2.0;

    double value = -c3*std::log10(sandGrainEquivalentRoughness/(c1*diameter) + c2/(reynoldsNumber*std::sqrt(f)));
    return value; // -inf if Re is 0
}

double utils::details::colebrookWhiteDerivative(
        const double f,
        const double sandGrainEquivalentRoughness,
        const double diameter,
        const double reynoldsNumber)
{
    double c1 = 3.7;
    double c2 = 2.51;
    double c3 = 2.0;

    double bracket = sandGrainEquivalentRoughness/(c1*diameter) + c2/(reynoldsNumber*std::sqrt(f)); // inf if Re is 0
    double value = c2*c3/(reynoldsNumber*f*std::sqrt(f)*std::log(10)*bracket*2.0); // -nan if Re is 0
    return value;
}

double utils::calculateHaalandFrictionFactor(
        const double sandGrainEquivalentRoughness,
        const double diameter,
        const double reynoldsNumber)
{
    // Haalands formula
    return 1.0/std::pow(-1.8*std::log10(
                    std::pow(sandGrainEquivalentRoughness/(3.7*diameter), 1.11)
                    + 6.9/reynoldsNumber)
                , 2);
}

double utils::details::KIOidealGasCP(
        const double specificGravity,
        const double temperature)
{
    // Eq. (19) has a correlation of regression of 0.9999 and a maximum relative
    // error of 0.04. It is applicable to samples with specific gravity in the
    // range of [0.55, 1.0] and temperature in the range of [100, 1500] K.

    const double a_1 = -1.09602e1; // kJ/kmol K -- see nomenclature
    const double a_2 =  2.59033e1; // kJ/kmol K
    const double b_1 =  2.1517e-1; // kJ/kmol K^2
    const double b_2 = -6.8687e-2; // kJ/kmol K^2
    const double c_1 = -1.3337e-4; // kJ/kmol K^3
    const double c_2 =  8.6387e-5; // kJ/kmol K^3
    const double d_1 =  3.1474e-8; // kJ/kmol K^4
    const double d_2 = -2.8396e-8; // kJ/kmol K^4

    const double T = temperature; // think Kelvin is okay for this
    const double cp_ideal =
            + (a_1*specificGravity + a_2)
            + (b_1*specificGravity + b_2)*T
            + (c_1*specificGravity + c_2)*std::pow(T, 2.0)
            + (d_1*specificGravity + d_2)*std::pow(T, 3.0);

    return cp_ideal; // [kJ/kmol K]
}

double utils::details::KIOdimensionlessResidualCP(
        const double reducedPressure,
        const double reducedTemperature)
{
    // Eq. (27) has a regression coefficient of 0.997 and maximum relative error
    // of 0.05 compared to the Starling Carnahan equation of state. It is
    // applicable when the pseudo reduced pressure and temperature fall within
    // [0.01, 15] and [1.1, 3], respectively.

    // NB: According to the authors fig. 10 was created by solving the integral
    // in eq. (23), not from eq. (27)

    const double a_1 =  4.80828;
    const double a_2 = -4.01563;
    const double a_3 = -0.0700681;
    const double a_4 =  0.0567;
    const double a_5 =  2.36642;
    const double a_6 = -3.82421;
    const double a_7 =  7.71784;

    const double t = 1.0/reducedTemperature; // inverse reduced temperature
    const double Pprt = reducedPressure*t;
    const double LHS =
            (
                1.0 + utils::pow2(a_1*std::exp(a_2*utils::pow2(1.0 - t))*Pprt)
            )/(
                a_7 + a_6*Pprt + a_5*Pprt*Pprt + a_4*Pprt*Pprt*Pprt
            );
    const double RHS =
            (
                utils::pow2(a_1*std::exp(a_2*utils::pow2(1.0 - t))*Pprt)
                *a_3*utils::pow6(Pprt)
            )/(
                utils::pow3(a_7 + a_6*Pprt + a_5*Pprt*Pprt + a_4*Pprt*Pprt*Pprt)
            );
    const double cp_residual = LHS - RHS;

    // dimensionless -- multiply with R to get appropriate unit
    return cp_residual;
}

// local function
double findZvalueJKH(
        const double A,
        const double B)
{
    if ((false))
    {
        // use Cardano's method
        const double a3 = 1;
        const double a2 = -1;
        const double a1 = (A - B - B*B);
        const double a0 = -A*B;
        const double b0 = a0/a3;
        const double b1 = a1/a3;
        const double b2 = a2/a3;
        const double p = b1 - b2*b2/3.0;
        const double q = b0 - b1*b2/3.0 + 2.0*b2*b2*b2/27.0;
        const double d = std::pow(p/3.0, 3.0) + std::pow(q/2.0, 2.0);

        double Z;
        if (d > 0)
        {
            // only one real root
            const double a = -q/2.0 + std::sqrt(d);
            const double b = -q/2.0 - std::sqrt(d);
            Z = std::pow(a, 1.0/3.0) + std::pow(b, 1.0/3.0) - b2/3.0;
        }
        else
        {
            // three real solutions
            throw std::runtime_error("not implemented");
    //        const double phi = std::acos(-q/2.0*std::sqrt(-27.0/(p*p*p)));
    //        const double x0 = 2.0*std::sqrt(-p/3.0)*std::cos((phi + 2.0*constants::pi*0)/3.0) - b2/3.0;
    //        const double x1 = 2.0*std::sqrt(-p/3.0)*std::cos((phi + 2.0*constants::pi*1)/3.0) - b2/3.0;
    //        const double x2 = 2.0*std::sqrt(-p/3.0)*std::cos((phi + 2.0*constants::pi*2)/3.0) - b2/3.0;
    //        const double x3 = 2.0*std::sqrt(-p/3.0)*std::cos((phi + 2.0*constants::pi*3)/3.0) - b2/3.0;
        }

        return Z;
    }
    else
    {
        // just do 50 steps of Newton's method
        double Z = 1; // starting point
        for (arma::uword i = 0; i < 50; i++)
        {
            const double f = Z*Z*Z - Z*Z + (A - B - B*B)*Z - A*B;
            const double fdiv = 3.0*Z*Z - 2.0*Z + (A - B - B*B);
            Z = Z - f/fdiv;
        }

        return Z;
    }
}

double utils::calculateIsobaricHeatCapacityJKH(
        const Composition& comp,
        const double pressure,
        const double temperature,
        const double compressibility)
{
    // Jarrahian et al
    // "On the isobaric specific heat capacity of natural gas"
    // 2014, Fluid Phase Equilibria 384
    // Evaluates equation 12
    // Falid for a wide range of pressures (0.1–40 MPa) and temperatures (250–414 K)

    const arma::vec molarMass_i = arma::vec({ 16.042,     30.068,     44.094,     58.12,      58.12,      72.146,     72.146,     86.172,     28.016,     44.01});            // [g/mol]
    const double molarMass = sum(comp.vec() % molarMass_i);
    const double specificGravity = molarMass/constants::molarMassOfAir; // [-]

    const double H2S = 0;
    const double H2 = 0;
    const double H2O = 0;
    const double Cp_molar = details::JKHdimensionlessCP(
                specificGravity, H2S, comp.CO2(), comp.N2(), H2, H2O,
                pressure, temperature, compressibility); // [J/mol K]

    const double Cp = Cp_molar/(molarMass/1000.0); // [J/kg K]

    return Cp; // [J/kg K]
}

double utils::details::JKHdimensionlessCP(
        const double specificGravity,
        const double H2S,
        const double CO2,
        const double N2,
        const double H2,
        const double H2O,
        const double pressure,
        const double temperature,
        const double compressibility)
{
    // Jarrahian et al
    // "On the isobaric specific heat capacity of natural gas"
    // 2014, Fluid Phase Equilibria 384
    // Evaluates equation 12
    // Falid for a wide range of pressures (0.1–40 MPa) and temperatures (250–414 K)

    const arma::vec Ji {
        1.19253457299316E-01,
        -2.87407186474540E-01,
        -4.89941222485594E-01,
        -2.36455629352375E-01,
        1.55395604093743E+00,
        -1.38856760022401E-01,
        7.30253696494153E-01,
        -1.18427135730060E-01
    };

    const arma::vec Ki {
        3.75487759217361E+00,
        -3.40472432080060E+00,
        -9.77003391871320E+00,
        -9.47072440761246E+00,
        1.35860847471836E+01,
        -8.94430462055307E-01,
        1.96772476618319E+01,
        -2.99180329435285E+00
    };

    const arma::vec betai {
        2.38242747862715E-01,
        -3.51550147947942E-02,
        6.20467284042863E-01,
        -5.74517899428874E-03,
        -1.18383359572768E-01,
        8.18368533389717E-02
    };

    const arma::vec Ai {
        4.59471825354044E+01,
        9.90750496843086E+00,
        4.17935179794448E-01,
        7.09501951412871E-01,
        -9.02465547872749E+00
    };

    const arma::vec Bi {
        5.57638260250257E-01,
        6.34844709395108E-01,
        -2.68227041459472E-02,
        7.83864423900529E-02,
        -1.11717924190626E-03
    };

    // imperial
//    const double T = temperature*constants::kelvinToRankine; // convert from Kelvin to Rankine
//    const double P = pressure*0.000145037737730; // convert from Pascal to psi
//    const double R = 10.7316; // [ft3 psia/lb-mol °R]
//    const double R = 4.61915; // J/(mol Rankine)

    // SI
    const double R = constants::gasConstant; // [J/mol K] == [m3 Pa /mol K]
    const double T = temperature;
    const double P = pressure;

    const double J = Ji(0) + Ji(1)*H2S + Ji(2)*CO2 + Ji(3)*N2 + Ji(4)*H2 + Ji(4)*H2O
            + Ji(6)*specificGravity + Ji(7)*specificGravity*specificGravity;
    const double K = Ki(0) + Ki(1)*H2S + Ki(2)*CO2 + Ki(3)*N2 + Ki(4)*H2 + Ki(5)*H2O
            + Ki(6)*specificGravity + Ki(7)*specificGravity*specificGravity;

    // pseudocritical properties
    // These pseudocritical properties do not in any way reflect the
    // true critical pressures and critical temperatures of the gas mixture;
    // they are simply the parameters used in the proposed EOS. It is
    // important to note that with regard to oil-industry standards,
    // they return values for psia and degrees R (Jarrahian and Heydaryan, 2014)
    const double TPC_imperial = K*K/J; // psia
    const double PPC_imperial = K*K/(J*J); // Rankine
    const double psiToPascal = 6894.75729;
    const double rankineToKelvin = 5.0/9.0;
    const double TPC = TPC_imperial*rankineToKelvin;
    const double PPC = PPC_imperial*psiToPascal;
//    const double PPC = PPC_imperial*psiToPascal*0.6; // this gives almost reasonable results...

    // pseudoreduced temperature and pressure
    const double Tpr = T/TPC;
    const double Ppr = P/PPC;

    const double beta = betai(0) + betai(1)*std::log(Ppr) + betai(2)/Tpr
            + betai(3)*std::log(Ppr)*std::log(Ppr) + betai(4)/(Tpr*Tpr)
            + betai(5)*std::log(Ppr)/Tpr; // beta_5 or beta_6 ??? beta_6 if we look at eq. (18)

    const double A = 0.49694*beta*Ppr/(Tpr*Tpr);
    const double B = 0.09012*Ppr/Tpr;

    const double Z = compressibility == 0 ? findZvalueJKH(A, B) : compressibility;

    const double a = 0.49694*std::pow(R*TPC, 2.0)/PPC;
    const double b = 0.09012*R*TPC/PPC;

    const double dbetadT = - betai(2)*TPC/(T*T) - 2.0*betai(4)*TPC/(T*T*T)
            - betai(5)*std::log(Ppr)*TPC/(T*T);

    const double d2betadT2 = 2.0*betai(2)*TPC/(T*T*T) + 6.0*betai(4)*TPC/(T*T*T*T)
            + 2.0*betai(5)*std::log(Ppr)*TPC/(T*T*T);

    const double M = Z*(Z + B)/(Z - B);
    const double N = a*B/(R*b)*dbetadT;

    const double C0P = Ai(0) + Ai(1)*Tpr + Ai(2)*Tpr*Tpr + Ai(3)/specificGravity
            + Ai(4)/(specificGravity*specificGravity);

    const double CcorrDL = Bi(0) + Bi(1)*std::log(Tpr) + Bi(2)*Ppr
            + Bi(3)*std::pow(std::log(Tpr), 2.0) + Bi(4)*Ppr*std::log(Tpr);

    const double CpRes =
            a*T/b*d2betadT2*std::log((Z + B)/Z)
            + (R*std::pow(M - N, 2.0))/(M*M - A*(2.0*Z + B))
            - R;

    const double Cp = CcorrDL*CpRes + C0P; // [J/mol K]

    return Cp; // [J/mol K]
}


double utils::details::JKHidealGasCP(
        const Composition& comp,
        const double specificGravity,
        const double temperature)
{
    const double H2S = 0;
    const double H2 = 0;
    const double H2O = 0;
    return utils::details::JKHidealGasCP(
                specificGravity, H2S, comp.CO2(), comp.N2(), H2, H2O, temperature);
}

double utils::details::JKHidealGasCP(
        const double specificGravity,
        const double H2S,
        const double CO2,
        const double N2,
        const double H2,
        const double H2O,
        const double temperature)
{
    // Jarrahian et al
    // "On the isobaric specific heat capacity of natural gas"
    // 2014, Fluid Phase Equilibria 384
    // Evaluates equation 12
    // Falid for a wide range of pressures (0.1–40 MPa) and temperatures (250–414 K)

    const arma::vec Ji {
        1.19253457299316E-01,
        -2.87407186474540E-01,
        -4.89941222485594E-01,
        -2.36455629352375E-01,
        1.55395604093743E+00,
        -1.38856760022401E-01,
        7.30253696494153E-01,
        -1.18427135730060E-01
    };

    const arma::vec Ki {
        3.75487759217361E+00,
        -3.40472432080060E+00,
        -9.77003391871320E+00,
        -9.47072440761246E+00,
        1.35860847471836E+01,
        -8.94430462055307E-01,
        1.96772476618319E+01,
        -2.99180329435285E+00
    };

    const arma::vec Ai {
        4.59471825354044E+01,
        9.90750496843086E+00,
        4.17935179794448E-01,
        7.09501951412871E-01,
        -9.02465547872749E+00
    };

    // SI
    const double T = temperature;

    const double J = Ji(0) + Ji(1)*H2S + Ji(2)*CO2 + Ji(3)*N2 + Ji(4)*H2 + Ji(4)*H2O
            + Ji(6)*specificGravity + Ji(7)*specificGravity*specificGravity;
    const double K = Ki(0) + Ki(1)*H2S + Ki(2)*CO2 + Ki(3)*N2 + Ki(4)*H2 + Ki(5)*H2O
            + Ki(6)*specificGravity + Ki(7)*specificGravity*specificGravity;

    // pseudocritical properties
    // These pseudocritical properties do not in any way reflect the
    // true critical pressures and critical temperatures of the gas mixture;
    // they are simply the parameters used in the proposed EOS. It is
    // important to note that with regard to oil-industry standards,
    // they return values for psia and degrees R (Jarrahian and Heydaryan, 2014)
    const double TPC_imperial = K*K/J; // psia
    const double rankineToKelvin = 5.0/9.0;
    const double TPC = TPC_imperial*rankineToKelvin;

    // pseudoreduced temperature
    const double Tpr = T/TPC;

    const double C0P = Ai(0) + Ai(1)*Tpr + Ai(2)*Tpr*Tpr + Ai(3)/specificGravity
            + Ai(4)/(specificGravity*specificGravity);

    return C0P; // [J/mol K]
}

double utils::details::calculateHeatCapacityConstantPressureKIO(
        const Composition& composition,
        const double H2S,
        const double pressure,
        const double temperature)
{
    //                         C1          C2          C3          iC4         nC4         iC5         nC5         C6          N2          CO2         H2S
    // critical pressure [Pa]
    arma::vec Pc = arma::vec({ 45.96,      48.839,     42.5,       36.48,      37.96,      33.81,      33.69,      27.34,      33.99,      73.825,     89.7})*1e5;
    // critical temperature [K]
    arma::vec Tc = arma::vec({ 190.69,     305.39,     369.89,     408.13,     425.19,     460.37,     469.49,     507.29,     126.15,     304.15,     373.54});
    // Molar mass of the different gas components [g/mol]
    //                                            C1          C2          C3          iC4         nC4         iC5         nC5         C6          N2          CO2         H2S
    arma::vec molarMassOfComponents = arma::vec({ 16.042,     30.068,     44.094,     58.12,      58.12,      72.146,     72.146,     86.172,     28.016,     44.01,      34.081});

    if (H2S > 1)
    {
        throw std::runtime_error("H2S fraction > 1");
    }

    // make new comp with H2S inserted at the end
    arma::vec comp = arma::join_cols(composition.vec(), arma::vec({H2S}));
    comp /= arma::sum(comp); // normalize

    const double molarMass = sum(comp % molarMassOfComponents);
    const double criticalPressureOfMixture = sum(comp % Pc);
    const double criticalTemperatureOfMixture = sum(comp % Tc);

    // correlation from Kareem, Iwalewa, Omeke (2014)
    // valid for pseudo-reduced temperatures in the range [1.1, 3],
    // pseudo-reduced pressures in the range [0.01, 15] and specific gravity
    // in the range of [0.55, 1.0]

    const double specificGravity = molarMass/constants::molarMassOfAir;
    const double cp_ideal = details::KIOidealGasCP(specificGravity, temperature); // [J/mol K]

    const double reducedPressure = pressure/criticalPressureOfMixture;
    const double reducedTemperature = temperature/criticalTemperatureOfMixture;
    const double cp_residual = details::KIOdimensionlessResidualCP(reducedPressure, reducedTemperature); // dimensionless

    const double heatCapacityConstantPressure =
            (
                cp_ideal
                + constants::gasConstant*cp_residual // multiply with R to get dimension [J/mol K]
             );

    return heatCapacityConstantPressure; // [J/mol K]
}
