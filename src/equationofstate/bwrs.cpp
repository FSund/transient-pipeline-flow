#include "equationofstate/bwrs.hpp"

#include <memory>
#include <cmath>

#include "utilities/utilities.hpp"
#include "utilities/errors.hpp"
#include "utilities/physics.hpp"
#include "utilities/stringbuilder.hpp"
#include "constants.hpp"

using utils::pow2;
using utils::pow3;
using utils::pow4;
using utils::pow5;
using utils::pow6;

using arma::vec;
using arma::mat;
using arma::zeros;
using arma::uvec;
using arma::uword;
using std::string;

BWRS::BWRS(const vec& composition, const string& parameterSet):
    EquationOfStateBase(composition)
{
    loadParametersAndCriticalProperties(parameterSet);
}

BWRS BWRS::fromFilePaths(
        const vec& composition,
        const string& ABparameterFile,
        const string& binaryInteractionTableFile,
        const string& criticalProperties)
{
    BWRS eos(composition);
    // NB: important to load parameter files before critical properties, since
    // loading new critical properties triggers re-evaluation of some other
    // properties, which depend on the AB and kij parameters
    eos.loadParameterFiles(ABparameterFile, binaryInteractionTableFile); // this can crash/throw error
    eos.loadCriticalProperties(criticalProperties);

    return eos;
}

void BWRS::loadCriticalProperties(const std::string name)
{
    if (name == "Calsep")
    {
        this->loadJFHCriticalProperties();
    }
    else if (name == "JFH")
    {
        this->loadCalsepCriticalProperties();
    }
    else if (name == "Starling")
    {
        this->loadStarlingCriticalProperties();
    }
    else
    {
        throw std::invalid_argument("name");
    }
}

void BWRS::loadParametersAndCriticalProperties(const std::string parameterSet)
{
    // NB! Need to load AB and kij before critical properties, since
    // loading new critical properties triggers re-evaluation of some other
    // properties, which depend on the AB and kij parameters
    if (parameterSet == "Calsep")
    {
        const string ABparameterFile = string(TRANSFLOW_RESOURCE_PATH) + "/equationofstate/bwrs/calsepABparameters.csv";
        const string binaryInteractionTableFile = string(TRANSFLOW_RESOURCE_PATH) + "/equationofstate/bwrs/calsepBinaryInteraction.csv";
        this->loadParameterFiles(ABparameterFile, binaryInteractionTableFile);
        this->loadStarlingCriticalProperties();
    }
    else if (parameterSet == "JFH")
    {
        const string ABparameterFile = string(TRANSFLOW_RESOURCE_PATH) + "/equationofstate/bwrs/calsepABparameters.csv";
        const string binaryInteractionTableFile = string(TRANSFLOW_RESOURCE_PATH) + "/equationofstate/bwrs/JFH_binaryInteraction.csv";
        this->loadParameterFiles(ABparameterFile, binaryInteractionTableFile);
        this->loadCalsepCriticalProperties();
    }
    else if (parameterSet == "Starling")
    {
        const string ABparameterFile = string(TRANSFLOW_RESOURCE_PATH) + "/equationofstate/bwrs/StarlingABparameters.csv";
        const string binaryInteractionTableFile = string(TRANSFLOW_RESOURCE_PATH) + "/equationofstate/bwrs/StarlingBinaryInteraction.csv";
        this->loadParameterFiles(ABparameterFile, binaryInteractionTableFile);
        this->loadStarlingCriticalProperties();
    }
    else
    {
        throw std::invalid_argument("parameterSet");
    }
}

// TODO: Document the formulas and derivatives below a bit better.
// Perhaps via an attached document? Think I have the notes somewhere.
vec BWRS::evaluate(
        const double pressure,
        const double temperature) const
{
    const double rho_m = findMolarDensity(pressure, temperature);

    const double gasConst = m_R; // gas constant

    const double pressure2 = pressure*pressure;
    const double pressure3 = pressure2*pressure;
    const double pressure4 = pressure3*pressure;
    const double pressure5 = pressure4*pressure;

    const double temperature2 = temperature*temperature;
    const double temperature3 = temperature2*temperature;
    const double temperature4 = temperature3*temperature;
    const double temperature5 = temperature4*temperature;
    const double temperature6 = temperature5*temperature;
    const double temperature7 = temperature6*temperature;
    const double temperature8 = temperature7*temperature;

    const double rho_m2 = rho_m*rho_m;
    const double rho_m5 = rho_m2*rho_m2*rho_m;

    const double gasConst2 = gasConst*gasConst;
    const double gasConst3 = gasConst2*gasConst;
    const double gasConst5 = gasConst3*gasConst*gasConst;
    const double gasConst6 = gasConst5*gasConst;

    const double gasConstT = gasConst*temperature;
    const double gasConstT2 = gasConstT*gasConstT;
    const double gasConstT3 = gasConstT2*gasConstT;
    const double gasConstT5 = gasConstT3*gasConstT*gasConstT;
    const double gasConstT6 = gasConstT5*gasConstT;

    vec Z_factor = zeros<vec>(6); // output vector

    // Z factor
    Z_factor(0) = pressure/(rho_m*gasConst*temperature);

    // dzdT (derivative of Z w.r.t. temperature, at constant pressure)
    const double col = Z_factor(0);
    const double col2 = col*col;
    const double col3 = col2*col;
    const double col4 = col3*col;
    const double col5 = col4*col;

    const double colgasConstT2 = pow((col*gasConst*temperature), 2);

    const double expo = exp(-m_GAMMA*rho_m2);
    const double nom1 =
          (         -m_B0/(gasConst*temperature2) + 2*m_A0/(gasConst2*temperature3)                                         + 4*m_C0/(gasConst2*temperature5) - 5*m_D0/(gasConst2*temperature6) + 6*m_E0/(gasConst2*temperature7)) *col4*pressure
        + (                                       - 2*m_b/(gasConst2*temperature3) +  3*m_a/(gasConst3*temperature4) +  4*m_d/(gasConst3*temperature5)                                                                                ) *col3*pressure2
        + m_ALPHA*(-6*m_a/(gasConst6*temperature7) - 7*m_d/(gasConst6*temperature8)) *pressure5
        + col3*m_c*pressure2/gasConst3 * (-5/temperature6 - 7*m_GAMMA*pressure2/(col2*gasConst2*temperature8)) * expo
        + col3*m_c*pressure2/gasConst3 * ( 1/temperature5 +   m_GAMMA*pressure2/(col2*gasConst2*temperature7)) * expo*2*m_GAMMA*pressure2/(colgasConstT2*temperature);

    const double den1 =
            6*col5
            - 5*col4
            - ( m_B0/(gasConst*temperature)         - m_A0/(gasConst2*temperature2)                                     - m_C0/(gasConst2*temperature4) + m_D0/(gasConst2*temperature5) - m_E0/(gasConst2*temperature6) )*4*col3*pressure
            - (                                m_b/(gasConst2*temperature2) - m_a/(gasConst3*temperature3) - m_d/(gasConst3*temperature4)                                                                                                            )*3*col2*pressure2
            - 3*col2*m_c*pressure2/(gasConst3*temperature5)*expo
            - m_GAMMA*pressure4*m_c/(gasConst5*temperature7)*expo
            - col3*m_c*pressure2/gasConst3*(1/temperature5 + m_GAMMA*pressure2/(col2*gasConst2*temperature7))*expo*2*m_GAMMA*pressure2/(colgasConstT2*col);
    Z_factor(1) = nom1/den1;

    // dzdp_T (derivative of Z w.r.t. pressure, at constant temperature)
    const double nom2 = (m_B0*gasConst*temperature - m_A0 - m_C0/temperature2 + m_D0/temperature3 - m_E0/temperature4) * col4/gasConstT2 + (m_b*gasConst*temperature - m_a - m_d/temperature)*col3*2*pressure/gasConstT3
            + m_ALPHA*(m_a + m_d/temperature)*5*pressure4/gasConstT6 + expo*(m_c*2*pressure*col3/(temperature2*gasConstT3) + 4*pressure3*m_c*col*m_GAMMA/(temperature2*gasConstT5))
            - m_GAMMA*2*pressure/colgasConstT2*expo*(m_c*pressure2*col3/(temperature2*gasConstT3) + m_c*pressure4*col*m_GAMMA/(temperature2*gasConstT5));
    const double den2 = 6*col5 - 5*col4 - (m_B0*gasConst*temperature - m_A0 - m_C0/temperature2 + m_D0/temperature3 - m_E0/temperature4)*4*col3*pressure/gasConstT2
            - (m_b*gasConst*temperature - m_a - m_d/temperature)*3*col2*pressure2/gasConstT3 - expo*(3*col2*m_c*pressure2/(temperature2*gasConstT3) + m_c*pressure4*m_GAMMA/(temperature2*gasConstT5))
            - expo*2*m_GAMMA*pressure2/(col3*gasConstT2)*(m_c*pressure2*col3/(temperature2*gasConstT3) + m_c*pressure4*col*m_GAMMA/(temperature2*gasConstT5));
    Z_factor(2) = nom2/den2;

    // dzdT_rho (derivative of Z w.r.t. temperature, at constant density)
    Z_factor(3) = m_A0*rho_m/(gasConst*temperature2)
            + 3*m_C0*rho_m/(gasConst*temperature4)
            - 4*m_D0*rho_m/(gasConst*temperature5)
            + 5*m_E0*rho_m/(gasConst*temperature6)
            + m_a*rho_m2/(gasConst*temperature2)
            + 2*m_d*rho_m2/(gasConst*temperature3)
            - m_ALPHA*m_a*rho_m5/(gasConst*temperature2)
            - 2*m_ALPHA*m_d*rho_m5/(gasConst*temperature3)
            - 3*m_c*rho_m2/(gasConst*temperature4)*(1 + m_GAMMA*rho_m2)*(exp( - m_GAMMA*rho_m2));

    if (m_useConstantHeatCapacities)
    {
        Z_factor(4) = 3000; // cp
        Z_factor(5) = 1750; // cv
//        Z_factor(5) = this->calculateHeatCapacityConstantPressureTGNet(m_molarMassOfMixture, pressure, temperature);
    }
    else
    {
        // the incorrect(?) Langelandsvik equation actually seems to be most
        // consistent with gerg c_p
        const double cp = utils::calculateHeatCapacityConstantPressureLangelandsvik(m_molarMassOfMixture, pressure, temperature);

        // both JFH and TGNet versions are bad in their own ways...
        const double cv = utils::calculateHeatCapacityConstantVolumeJFH(pressure);
//        const double cv = utils::calculateHeatCapacityConstantVolumeTGNet(m_molarMassOfMixture, pressure, temperature);

        Z_factor(4) = cp;
        Z_factor(5) = cv;
    }

    return Z_factor;
}

double BWRS::calculateCompressibility(
        const double pressure,
        const double temperature) const
{
    const double rho_m = findMolarDensity(pressure, temperature);
    double Z = pressure/(rho_m*m_R*temperature);

    return Z;
}

bool BWRS::setComposition(const vec& composition, const bool force)
{
    if (EquationOfStateBase::setComposition(composition, force))
    {
        findNonZeroComponents();

        m_criticalPressureOfMixture = sum(m_composition % m_pc);
        m_criticalTemperatureOfMixture = sum(m_composition % m_Tc);

        calculateCoefficients();

        return true; // composition changed
    }
    else
    {
        return false; // composition not changed
    }
}

void BWRS::loadJFHCriticalProperties()
{
    // parameters from JFH code
    //                  C1          C2          C3          iC4         nC4         iC5         nC5         C6          N2          CO2
    m_Tc        = vec({ 190.56,     305.32,     369.83,     407.82,     425.13,     460.35,     469.7,      507.82,     126.19,     304.13});           // critical temperature
    m_rhoc      = vec({ 10.14*16.04,6.87*30.07, 5*44.1,     3.86*58.12, 3.92*58.12, 3.27*72.15, 3.21*72.15, 2.71*86.18, 11.18*28.01,10.62*44.01});      // critical density
    m_w         = vec({ 0.0115,     0.099,      0.153,      0.1756,     0.19,       0.22,       0.25,       0.281,      0.037,      0.22});             // accentric factor
    m_pc        = vec({ 45.96,      48.839,     42.5,       36.48,      37.96,      33.81,      33.69,      30.2,       33.99,      73.825})*1e5;       // critical pressure
    m_molarMass = vec({ 16.04,      30.07,      44.1,       58.12,      58.12,      72.15,      72.15,      86.18,      28.13,      44.01});     // molar weight of components [g/mol]

    m_expW = exp(-3.8*m_w);
    m_R = constants::gasConstant; // gas constant [m3 Pa / K mol]

    setComposition(m_composition); // to update m_molarMassOfMixture if m_molarMass has changed, also updates m_criticalPressureOfMixture and m_criticalTemperatureOfMixture

    std::cout << "WARNING: Using Calsep critical parameters, these have been shown to give bad results with the new molar-density based solver." << std::endl;
}

void BWRS::loadCalsepCriticalProperties()
{
    // parameters from Calsep report, with averaged C6+ critical properties
    //                  C1          C2          C3          iC4         nC4         iC5         nC5         C6          N2          CO2
    m_Tc        = vec({ 190.56,     305.32,     369.83,     407.82,     425.13,     460.35,     469.7,      530.3,      126.19,     304.13});           // critical temperature
    m_rhoc      = vec({ 1.014e+04,  6.87e+03,   5.0+03,     3.86e+03,   3.92e+03,   3.27e+03,   3.21e+03,   2.82664e+3, 1.11323e+4, 1.062e+04});        // critical molar density [mol / m3] -- using rho_c from line above, converted using molar mass ([16.04, 30.07, ...])
    m_w         = vec({ 0.0115,     0.099,      0.153,      0.1756,     0.19,       0.22,       0.25,       0.339,      0.037,      0.22});             // accentric factor
    m_pc        = vec({ 45.96,      48.839,     42.5,       36.48,      37.96,      33.81,      33.69,      27.34,      33.99,      73.825})*1e5;       // critical pressure
    m_molarMass = vec({ 16.04,      30.07,      44.1,       58.12,      58.12,      72.15,      72.15,      86.18,      28.13,      44.01});     // molar weight of components [g/mol]

    m_expW = exp(-3.8*m_w);
    m_R = constants::gasConstant; // gas constant [m3 Pa / K mol]

    setComposition(m_composition); // to update m_molarMassOfMixture if m_molarMass has changed, also updates m_criticalPressureOfMixture and m_criticalTemperatureOfMixture

    std::cout << "WARNING: Using Calsep critical parameters, these have been shown to give bad results with the new molar-density based solver." << std::endl;
}

void BWRS::loadStarlingCriticalProperties()
{
    // parameters from Starling book
    // these are converted from Fahrenheit (critical temperature), lb-mol/ft^3 (critical density) and psia (pressure) -- accentric factor has no unit
    //                  C1          C2          C3          iC4         nC4         iC5         nC5         C6          N2          CO2
    m_Tc        = vec({ 190.69,     305.39,     369.89,     408.13,     425.19,     460.37,     469.49,     507.29,     126.15,     304.15});           // critical temperature
    m_rhoc      = vec({ 1.00500e+4, 6.75659e+3, 4.99936e+3, 3.80118e+3, 3.92132e+3, 3.24694e+3, 3.21491e+3, 2.71673e+3, 1.10992e+4, 1.06379e+4});       // critical molar density [mol/m3]
    m_w         = vec({ 0.013,      0.1018,     0.157,      0.183,      0.197,      0.226,      0.252,      0.302,      0.035,      0.21});             // accentric factor
    m_pc        = vec({ 45.96,      48.839,     42.5,       36.48,      37.96,      33.81,      33.69,      27.34,      33.99,      73.825})*1e5;       // critical pressure
    m_molarMass = vec({ 16.042,     30.068,     44.094,     58.12,      58.12,      72.146,     72.146,     86.172,     28.016,     44.01}); // molar mass of components

    m_expW = exp(-3.8*m_w);
    m_R = 8.3160; // gas constant from Starling [m3 Pa / K mol]

    setComposition(m_composition); // to update m_molarMassOfMixture if m_molarMass has changed, also updates m_criticalPressureOfMixture and m_criticalTemperatureOfMixture
}

void BWRS::loadGasscoParameters()
{
    const string ABparameterFile = string(TRANSFLOW_RESOURCE_PATH) + "/equationofstate/bwrs/calsepABparameters.csv";
    const string binaryInteractionTableFile = string(TRANSFLOW_RESOURCE_PATH) + "/equationofstate/bwrs/JFH_binaryInteraction.csv";

    loadParameterFiles(ABparameterFile, binaryInteractionTableFile);
    calculateCoefficients();
}

void BWRS::loadCalsepParameters()
{
    const string ABparameterFile = string(TRANSFLOW_RESOURCE_PATH) + "/equationofstate/bwrs/calsepABparameters.csv";
    const string binaryInteractionTableFile = string(TRANSFLOW_RESOURCE_PATH) + "/equationofstate/bwrs/calsepBinaryInteraction.csv";

    loadParameterFiles(ABparameterFile, binaryInteractionTableFile);
    calculateCoefficients();
}

void BWRS::loadStarlingParameters()
{
    const string ABparameterFile = string(TRANSFLOW_RESOURCE_PATH) + "/equationofstate/bwrs/StarlingABparameters.csv";
    const string binaryInteractionTableFile = string(TRANSFLOW_RESOURCE_PATH) + "/equationofstate/bwrs/StarlingBinaryInteraction.csv";

    loadParameterFiles(ABparameterFile, binaryInteractionTableFile);
    calculateCoefficients();
}

void BWRS::calculateCoefficients()
{
    const vec& x = m_composition;
    const vec& w = m_w;
    const double R = m_R; // gas constant [m3 Pa / K mol]
    const vec& Tc = m_Tc;
    const vec& rhoc = m_rhoc;
    const vec& Ai = m_Ai;
    const vec& Bi = m_Bi;
    const mat& k = m_binaryInteractionParameterTable;

    double A0 = 0;
    double C0 = 0;
    double D0 = 0;
    double E0 = 0;
    for (uword i : m_indices)
    {
        double A0i =  (Ai(1) + Bi(1)*w(i)) *R       *Tc(i) /     rhoc(i);
        double C0i =  (Ai(2) + Bi(2)*w(i)) *R  *pow3(Tc(i))/     rhoc(i);
        double D0i =  (Ai(8) + Bi(8)*w(i)) *R  *pow4(Tc(i))/     rhoc(i);
        double E0i = (Ai(10) + Bi(10)*w(i)*m_expW(i)) *R  *pow5(Tc(i))/     rhoc(i);
        for (uword j : m_indices)
        {
            double A0j =  (Ai(1) + Bi(1)*w(j)) *R       *Tc(j) /     rhoc(j);
            double C0j =  (Ai(2) + Bi(2)*w(j)) *R  *pow3(Tc(j))/     rhoc(j);
            double D0j =  (Ai(8) + Bi(8)*w(j)) *R  *pow4(Tc(j))/     rhoc(j);
            double E0j = (Ai(10) + Bi(10)*w(j)*m_expW(j)) *R  *pow5(Tc(j))/     rhoc(j);

            A0 += x(i)*x(j)    *(1.0 - k(i,j))*sqrt(A0i)*sqrt(A0j);
            C0 += x(i)*x(j)*pow3(1.0 - k(i,j))*sqrt(C0i)*sqrt(C0j);
            D0 += x(i)*x(j)*pow4(1.0 - k(i,j))*sqrt(D0i)*sqrt(D0j);
            E0 += x(i)*x(j)*pow5(1.0 - k(i,j))*sqrt(E0i)*sqrt(E0j);
        }
    }

    double B0 = 0;
    double a = 0;
    double b = 0;
    double c = 0;
    double d = 0;
    double alpha = 0;
    double gamma = 0;
    for (uword i : m_indices)
    {
        double B0i =    (Ai(0) + Bi(0)*w(i)) *      1.0      /     rhoc(i);
        double ai =     (Ai(5) + Bi(5)*w(i)) *R       *Tc(i) /pow2(rhoc(i));
        double bi =     (Ai(4) + Bi(4)*w(i)) *1.0            /pow2(rhoc(i));
        double ci =     (Ai(7) + Bi(7)*w(i)) *R  *pow3(Tc(i))/pow2(rhoc(i));
        double di =     (Ai(9) + Bi(9)*w(i)) *R  *pow2(Tc(i))/pow2(rhoc(i));
        double alphai = (Ai(6) + Bi(6)*w(i)) *1.0            /pow3(rhoc(i));
        double gammai = (Ai(3) + Bi(3)*w(i)) *1.0            /pow2(rhoc(i));

        B0 += x(i)*B0i;

        a += x(i)*cbrt(ai);
        b += x(i)*cbrt(bi);
        c += x(i)*cbrt(ci);
        d += x(i)*cbrt(di);

        alpha += x(i)*cbrt(alphai);
        gamma += x(i)*sqrt(gammai);
    }
    a = pow3(a);
    b = pow3(b);
    c = pow3(c);
    d = pow3(d);
    alpha = pow3(alpha);
    gamma = pow2(gamma);

    // all these coefficients only depend on the gas composition, so they are
    // calculated once here and store for later use when evaluating BWRS at
    // different temperatures and pressures
    m_A0 = A0;
    m_B0 = B0;
    m_C0 = C0;
    m_D0 = D0;
    m_E0 = E0;
    m_a = a;
    m_b = b;
    m_c = c;
    m_d = d;
    m_ALPHA = alpha;
    m_GAMMA = gamma;
}

void BWRS::loadParameterFiles(
        const string& ABparameterFile,
        const string& binaryInteractionTableFile)
{
    // load binary interaction parameters and Ai, Bi
    mat loading;

    loading.load(ABparameterFile, arma::file_type::csv_ascii);
    m_Ai = loading.col(0);
    m_Bi = loading.col(1);

    loading.load(binaryInteractionTableFile, arma::file_type::csv_ascii);
    m_binaryInteractionParameterTable = loading;
}

void BWRS::findNonZeroComponents()
{
    m_indices = {};
    for (uword i = 0; i < m_composition.n_elem; i++)
    {
        if (m_composition(i) > 2.0*std::numeric_limits<double>::epsilon()) // x should be 0 <= x <= 1, so this comparison should be okay
        {
            m_indices = join_vert(m_indices, uvec({i}));
        }
    }
}

void BWRS::enableConstantHeatCapacities()
{
    m_useConstantHeatCapacities = true;
}

double BWRS::findDensity(
        const double pressure,
        const double temperature,
        const double tolerance) const
{
    const double rho_m = findMolarDensity(pressure, temperature, tolerance);

    // convert from reduced density [mol / m3] to actual density [kg / m3]
    // molar mass of mixture has unit [g/mol]
    const double density = rho_m*m_molarMassOfMixture/1000.0;

    m_density = density; // cache density (mutable)

    return density;
}

// private
double BWRS::findMolarDensity(
        const double pressure,
        const double temperature,
        const double tolerance) const
{
    // use Newtons method to find the density of the gas

    // it should not be necessary to check this here...
    if (temperature <= 0)
    {
        return 0;
    }

    if (pressure <= 0)
    {
        return 0;
    }

    double initialDensity;
    if (m_density == 0)
    {
        // use density of ideal gas as starting point
        const double specificGasConstant = constants::gasConstant/(m_molarMassOfMixture/1000.0 /*g->kg*/); // J/(kg*K)
        const double idealGasDensity = pressure/(specificGasConstant*temperature);
        initialDensity = idealGasDensity;
    }
    else
    {
        // use cached density
        initialDensity = m_density;
    }

    // convert density guess from actual density [kg / m3] to reduced density [mol / m3]
    // (molar mass of mixture has unit [g/mol])
    const double molarDensityGuess = initialDensity/(m_molarMassOfMixture/1000.0);

    double current = molarDensityGuess;
    double previous;

    uword counter = 0;
    uword nMax = 1000;
    while (counter < nMax)
    {
        double temp = temperature;
        double temp2 = temp*temp;
        double temp3 = temp2*temp;
        double temp4 = temp2*temp2;

        double current2 = current*current;
        double current3 = current2*current;
        double current4 = current2*current2;
        double current5 = current2*current3;
        double current6 = current3*current3;

        double expGammaCurrent2 = exp(-m_GAMMA*current2);

        double f = - pressure + current*m_R*temp
                + ( m_B0*m_R*temperature - m_A0 - m_C0/temp2 + m_D0/temp3 - m_E0/temp4 )*current2
                + ( m_b*m_R*temp - m_a - m_d/temp )*current3
                + m_ALPHA*( m_a + m_d/temp )*current6
                + ( (m_c*current3)/temp2 )*( 1 + m_GAMMA*current2 )*( expGammaCurrent2 );

        double fdiv = m_R*temp
                + (m_B0*m_R*temp - m_A0 - m_C0/temp2 + m_D0/temp3 - m_E0/temp4)*2.0*current
                + (m_b*m_R*temp - m_a - m_d/temp)*3.0*current2
                + m_ALPHA*(m_a + m_d/temp)*6.0*current5
                + m_c*current2/temp2*(3.0 + 3.0*m_GAMMA*current2 - 2.0*m_GAMMA*m_GAMMA*current4)*exp(-m_GAMMA*current2);

        previous = current;
        double change = f/fdiv;
        current -= change;

        if (current < 0) // avoid negative
        {
            while (current < 0)
            {
                change /= 2;
                current = previous - change;
            }
            counter++;
            continue; // to skip check below and instead do another iteration
        }

        // check relative change
        if (std::abs(change/previous) < tolerance)
        {
            break;
        }
        counter++;
    }
    if (counter >= nMax)
    {
        throw utils::no_convergence_error(
                    utils::stringbuilder()
                        << "BWRS::findMolarDensity(): Newton's method for finding density did not converge (after "
                        << counter
                        << " iterations).",
                    counter
                    );
    }
    const double rho_m = current;

    return rho_m;
}

//double BWRS::getDerivativeOfZwrtDensityAtConstantTemperature(const double density, const double temperature) const
//{
//    return
//        - m_A0/(m_R*temperature) + m_B0 - m_C0/(m_R*pow(temperature, 3))
//        + m_D0/(m_R*pow(temperature, 4)) - m_E0/(m_R*pow(temperature, 5))
//        + 2*m_b*density + 5*m_a*m_ALPHA*pow(density, 4)/(m_R*temperature)
//        - 2*m_a*density/(m_R*temperature) + 5*m_ALPHA*m_d*pow(density, 4)/(m_R*pow(temperature, 2))
//        - 2*m_d*density/(m_R*pow(temperature, 2))
//        - 2*m_c*pow(m_GAMMA, 2)*pow(density, 5)*exp(-m_GAMMA*pow(density, 2))/(m_R*pow(temperature, 3))
//        + 2*m_c*m_GAMMA*pow(density, 3)*exp(-m_GAMMA*pow(density, 2))/(m_R*pow(temperature, 3))
//        + 2*m_c*density*exp(-m_GAMMA*pow(density, 2))/(m_R*pow(temperature, 3));
//}
