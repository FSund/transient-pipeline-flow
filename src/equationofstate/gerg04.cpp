#include "equationofstate/gerg04.hpp"

#include <cmath>

#include "utilities/utilities.hpp"
#include "utilities/errors.hpp"
#include "utilities/physics.hpp"
#include "utilities/stringbuilder.hpp"
#include "constants.hpp"

using utils::pow2;
using utils::pow3;
using utils::pow5;

using arma::vec;
using arma::uvec;
using arma::mat;
using arma::cube;
using arma::zeros;
using arma::uword;
using std::string;

const string path = string(TRANSFLOW_RESOURCE_PATH) + "/equationofstate/gerg04/";
const arma::file_type fileType = arma::hdf5_binary;

mat GERG04::betav = utils::loadMatFromFile(path + "betav.h5", fileType); // beta density interaction coefficient
mat GERG04::betat = utils::loadMatFromFile(path + "betat.h5", fileType); // beta temperature coefficient
mat GERG04::gammav = utils::loadMatFromFile(path + "gammav.h5", fileType); // gamma density interaction coefficient
mat GERG04::gammat = utils::loadMatFromFile(path + "gammat.h5", fileType); //gamma temperature interaction coefficient

// Kpol coefficients
mat GERG04::noipol = utils::loadMatFromFile(path + "noipol.h5", fileType);
mat GERG04::doipol = utils::loadMatFromFile(path + "doipol.h5", fileType);
mat GERG04::toipol = utils::loadMatFromFile(path + "toipol.h5", fileType);

// noipol, doipol, toipol (1, 1:6) -- C1 methane
// noipol, doipol, toipol (2, 1:6) -- nitrogen
// noipol, doipol, toipol (3, 1:4) -- CO2
// noipol, doipol, toipol (4, 1:6) -- C2 ethane
// noipol, doipol, toipol (5, 1:6) -- C3 propane
// noipol, doipol, toipol (6, 1:6) -- nC4 n-buthane
// noipol, doipol, toipol (7, 1:6) -- iC4 i-buthane
// noipol, doipol, toipol (8, 1:6) -- nC5 n-pentane
// noipol, doipol, toipol (9, 1:6) -- iC5 i-pentane
// noipol, doipol, toipol (10, 1:6) -- nC6 n-hexane

//// Kexp coefficient
mat GERG04::noiexp = utils::loadMatFromFile(path + "noiexp.h5", fileType);
mat GERG04::doiexp = utils::loadMatFromFile(path + "doiexp.h5", fileType);
mat GERG04::coiexp = utils::loadMatFromFile(path + "coiexp.h5", fileType);
mat GERG04::toiexp = utils::loadMatFromFile(path + "toiexp.h5", fileType);

// noiexp, doiexp, coiexp, toiexp (1, 7:24) -- C1 methane
// noiexp, doiexp, coiexp, toiexp (2, 7:24) -- nitrogen
// noiexp, doiexp, coiexp, toiexp (3, 7:24) -- CO2
// noiexp, doiexp, coiexp, toiexp (4, 7:24) -- C2 ethane
// noiexp, doiexp, coiexp, toiexp (5, 7:12) -- C3 propane
// noiexp, doiexp, coiexp, toiexp (6, 7:12) -- nC4 n-butane
// noiexp, doiexp, coiexp, toiexp (7, 7:12) -- iC4 i-butane
// noiexp, doiexp, coiexp, toiexp (8, 7:12) -- nC5 n-pentane
// noiexp, doiexp, coiexp, toiexp (9, 7:12) -- iC5 i-pentane
// noiexp, doiexp, coiexp, toiexp (10, 7:12) -- nC6 n-hexane

mat GERG04::Fij = utils::loadMatFromFile(path + "Fij.h5", fileType); // Fij interaction coefficient

cube GERG04::nijpol = utils::loadCubeFromFile(path + "nijpol.h5", fileType);
cube GERG04::nijexp = utils::loadCubeFromFile(path + "nijexp.h5", fileType);
cube GERG04::dijpol = utils::loadCubeFromFile(path + "dijpol.h5", fileType);
cube GERG04::dijexp = utils::loadCubeFromFile(path + "dijexp.h5", fileType);
cube GERG04::tijpol = utils::loadCubeFromFile(path + "tijpol.h5", fileType);
cube GERG04::tijexp = utils::loadCubeFromFile(path + "tijexp.h5", fileType);

cube GERG04::nuijexp = utils::loadCubeFromFile(path + "nuijexp.h5", fileType);
cube GERG04::epijexp = utils::loadCubeFromFile(path + "epijexp.h5", fileType);
cube GERG04::beijexp = utils::loadCubeFromFile(path + "beijexp.h5", fileType);
cube GERG04::gaijexp = utils::loadCubeFromFile(path + "gaijexp.h5", fileType);

// other coefficients
mat GERG04::noik = utils::loadMatFromFile(path + "noik.h5", fileType);
mat GERG04::voik = utils::loadMatFromFile(path + "voik.h5", fileType);

// constant terms
cube GERG04::nijpol_times_tijpol = GERG04::nijpol % GERG04::tijpol;
cube GERG04::nijpol_times_tijpol_times_tijpol_minus_one =
        (GERG04::nijpol % GERG04::tijpol) % (GERG04::tijpol - 1.0);

cube GERG04::nijexp_times_tijexp = GERG04::nijexp % GERG04::tijexp;
cube GERG04::nijexp_times_tijexp_times_tijexp_minus_one =
        GERG04::nijexp % GERG04::tijexp % (GERG04::tijexp - 1.0);

cube GERG04::nijpol_times_dijpol =
        GERG04::nijpol % GERG04::dijpol;
cube GERG04::nijpol_times_dijpol_times_tijpol =
        GERG04::nijpol % GERG04::dijpol % GERG04::tijpol;
cube GERG04::nijpol_times_dijpol_times_dijpol_minus_one =
        GERG04::nijpol % GERG04::dijpol % (GERG04::dijpol - 1.0);

GERG04::GERG04(const vec& composition):
    EquationOfStateBase(composition)
{
    GERG04::setComposition(composition);
}

void GERG04::setNonZeroComponents()
{
    m_indices = arma::find(X > 2.0*std::numeric_limits<double>::epsilon());
    m_firstIndices = m_indices.head(m_indices.n_elem - 1); // [0, N-1]
    m_lastIndices = m_indices.tail(m_indices.n_elem - 1); // [1, N]
}

void GERG04::calculateCoefficients()
{
    const double R = 8.314472; // [J/mol K] gas constant from TM15 page XII
    Ra = R*1000/m_molarMassOfMixture; // gas constant of mixture

    // find reducing function for mixture density and temperature
    // using eq. (7.9) and eq. (7.10)
    rhored = 0; // (inverse) reducing function for mixture density == 1/\rho_r(\bar x)
    tred = 0; // reducing function for mixture temperature == T_r(\bar x)
    for (uword i : m_indices)
    {
        rhored += pow2(X(i))/rhoc(i);
        tred += pow2(X(i))*Tc(i);
    }
    for (uword i : m_firstIndices)
    {
        for (uword j : m_lastIndices)
        {
            rhored +=
                2.0*X(i)*X(j)*betav(i,j)*gammav(i,j)*(X(i) + X(j))/(pow2(betav(i,j))*X(i) + X(j))*(1.0/8.0)*pow3(1.0/(pow(rhoc(i), 1.0/3.0)) + 1.0/(pow(rhoc(j), 1.0/3.0)));
            tred +=
                2.0*X(i)*X(j)*betat(i,j)*gammat(i,j)*(X(i) + X(j))/(pow2(betat(i,j))*X(i) + X(j))*sqrt(Tc(i)*Tc(j));
        }
    }
}

bool GERG04::setComposition(const vec& composition, const bool force)
{
    if (EquationOfStateBase::setComposition(composition, force))
    {
        const vec& comp = m_composition;

        //   ch4,       n2,       c02,      c2h6,     c3h8,    nc4h10,   ic4h10,   nc5h12,   ic5h12   nc6h14
        this->X = {comp(0),  comp(8),  comp(9),  comp(1),  comp(2),  comp(4),  comp(3),  comp(6),  comp(5), comp(7)};

        setNonZeroComponents();

        calculateCoefficients();

        return true; // composition changed
    }
    else
    {
        return false; // composition not changed
    }
}

vec GERG04::evaluate(const double pressure, const double temperature) const
{
    vec Z_all = evaluateAllProperties(pressure, temperature);
    vec Z_factor = zeros<vec>(6);
    Z_factor(0) = Z_all(0); // Z
    Z_factor(1) = Z_all(1);
    Z_factor(2) = Z_all(2);
    Z_factor(3) = Z_all(3);
    if (false)
    {
        Z_factor(4) = utils::calculateHeatCapacityConstantPressureTGNet(m_molarMassOfMixture, pressure, temperature); // c_p
        Z_factor(5) = utils::calculateHeatCapacityConstantVolumeJFH(pressure); // c_v
    }
    else
    {
        Z_factor(4) = Z_all(8); // c_p
        Z_factor(5) = Z_all(6); // c_v
    }

    return Z_factor;
}

double GERG04::calculateCompressibility(const double pressure, const double temperature) const
{
    double density = findDensity(pressure, temperature);
    double Z = pressure/(density*Ra*temperature);

    return Z;
}

double GERG04::findDensity(const double pressure, const double temperature) const
{
    double density = 0;
    double aroidelta = 0;
    double arijdelta = 0;
    double aroideltadelta = 0;
    double arijdeltadelta = 0;

    return findDensity(pressure, temperature, density, aroidelta, arijdelta, aroideltadelta, arijdeltadelta);
}

double GERG04::findDensity(
        const double pressure,
        const double temperature,
        double& density, // output
        double& aroidelta, // output
        double& arijdelta, // output
        double& aroideltadelta, // output
        double& arijdeltadelta // output
        ) const
{
    // compressibility factor calculations
    // find density using Newton's method

    // final alpha coefficient
    density = 0;
    aroidelta = 0;
    arijdelta = 0;
    aroideltadelta = 0;
    arijdeltadelta = 0;

    const double tred_temperature = tred/temperature;

    double initialDensity;
    if (m_density <= 0)
    {
        // use density of ideal gas as starting point
        const double specificGasConstant = constants::gasConstant/(m_molarMassOfMixture/1000.0 /*g->kg*/); // J/(kg*K)
        const double idealGasDensity = pressure/(specificGasConstant*temperature);
        if (idealGasDensity <= 0)
            throw std::runtime_error("ideal gas density <= 0");
        initialDensity = idealGasDensity;
    }
    else
    {
        // use cached density
        initialDensity = m_density;
    }

    const bool noInitialDensity = true;
    double previousDensity = initialDensity;
    double pow2_density = previousDensity*previousDensity;
    double density_rhored = previousDensity*rhored;
    double density_rhored_pow_minusOne = 1.0/density_rhored;
    double density_rhored_pow_minusTwo = 1.0/pow2(density_rhored);

    uword firstSORlimit = 10;
    if (noInitialDensity) firstSORlimit = 50;
    uword secondSORlimit = 50;
    if (noInitialDensity) secondSORlimit = 200;
    uword maxIterations = 300;
    if (noInitialDensity) maxIterations = 500;

    vec densityIterations = zeros<vec>(maxIterations + 2);
    densityIterations(0) = previousDensity;

    bool run=true;
    uword counter = 0;
    while (run)
    {
        this->evaluateAlpha_roi_deltas(tred_temperature, density_rhored, density_rhored_pow_minusOne, density_rhored_pow_minusTwo, aroidelta, arijdelta, aroideltadelta, arijdeltadelta);
        const double diff = (previousDensity + pow2_density*rhored*(aroidelta + arijdelta) - pressure/(Ra*temperature))/(1 + 2*density_rhored*(aroidelta + arijdelta) + pow2_density*rhored*(aroideltadelta + arijdeltadelta));

        if (counter > firstSORlimit)
        {
            // use over-relaxation
            density = previousDensity - 1.5*diff;
        }
        else if (counter > secondSORlimit)
        {
            // use over-relaxation
            density = previousDensity - 1.9*diff;
        }
        else
        {
            density = previousDensity - diff;
        }

        if (std::abs(diff/previousDensity) < 0.0001) // relative difference less than 0.01 %
        {
            run = false;
        }

        previousDensity = density;
        pow2_density = previousDensity*previousDensity;
        density_rhored = previousDensity*rhored;
        density_rhored_pow_minusOne = 1.0/density_rhored;
        density_rhored_pow_minusTwo = 1.0/pow2(density_rhored);
        counter++;

        densityIterations(counter) = previousDensity;

        if (counter > maxIterations)
        {
            throw utils::no_convergence_error(
                        utils::stringbuilder()
                        << "GERG04::findDensity(): Newton's method for finding density did not converge (after "
                        << counter
                        << " iterations).",
                        counter
                        );
        }
    }

    m_density = density; // cache density (mutable)

    return density;
}

void GERG04::evaluateAlpha_roi_deltas(
        const double tred_temperature,
        const double start_rhored,
        const double start_rhored_pow_minusOne,
        const double start_rhored_pow_minusTwo,
        double& aroidelta, // output
        double& arijdelta, // output
        double& aroideltadelta, // output
        double& arijdeltadelta // output
        ) const
{
    aroidelta=0;
    arijdelta=0;
    aroideltadelta=0;
    arijdeltadelta=0;
    for (uword i : m_indices)
    {
        double aroidelta_part = 0;
        double aroideltadelta_part = 0;
        for (uword k = 0; k < 24; k++)
        {
            const double pow_t_toipol = pow(tred_temperature, toipol(i,k));
            const double noipoldoipol = noipol(i,k)*doipol(i,k);
            const double powstartrhoredcoiexp = pow(start_rhored, coiexp(i,k));
            const double doiexcoiexppow = doiexp(i,k)-coiexp(i,k)*powstartrhoredcoiexp;
            const double pow_start_rhored_doipol = pow(start_rhored, doipol(i,k));
            const double pow_start_rhored_doiexp = pow(start_rhored, doiexp(i,k));

            aroidelta_part +=
                    +(noipoldoipol*pow_start_rhored_doipol*start_rhored_pow_minusOne*pow_t_toipol)
                    +(noiexp(i,k)*pow_start_rhored_doiexp*start_rhored_pow_minusOne*(doiexcoiexppow)*pow(tred_temperature, toiexp(i,k))*exp(-(start_rhored)*coiexp(i,k)));
            aroideltadelta_part +=
                    +noipoldoipol*(doipol(i,k)-1)*pow_start_rhored_doipol*start_rhored_pow_minusTwo*pow_t_toipol
                    +noiexp(i,k)*pow_start_rhored_doiexp*start_rhored_pow_minusTwo*
                    (
                        (doiexcoiexppow)
                        *(doiexp(i,k)-1-coiexp(i,k)*powstartrhoredcoiexp)-pow2(coiexp(i,k))*powstartrhoredcoiexp
                    // this line contains an error when comparing to the
                    // expression in Table 7.7 in TM15 --
                    // (tred_temperature) should be replaced by
                    // pow(tred_temperature, toiexp(i,k))
                    // but testing this leads to very bad results, so I'm
                    // not sure what is going on
                    )*(tred_temperature)*exp(-powstartrhoredcoiexp);
        }
        aroidelta += aroidelta_part*X(i);
        aroideltadelta += aroideltadelta_part*X(i);
    }
    for (uword i : m_firstIndices)
    {
        for (uword j : m_lastIndices)
        {
            double arijdelta_part = 0;
            double arijdeltadelta_part = 0;
            for (uword k = 0; k < 20; k++)
            {
                const double powtredtijpol = pow(tred_temperature, tijpol(i,j,k));
                const double powtredtijexp = pow(tred_temperature, tijexp(i,j,k));
                const double pow_start_rhored_dijexp = pow(start_rhored, dijexp(i,j,k));
                const double pow_start_rhored_dijpol = pow(start_rhored, dijpol(i,j,k));
                const double WHAT = exp(-nuijexp(i,j,k)*pow2(start_rhored - epijexp(i,j,k)) - beijexp(i,j,k)*(start_rhored - gaijexp(i,j,k)));
                const double HEI = dijexp(i,j,k)/(start_rhored) - 2*nuijexp(i,j,k)*(start_rhored - epijexp(i,j,k)) - beijexp(i,j,k);

                arijdelta_part +=
                        nijpol_times_dijpol(i,j,k)*pow_start_rhored_dijpol*start_rhored_pow_minusOne*powtredtijpol
                        +nijexp(i,j,k)*pow_start_rhored_dijexp*powtredtijexp*WHAT*HEI;

                arijdeltadelta_part +=
                        nijpol_times_dijpol_times_dijpol_minus_one(i,j,k)*pow_start_rhored_dijpol*start_rhored_pow_minusTwo*powtredtijpol
                        +nijexp(i,j,k)*pow_start_rhored_dijexp*powtredtijexp*WHAT*(HEI*HEI-dijexp(i,j,k)/pow2(start_rhored)-2*nuijexp(i,j,k));
            }
            arijdelta += arijdelta_part*X(i)*X(j)*Fij(i,j);
            arijdeltadelta += arijdeltadelta_part*X(i)*X(j)*Fij(i,j);
        }
    }
}

vec GERG04::evaluateAllProperties(const double pressure, const double temperature) const
{
    // final alpha coefficient
    double value = 0;
    double aroidelta = 0;
    double arijdelta = 0;
    double aroideltadelta = 0;
    double arijdeltadelta = 0;

    findDensity(pressure, temperature,
                value, aroidelta, arijdelta, aroideltadelta, arijdeltadelta // output parameters
                );

    const double value_rhored = value*rhored;
    const double tred_temperature = tred/temperature;
    const double value_rhored_pow_minusOne = pow(value_rhored,-1);
    const double value_rhored_pow_minusTwo = pow(value_rhored,-2);
//    const double value_rhored_pow_plusOne = pow(value_rhored,+1);
    const double value_rhored_pow_plusOne = value_rhored;
    const double tred_temperatured_pow_minusOne = pow(tred_temperature,-1);
    const double tred_temperatured_pow_minusTwo = pow(tred_temperature,-2);


    // calculating partial derivatives
    double dzdTp=0;
    double dzdTp2=0;
    double dzdTpij=0;
    double dzdTpij2=0;
    double dzdp=0;
    double dzdT_rho=0;
    double dzdT_rho2=0;
    for (uword i : m_indices)
    {
        for (uword k = 0; k < 24; k++)
        {
            const double pow_value_rhored_minusOne = pow(value_rhored, doipol(i,k))*value_rhored_pow_minusOne;
            const double pow_tred_temperature_toipol = pow(tred_temperature, toipol(i,k));
            const double pow_value_rhored_doiexp_minusOne = pow(value_rhored, doiexp(i,k))*value_rhored_pow_minusOne;
            const double pow_tred_temperature_toiexp = pow(tred_temperature, toiexp(i,k));
            const double pow_value_rhored_coiexp = pow(value_rhored, coiexp(i,k));
            const double pow_value_rhored_doiexp_coiexp = pow(value_rhored, doiexp(i,k)+coiexp(i,k)-1);

            //dzdT terms
            dzdTp=
                    dzdTp+X(i)*(noipol(i,k)*doipol(i,k)*(doipol(i,k)-1)*pow_value_rhored_minusOne*pow_tred_temperature_toipol)+X(i)*noiexp(i,k)
                    *(doiexp(i,k)-1)*pow_value_rhored_doiexp_minusOne*doiexp(i,k)*pow_tred_temperature_toiexp*exp(-pow_value_rhored_coiexp)
                    -X(i)*noiexp(i,k)*coiexp(i,k)*(doiexp(i,k)+coiexp(i,k)-1)*pow_value_rhored_doiexp_coiexp*pow_tred_temperature_toiexp
                    *exp(-pow_value_rhored_coiexp)-X(i)*noiexp(i,k)*pow_value_rhored_doiexp_minusOne*(doiexp(i,k)-coiexp(i,k)*pow_value_rhored_coiexp)
                    *pow_tred_temperature_toiexp*coiexp(i,k)*pow_value_rhored_coiexp*exp(-pow_value_rhored_coiexp);
            //right hand side terms
            dzdTp2=
                    dzdTp2-X(i)*noipol(i,k)*doipol(i,k)*pow_value_rhored_minusOne*(doipol(i,k)+toipol(i,k)-1)*pow_tred_temperature_toipol-X(i)*noiexp(i,k)
                    *doiexp(i,k)*(doiexp(i,k)+toiexp(i,k)-1)*pow_value_rhored_doiexp_minusOne*pow_tred_temperature_toiexp
                    *exp(-pow_value_rhored_coiexp)+X(i)*noiexp(i,k)*coiexp(i,k)*(doiexp(i,k)+coiexp(i,k)+toiexp(i,k)-1)*pow_value_rhored_doiexp_coiexp
                    *pow_tred_temperature_toiexp*(exp(-pow_value_rhored_coiexp))+X(i)*noiexp(i,k)*pow_value_rhored_doiexp_minusOne
                    *(doiexp(i,k)-coiexp(i,k)*pow_value_rhored_coiexp)*pow_tred_temperature_toiexp*coiexp(i,k)*pow_value_rhored_coiexp
                    *exp(-pow_value_rhored_coiexp);
            //dzdp terms
            dzdp=
                    dzdp+X(i)*noipol(i,k)*doipol(i,k)*(doipol(i,k)-1)*pow_value_rhored_minusOne*pow_tred_temperature_toipol+X(i)*noiexp(i,k)
                    *(doiexp(i,k)-1)*pow_value_rhored_doiexp_minusOne*doiexp(i,k)*pow_tred_temperature_toiexp*exp(-pow_value_rhored_coiexp)
                    -X(i)*noiexp(i,k)*coiexp(i,k)*(doiexp(i,k)+coiexp(i,k)-1)*pow_value_rhored_doiexp_coiexp*pow_tred_temperature_toiexp
                    *exp(-pow_value_rhored_coiexp)-X(i)*noiexp(i,k)*pow_value_rhored_doiexp_coiexp
                    *(doiexp(i,k)-coiexp(i,k)*pow_value_rhored_coiexp)*pow_tred_temperature_toiexp*coiexp(i,k)*exp(-pow_value_rhored_coiexp);
            //dzdT_rho terms
            dzdT_rho=
                    dzdT_rho+X(i)*(-noipol(i,k)*doipol(i,k)*toipol(i,k)*pow_value_rhored_minusOne*pow_tred_temperature_toipol)/temperature
                    -X(i)*(noiexp(i,k)*toiexp(i,k)*pow_value_rhored_doiexp_minusOne*(doiexp(i,k)-coiexp(i,k)*pow_value_rhored_coiexp)
                           *pow_tred_temperature_toiexp/temperature*exp(-(value_rhored)*coiexp(i,k)));
        }
    }

    const cube pow2_value_rhored_epijexp = (value_rhored - epijexp) % (value_rhored - epijexp);
    const cube exp_nijexp_pow2_value_rhored_epijexp_beijexp_value_rhored_gaijexp_cube = exp(-nuijexp % pow2_value_rhored_epijexp - beijexp % (value_rhored - gaijexp));
    for (uword i : m_firstIndices)
    {
        for (uword j : m_lastIndices)
        {
            const double XXFij = X(i)*X(j)*Fij(i,j);
            for (uword k = 0; k < 20; k++)
            {
                const double pow_value_rhored_dijpol = pow(value_rhored, dijpol(i,j,k)-1);
                const double pow_tred_temperature_dijpol = pow(tred_temperature, tijpol(i,j,k));
                const double pow_tred_temperature_tijexp = pow(tred_temperature, tijexp(i,j,k));

                const double pow_value_rhored_dijexp = pow(value_rhored, dijexp(i,j,k));
                const double pow_value_rhored_dijexp_minusOne = pow_value_rhored_dijexp*value_rhored_pow_minusOne;
                const double pow_value_rhored_dijexp_plusOne = pow_value_rhored_dijexp*value_rhored_pow_plusOne;
//                double exp_nijexp_pow2_value_rhored_epijexp_beijexp_value_rhored_gaijexp = exp(-nuijexp(w,k,l)*pow2(value_rhored-epijexp(w,k,l))-beijexp(w,k,l)*(value_rhored-gaijexp(w,k,l)));

                //dzdT_p terms
                dzdTpij +=
                        +XXFij
                            *nijpol(i,j,k)*dijpol(i,j,k)*(dijpol(i,j,k)-1)*pow_value_rhored_dijpol*pow_tred_temperature_dijpol
                        +XXFij
                            *nijexp(i,j,k)*pow_value_rhored_dijexp_minusOne*(dijexp(i,j,k)-1)*pow_tred_temperature_tijexp*dijexp(i,j,k)
                            *exp_nijexp_pow2_value_rhored_epijexp_beijexp_value_rhored_gaijexp_cube(i,j,k)
                        -XXFij
                            *2*nijexp(i,j,k)
                            *nuijexp(i,j,k)*pow_value_rhored_dijexp_plusOne*(dijexp(i,j,k)+1)*pow_tred_temperature_tijexp
                            *exp_nijexp_pow2_value_rhored_epijexp_beijexp_value_rhored_gaijexp_cube(i,j,k)
                        +XXFij
                            *2*nijexp(i,j,k)
                            *nuijexp(i,j,k)*dijexp(i,j,k)*pow_value_rhored_dijexp*pow_tred_temperature_tijexp*epijexp(i,j,k)
                            *exp_nijexp_pow2_value_rhored_epijexp_beijexp_value_rhored_gaijexp_cube(i,j,k);
                //right hand side terms
                dzdTpij2 +=
                        -XXFij
                            *nijpol(i,j,k)*dijpol(i,j,k)*(dijpol(i,j,k)+tijpol(i,j,k)-1)*pow_value_rhored_dijpol
                            *pow_tred_temperature_dijpol
                        -XXFij
                            *nijexp(i,j,k)*pow_value_rhored_dijexp_minusOne*tijexp(i,j,k)
                            *pow_tred_temperature_tijexp*dijexp(i,j,k)
                            *exp_nijexp_pow2_value_rhored_epijexp_beijexp_value_rhored_gaijexp_cube(i,j,k)
                        +XXFij
                            *2*nijexp(i,j,k)*nuijexp(i,j,k)*pow_value_rhored_dijexp_plusOne*tijexp(i,j,k)
                            *pow_tred_temperature_tijexp
                            *exp_nijexp_pow2_value_rhored_epijexp_beijexp_value_rhored_gaijexp_cube(i,j,k)
                        -XXFij
                            *2*nijexp(i,j,k)*nuijexp(i,j,k)*pow_value_rhored_dijexp*tijexp(i,j,k)*pow_tred_temperature_tijexp
                            *epijexp(i,j,k)*exp_nijexp_pow2_value_rhored_epijexp_beijexp_value_rhored_gaijexp_cube(i,j,k);
                //dzdT_rho terms;
                dzdT_rho2 +=
                        +XXFij
                        *(
                            -nijpol(i,j,k)*toipol(i,j)*dijpol(i,j,k)*pow_value_rhored_dijpol*pow_tred_temperature_dijpol/temperature-nijexp(i,j,k)
                            *tijexp(i,j,k)*pow_value_rhored_dijexp*pow_tred_temperature_tijexp/temperature
                            *exp(-nuijexp(i,j,k)*pow2((value_rhored)-epijexp(i,j,k))-beijexp(i,j,k)*(value_rhored-gaijexp(i,j,k)))
                            *(dijexp(i,j,k)/(value_rhored)-2*nuijexp(i,j,k)*(value_rhored-epijexp(i,j,k))-beijexp(i,j,k))
                        );
            }
        }
    }
    double a01=0;
    double a0=0;
    double cv1=0;
    for (uword i : m_indices)
    {
       a01=a01+X(i)*(Tc(i)/tred)*(8.314510/8.314472)
               *(
                   noik(i, 2-1)+noik(i, 3-1)*temperature/Tc(i)+noik(i, 4-1)*voik(i, 4-1)/tanh(voik(i, 4-1)*Tc(i)/temperature)+noik(i, 6-1)
                   *voik(i, 6-1)/tanh(voik(i, 6-1)*Tc(i)/temperature)-noik(i, 5-1)*voik(i, 5-1)*tanh(voik(i, 5-1)*Tc(i)/temperature)
                   -noik(i, 7-1)*voik(i, 7-1)*tanh(voik(i, 7-1)*Tc(i)/temperature)
                );
       //a0=a0+X(w)*((8.314510/8.314472)*(log(value/rhoc(w))+noik(w, 1-1)+noik(w, 2-1)*Tc(w)/T(i)+noik(w, 3-1)*log(Tc(w)/T(i))+noik(w, 4-1)*log(sinh(voik(w, 4-1)*Tc(w)/T(i)))+noik(w, 6-1)*log(sinh(voik(w, 6-1)*Tc(w)/T(i)))-noik(w, 5-1)*log(cosh(voik(w, 5-1)*Tc(w)/T(i)))-noik(w, 7-1)*log(cosh(voik(w, 7-1)*Tc(w)/T(i))))+log(X(w)));
//       a0=a0+X(i)*
//                (
//                   log(value/rhoc(i))+(8.314510/8.314472)
//                   *(
//                       noik(i, 1-1)+noik(i, 2-1)*Tc(i)/temperature+noik(i, 3-1)*log(Tc(i)/temperature)+noik(i, 4-1)*log(abs(sinh(voik(i, 4-1)*Tc(i)/temperature)))
//                       +noik(i, 6-1)*log(abs(sinh(voik(i, 6-1)*Tc(i)/temperature)))-noik(i, 5-1)*log(cosh(voik(i, 5-1)*Tc(i)/temperature))
//                       -noik(i, 7-1)*log(cosh(voik(i, 7-1)*Tc(i)/temperature))
//                    )+log(X(i))
//                );
       a0 += X(i)*
           (
               // a^o_oi
               (8.314510/8.314472)
               *(
                   log(value/rhoc(i))
                   +noik(i, 1-1)
                   +noik(i, 2-1)*Tc(i)/temperature
                   +noik(i, 3-1)*log(Tc(i)/temperature)
                   +noik(i, 4-1)*log(std::abs(sinh(voik(i, 4-1)*Tc(i)/temperature)))
                   +noik(i, 6-1)*log(std::abs(sinh(voik(i, 6-1)*Tc(i)/temperature)))
                   -noik(i, 5-1)*log(cosh(voik(i, 5-1)*Tc(i)/temperature))
                   -noik(i, 7-1)*log(cosh(voik(i, 7-1)*Tc(i)/temperature))
                )
                +log(X(i))
            );
       cv1=cv1+X(i)*pow2(Tc(i)/tred)*(8.314510/8.314472)
               *(
                   -noik(i, 3-1)*pow2(temperature/Tc(i))-noik(i, 4-1)*pow2(voik(i, 4-1))/pow2(sinh(voik(i, 4-1)*Tc(i)/temperature))
                   -noik(i, 6-1)*pow2(voik(i, 6-1))/pow2(sinh(voik(i, 6-1)*Tc(i)/temperature))-noik(i, 5-1)*pow5(voik(i, 5-1))/pow2(cosh(voik(i, 5-1)*Tc(i)/temperature))
                   -noik(i, 7-1)*pow2(voik(i, 7-1))/pow2(cosh(voik(i, 7-1)*Tc(i)/temperature))
                );
    }

    double ar1=0;
    double ar2=0;
    double cv2=0;
    double cp11=0;
    double cp111=0;
    for (uword i : m_indices)
    {
        for (uword k = 0; k < 24; k++)
        {
            const double pow_value_rhored_doipol = pow(value_rhored, doipol(i,k));
            const double pow_tred_temperature_toipol = pow(tred_temperature, toipol(i,k));
            const double pow_value_rhored_doiexp = pow(value_rhored, doiexp(i,k));
            const double pow_value_rhored_coiexp = pow(value_rhored, coiexp(i,k));
            const double pow_tred_temperature_toiexp = pow(tred_temperature, toiexp(i,k));
            ar1 +=
                    X(i)
                    *(
                        noipol(i,k)*toipol(i,k)*pow_value_rhored_doipol*pow_tred_temperature_toipol*tred_temperatured_pow_minusOne
                        +noiexp(i,k)*toiexp(i,k)*pow_value_rhored_doiexp*pow_tred_temperature_toiexp*tred_temperatured_pow_minusOne*exp(-pow_value_rhored_coiexp)
                    );
            ar2 +=
                    X(i)
                    *(
                        noipol(i,k)*pow_value_rhored_doipol*pow_tred_temperature_toipol
                        +noiexp(i,k)*pow_value_rhored_doiexp*pow_tred_temperature_toiexp*exp(-pow_value_rhored_coiexp)
                    );
            cv2 +=
                    +X(i)
                    *(
                        noipol(i,k)*toipol(i,k)*(toipol(i,k)-1)*pow_value_rhored_doipol*pow_tred_temperature_toipol*tred_temperatured_pow_minusTwo
                        +noiexp(i,k)*toiexp(i,k)*(toiexp(i,k)-1)*pow_value_rhored_doiexp*pow_tred_temperature_toiexp*tred_temperatured_pow_minusTwo
                            *exp(-pow_value_rhored_coiexp)
                    );
            cp11 +=
                    +X(i)
                    *(
                        noipol(i,k)*doipol(i,k)*toipol(i,k)*pow_value_rhored_doipol*value_rhored_pow_minusOne*pow_tred_temperature_toipol*tred_temperatured_pow_minusOne
                        +noiexp(i,k)*toiexp(i,k)*pow_value_rhored_doiexp*value_rhored_pow_minusOne*(doiexp(i,k)-coiexp(i,k)*pow_value_rhored_coiexp)
                            *pow_tred_temperature_toiexp*tred_temperatured_pow_minusOne*exp(-pow_value_rhored_coiexp)
                    );
            cp111 +=
                    +X(i)
                    *(
                        noipol(i,k)*doipol(i,k)*(doipol(i,k)-1)*pow_value_rhored_doipol*value_rhored_pow_minusTwo*pow_tred_temperature_toipol
                        +noiexp(i,k)*pow_value_rhored_doiexp*value_rhored_pow_minusTwo
                            *(
                                (doiexp(i,k)-coiexp(i,k)*pow_value_rhored_coiexp)*(doiexp(i,k)-1-coiexp(i,k)*pow_value_rhored_coiexp)
                                -pow2(coiexp(i,k))*pow_value_rhored_coiexp
                            )*pow_tred_temperature_toiexp*exp(-pow_value_rhored_coiexp)
                    );
        }
    }

    double ar11=0;
    double ar22=0;
    double cv3=0;
    double cp22=0;
    double cp222=0;
    for (uword i : m_firstIndices)
    {
        for (uword j : m_lastIndices)
        {
            for (uword k = 0; k < 20; k++)
            {
                const double pow_value_rhored_dijpol = pow(value_rhored, dijpol(i,j,k));
                const double pow_value_rhored_dijexp = pow(value_rhored, dijexp(i,j,k));

                const double pow_tred_temperature_tijpol = pow(tred_temperature, tijpol(i,j,k));
                const double pow_tred_temperature_tijexp = pow(tred_temperature, tijexp(i,j,k));
                const double pow_tred_temperature_tijexp_minusOne = pow_tred_temperature_tijexp*tred_temperatured_pow_minusOne;

                const double exp_nuijexp_pow2_value_rhored_epijexp_beijexp_value_rhored_gaijexp = exp(-nuijexp(i,j,k)*pow2(value_rhored-epijexp(i,j,k))-beijexp(i,j,k)*(value_rhored-gaijexp(i,j,k)));
                const double XiXjFij = X(i)*X(j)*Fij(i,j);

                ar11 += XiXjFij*
                        (
                            nijpol_times_tijpol(i,j,k)*pow_value_rhored_dijpol*pow_tred_temperature_tijpol*tred_temperatured_pow_minusOne
                            + nijexp_times_tijexp(i,j,k)*pow_value_rhored_dijexp*pow_tred_temperature_tijexp_minusOne
                            *exp_nuijexp_pow2_value_rhored_epijexp_beijexp_value_rhored_gaijexp
                        );
                ar22 += XiXjFij*
                        (
                            nijpol(i,j,k)*pow_value_rhored_dijpol*pow_tred_temperature_tijpol
                            +nijexp(i,j,k)*pow_value_rhored_dijexp*pow_tred_temperature_tijexp
                            *exp_nuijexp_pow2_value_rhored_epijexp_beijexp_value_rhored_gaijexp
                        );
                cv3 += XiXjFij*
                        (
                            nijpol_times_tijpol_times_tijpol_minus_one(i,j,k)*pow_value_rhored_dijpol*pow_tred_temperature_tijpol*tred_temperatured_pow_minusTwo
                            +nijexp_times_tijexp_times_tijexp_minus_one(i,j,k)*pow_value_rhored_dijexp*pow_tred_temperature_tijexp*tred_temperatured_pow_minusTwo
                            *exp_nuijexp_pow2_value_rhored_epijexp_beijexp_value_rhored_gaijexp
                        );
                cp22 += XiXjFij*
                        (
                            +nijpol_times_dijpol_times_tijpol(i,j,k)*pow_value_rhored_dijpol*value_rhored_pow_minusOne
                            +nijexp_times_tijexp(i,j,k)*pow_value_rhored_dijexp*pow_tred_temperature_tijexp_minusOne
                                *exp_nuijexp_pow2_value_rhored_epijexp_beijexp_value_rhored_gaijexp
                                *(dijexp(i,j,k)/(value_rhored)-2*nuijexp(i,j,k)*(value_rhored-epijexp(i,j,k))-beijexp(i,j,k))
                        );
                cp222 += XiXjFij*
                        (
                            +nijpol_times_dijpol_times_dijpol_minus_one(i,j,k)*pow_value_rhored_dijpol*value_rhored_pow_minusTwo*pow_tred_temperature_tijpol
                            +nijexp(i,j,k)*pow_value_rhored_dijexp
                                *pow_tred_temperature_tijexp
                                *exp_nuijexp_pow2_value_rhored_epijexp_beijexp_value_rhored_gaijexp
                                *(pow2(dijexp(i,j,k)/(value_rhored)-2*nuijexp(i,j,k)*(value_rhored-epijexp(i,j,k))-beijexp(i,j,k))-dijexp(i,j,k)/pow2(value_rhored)-2*nuijexp(i,j,k))
                        );
            }
        }
    }

    vec output = zeros<vec>(15);

    double adeltar = aroidelta + arijdelta;
    double rho = value;
    double delta = rho*rhored; // rhored = 1/\rho_r(\bar x)
    double tau = tred_temperature;

    double ataur = ar1+ar11;

    //Z factor
    output(0)=
            1+delta*adeltar;
    //dzdT_p
    output(1)=
            (-delta*adeltar/temperature+delta*(dzdTp2+dzdTpij2)/temperature)
            /(1+delta*adeltar/output(0)+delta*(dzdTp+dzdTpij)/output(0));
    //dzdp_T
    output(2)=
            (adeltar*delta/pressure+delta*dzdp/pressure)
            /(1+delta*adeltar/output(0)+delta*dzdp/output(0));
    //dzdpT_rho
    output(3)=
            delta*(dzdT_rho+dzdT_rho2);
    //entropy s
    double ar = (ar2+ar22);
    output(4)=
            Ra*(tau*(a01+ataur) - a0 - ar);
    //internal energy
    output(5)=
            Ra*temperature*tau*(a01+ataur);
    //heat capacity constant volume cv
    output(6)=
            -Ra*pow2(tau)*(cv1+cv2+cv3);
    //enthalpy
    output(7)=
            Ra*temperature*(1 + tau*(a01 + ataur) + delta*adeltar);
    //heat capacity constant pressure cp
    output(8)=
            -Ra*pow2(tau)*(cv1+cv2+cv3)
            +Ra*pow2(1+delta*adeltar-delta*tau*(cp11+cp22))
            /(1+2*delta*adeltar+pow2(delta)*(cp111+cp222));
    //gibbs free energi
    output(9)=
            Ra*temperature*(1+a0+ar1+ar11+delta*adeltar); // THIS MIGHT BE WRONG, should be a^r, not a_tau^r (== ar1+ar11)
    //joule thomson coefficient
    output(10)=
            -(1.0/(Ra*rho))*(delta*adeltar+pow2(delta)*(cp111+cp222)+delta*tau*(cp11+cp22))
            /(
                pow2(1+delta*adeltar-delta*tau*(cp11+cp22))
                -pow2(tau)*(cv1+cv2+cv3)*(1+2*delta*adeltar+pow2(delta)*(cp111+cp222))
            );
    //speed of sound
    output(11)=
            sqrt(
                Ra*temperature
                *(
                    1+2*delta*adeltar+pow2(delta)*(cp111+cp222)
                    -pow2(1+delta*adeltar-delta*tau*(cp11+cp22))
                    /(pow2(tau)*(cv1+cv2+cv3))
                )
            );
    //isothermal throttling coefficient
    output(12)=
            (
                1-(1+delta*adeltar-delta*tau*(cp11+cp22))
                /(1+2*delta*adeltar+pow2(delta)*(cp111+cp222))
            )/rho;
    //density
    output(13)=
            rho;
    //isentropic exponent
    output(14)=
            (1+2*delta*adeltar+pow2(delta)*(cp111+cp222))/(1+delta*adeltar)
            *(
                1
                -pow2(1+delta*adeltar-delta*tau*(cp11+cp22))
                /(pow2(tau)*(cv1+cv2+cv3)*(1+2*delta*adeltar+pow2(delta)*(cp111+cp222)))
            );

    return output;
}

double GERG04::findSpeedOfSound(const double temperature, const double density) const
{
    // find alpha coefficient
    double value = density;
    double aroidelta = 0;
    double arijdelta = 0;
    double aroideltadelta = 0;
    double arijdeltadelta = 0;

    double start = density;
//    double pow2_start = start*start;
    double start_rhored = start*rhored;
    double start_rhored_pow_minusOne = 1.0/start_rhored;
    double start_rhored_pow_minusTwo = 1.0/pow2(start_rhored);

    const double value_rhored = value*rhored;
    const double tred_temperature = tred/temperature;
    const double value_rhored_pow_minusOne = pow(value_rhored,-1);
    const double value_rhored_pow_minusTwo = pow(value_rhored,-2);
    const double tred_temperatured_pow_minusOne = pow(tred_temperature,-1);
    const double tred_temperatured_pow_minusTwo = pow(tred_temperature,-2);

    evaluateAlpha_roi_deltas(tred_temperature, start_rhored, start_rhored_pow_minusOne, start_rhored_pow_minusTwo, aroidelta, arijdelta, aroideltadelta, arijdeltadelta);

    double cv1=0;
    for (uword i : m_indices)
    {
       cv1=cv1+X(i)*pow2(Tc(i)/tred)*(8.314510/8.314472)
               *(
                   -noik(i, 3-1)*pow2(temperature/Tc(i))-noik(i, 4-1)*pow2(voik(i, 4-1))/pow2(sinh(voik(i, 4-1)*Tc(i)/temperature))
                   -noik(i, 6-1)*pow2(voik(i, 6-1))/pow2(sinh(voik(i, 6-1)*Tc(i)/temperature))-noik(i, 5-1)*pow5(voik(i, 5-1))/pow2(cosh(voik(i, 5-1)*Tc(i)/temperature))
                   -noik(i, 7-1)*pow2(voik(i, 7-1))/pow2(cosh(voik(i, 7-1)*Tc(i)/temperature))
                );
    }

    double cv2=0;
    double cp11=0;
    double cp111=0;
    for (uword i : m_indices)
    {
        for (uword k = 0; k < 24; k++)
        {
            const double pow_value_rhored_doipol = pow(value_rhored, doipol(i,k));
            const double pow_tred_temperature_toipol = pow(tred_temperature, toipol(i,k));
            const double pow_value_rhored_doiexp = pow(value_rhored, doiexp(i,k));
            const double pow_value_rhored_coiexp = pow(value_rhored, coiexp(i,k));
            const double pow_tred_temperature_toiexp = pow(tred_temperature, toiexp(i,k));

            cv2 +=
                    +X(i)
                    *(
                        noipol(i,k)*toipol(i,k)*(toipol(i,k)-1)*pow_value_rhored_doipol*pow_tred_temperature_toipol*tred_temperatured_pow_minusTwo
                        +noiexp(i,k)*toiexp(i,k)*(toiexp(i,k)-1)*pow_value_rhored_doiexp*pow_tred_temperature_toiexp*tred_temperatured_pow_minusTwo
                            *exp(-pow_value_rhored_coiexp)
                    );
            cp11 +=
                    +X(i)
                    *(
                        noipol(i,k)*doipol(i,k)*toipol(i,k)*pow_value_rhored_doipol*value_rhored_pow_minusOne*pow_tred_temperature_toipol*tred_temperatured_pow_minusOne
                        +noiexp(i,k)*toiexp(i,k)*pow_value_rhored_doiexp*value_rhored_pow_minusOne*(doiexp(i,k)-coiexp(i,k)*pow_value_rhored_coiexp)
                            *pow_tred_temperature_toiexp*tred_temperatured_pow_minusOne*exp(-pow_value_rhored_coiexp)
                    );
            cp111 +=
                    +X(i)
                    *(
                        noipol(i,k)*doipol(i,k)*(doipol(i,k)-1)*pow_value_rhored_doipol*value_rhored_pow_minusTwo*pow_tred_temperature_toipol
                        +noiexp(i,k)*pow_value_rhored_doiexp*value_rhored_pow_minusTwo
                            *(
                                (doiexp(i,k)-coiexp(i,k)*pow_value_rhored_coiexp)*(doiexp(i,k)-1-coiexp(i,k)*pow_value_rhored_coiexp)
                                -pow2(coiexp(i,k))*pow_value_rhored_coiexp
                            )*pow_tred_temperature_toiexp*exp(-pow_value_rhored_coiexp)
                    );
        }
    }

    double ar22=0;
    double cv3=0;
    double cp22=0;
    double cp222=0;
    for (uword i : m_firstIndices)
    {
        for (uword j : m_lastIndices)
        {
            for (uword k = 0; k < 20; k++)
            {
                const double pow_value_rhored_dijpol = pow(value_rhored, dijpol(i,j,k));
                const double pow_value_rhored_dijexp = pow(value_rhored, dijexp(i,j,k));

                const double pow_tred_temperature_tijpol = pow(tred_temperature, tijpol(i,j,k));
                const double pow_tred_temperature_tijexp = pow(tred_temperature, tijexp(i,j,k));
                const double pow_tred_temperature_tijexp_minusOne = pow_tred_temperature_tijexp*tred_temperatured_pow_minusOne;

                const double exp_nuijexp_pow2_value_rhored_epijexp_beijexp_value_rhored_gaijexp = exp(-nuijexp(i,j,k)*pow2(value_rhored-epijexp(i,j,k))-beijexp(i,j,k)*(value_rhored-gaijexp(i,j,k)));
                const double XiXjFij = X(i)*X(j)*Fij(i,j);

                ar22 += XiXjFij*
                        (
                            nijpol(i,j,k)*pow_value_rhored_dijpol*pow_tred_temperature_tijpol
                            +nijexp(i,j,k)*pow_value_rhored_dijexp*pow_tred_temperature_tijexp
                            *exp_nuijexp_pow2_value_rhored_epijexp_beijexp_value_rhored_gaijexp
                        );
                cv3 += XiXjFij*
                        (
                            nijpol_times_tijpol_times_tijpol_minus_one(i,j,k)*pow_value_rhored_dijpol*pow_tred_temperature_tijpol*tred_temperatured_pow_minusTwo
                            +nijexp_times_tijexp_times_tijexp_minus_one(i,j,k)*pow_value_rhored_dijexp*pow_tred_temperature_tijexp*tred_temperatured_pow_minusTwo
                            *exp_nuijexp_pow2_value_rhored_epijexp_beijexp_value_rhored_gaijexp
                        );
                cp22 += XiXjFij*
                        (
                            +nijpol_times_dijpol_times_tijpol(i,j,k)*pow_value_rhored_dijpol*value_rhored_pow_minusOne
                            +nijexp_times_tijexp(i,j,k)*pow_value_rhored_dijexp*pow_tred_temperature_tijexp_minusOne
                                *exp_nuijexp_pow2_value_rhored_epijexp_beijexp_value_rhored_gaijexp
                                *(dijexp(i,j,k)/(value_rhored)-2*nuijexp(i,j,k)*(value_rhored-epijexp(i,j,k))-beijexp(i,j,k))
                        );
                cp222 += XiXjFij*
                        (
                            +nijpol_times_dijpol_times_dijpol_minus_one(i,j,k)*pow_value_rhored_dijpol*value_rhored_pow_minusTwo*pow_tred_temperature_tijpol
                            +nijexp(i,j,k)*pow_value_rhored_dijexp
                                *pow_tred_temperature_tijexp
                                *exp_nuijexp_pow2_value_rhored_epijexp_beijexp_value_rhored_gaijexp
                                *(pow2(dijexp(i,j,k)/(value_rhored)-2*nuijexp(i,j,k)*(value_rhored-epijexp(i,j,k))-beijexp(i,j,k))-dijexp(i,j,k)/pow2(value_rhored)-2*nuijexp(i,j,k))
                        );
            }
        }
    }

    double adeltar = aroidelta + arijdelta;
    double rho = value;
    double delta = rho*rhored; // rhored = 1/\rho_r(\bar x)
    double tau = tred_temperature;

    //speed of sound
    double speedOfSound =
            sqrt(
                Ra*temperature
                *(
                    1+2*delta*adeltar+pow2(delta)*(cp111+cp222)
                    -pow2(1+delta*adeltar-delta*tau*(cp11+cp22))
                    /(pow2(tau)*(cv1+cv2+cv3))
                )
            );

    return speedOfSound;
}
