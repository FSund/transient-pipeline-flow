#pragma once

#include <string>
#include <cmath>
#include <exception>
#include <armadillo>

/*!
  * \defgroup typedefs Typedefs
  */

namespace utils
{
    inline double pow2(const double a);
    inline arma::vec pow2(const arma::vec& a);
    inline double pow3(const double a);
    inline double pow4(const double a);
    inline double pow5(const double a);
    inline double pow6(const double a);

    inline arma::vec cubeRoot(const arma::vec& input);

    inline arma::vec centerDifference(const arma::vec& input);

    inline arma::vec centerDifference(const arma::subview_col<arma::mat::elem_type>& input);

    inline arma::vec centerSum(const arma::vec& input);

    inline arma::vec centerSum(const arma::subview_col<arma::mat::elem_type>& input);

    inline arma::vec centerAverage(const arma::vec& input);

    inline arma::vec centerAverage(const arma::subview_col<arma::mat::elem_type>& input);

    inline arma::vec weightedAverage(const arma::vec& weightA, const arma::vec& weightB, const arma::vec& A, const arma::vec& B);

    inline double weightedAverage(const double weightA, const double weightB, const double A, const double B);

    inline bool isRelativeDifferenceWithinTolerance(
            const double oldValue,
            const double newValue,
            const double tolerance);

    inline bool areRelativeDifferencesWithinTolerance(
            const arma::vec& oldValues,
            const arma::vec& newValues,
            const double tolerance);

    inline double convertBurialToStandardDefinition(
            const double burial,
            const double innerDiameter,
            const arma::rowvec& width);

    inline double convertBurialToStandardDefinition(
            const double burial,
            const double innerDiameter,
            const arma::subview_row<arma::mat::elem_type>& width);

    inline arma::vec convertBurialToStandardDefinition(
            const arma::vec& burial,
            const arma::vec& innerDiameter,
            const arma::colvec& width);

    inline arma::vec convertBurialToStandardDefinition(
            const arma::vec& burial,
            const arma::vec& innerDiameter,
            const arma::mat& width);

    inline arma::vec smooth(const arma::vec& x, const double smoothness);

    inline arma::vec smoothTransition(const arma::vec& firstPart, const arma::vec& lastPart, double smoothness = 1.0, const double timeStep = 60);

    inline arma::vec createSmoothTransient(
            const double initialValue,
            const double finalValue,
            const arma::uword n,
            const double smoothness,
            const arma::uword dt);

    inline arma::vec findLogSpacedConcentricShellWidths(
            const double innerRadius,
            const double outerRadius,
            const arma::uword nShells);

    inline arma::cube loadCubeFromFile(const std::string& filename, const arma::file_type& fileType = arma::csv_ascii);
    inline arma::mat loadMatFromFile(const std::string& filename, const arma::file_type& fileType = arma::csv_ascii);
    inline arma::vec loadVecFromFile(const std::string& filename, const arma::file_type& fileType = arma::csv_ascii);
}

inline double utils::pow2(const double a)
{
    return a*a;
}

inline arma::vec utils::pow2(const arma::vec& a)
{
    return a%a;
}

inline double utils::pow3(const double a)
{
    return a*a*a;
}

inline double utils::pow4(const double a)
{
    return a*a*a*a;
}

inline double utils::pow5(const double a)
{
    return a*a*a*a*a;
}

inline double utils::pow6(const double a)
{
    return a*a*a*a*a*a;
}

inline arma::vec utils::cubeRoot(const arma::vec& input)
{
    arma::vec cubed = arma::zeros<arma::vec>(input.n_elem);
    for (arma::uword i = 0; i < input.n_elem; i++)
    {
        cubed(i) = std::cbrt(input(i));
    }
    return cubed;
}

arma::vec utils::centerDifference(const arma::vec& input)
{
    return input(arma::span(1, input.n_elem - 1)) - input(arma::span(0, input.n_elem - 2));
}

arma::vec utils::centerDifference(const arma::subview_col<arma::mat::elem_type>& input)
{
    return input.rows(1, input.n_elem - 1) - input.rows(0, input.n_elem - 2);
}

arma::vec utils::centerSum(const arma::vec& input)
{
    return input(arma::span(1, input.n_elem - 1)) + input(arma::span(0, input.n_elem - 2));
}

arma::vec utils::centerSum(const arma::subview_col<arma::mat::elem_type>& input)
{
    return input.rows(1, input.n_elem - 1) + input.rows(0, input.n_elem - 2);
}

arma::vec utils::centerAverage(const arma::vec& input)
{
    return centerSum(input)/2.0;
}

arma::vec utils::centerAverage(const arma::subview_col<arma::mat::elem_type>& input)
{
    return centerSum(input)/2.0;
}

double utils::weightedAverage(const double weightA, const double weightB, const double A, const double B)
{
    if (std::round(weightA + weightB) == 0)
    {
        throw std::invalid_argument("utils::weightedAverage(): both weights are zero");
    }

    return (A*weightA + B*weightB)/(weightA + weightB);
}

arma::vec utils::weightedAverage(const arma::vec& weightA, const arma::vec& weightB, const arma::vec& A, const arma::vec& B)
{
    if (!arma::all((weightA + weightB) != 0))
    {
        throw std::invalid_argument("utils::weightedAverage(): both weights are zero");
    }

    return (A % weightA + B % weightB)/(weightA + weightB);
}

bool utils::isRelativeDifferenceWithinTolerance(
        const double oldValue,
        const double newValue,
        const double tolerance)
{
    double relativeDifference;
    if (oldValue != 0)
    {
        relativeDifference = std::abs(newValue - oldValue)/std::abs(oldValue);
    }
    else if (newValue != 0)
    {
        relativeDifference = std::abs(newValue - oldValue)/std::abs(newValue);
    }
    else
    {
        // both values zero, difference also zero
        return true;
    }
    return relativeDifference < tolerance;
}

bool utils::areRelativeDifferencesWithinTolerance(
        const arma::vec& oldValues,
        const arma::vec& newValues,
        const double tolerance)
{
    for (arma::uword i = 0; i < oldValues.n_elem; i++)
    {
        double relativeDifference;
        if (oldValues(i) != 0)
        {
            relativeDifference = std::abs(newValues(i) - oldValues(i))/std::abs(oldValues(i));
        }
        else if (newValues(i) != 0)
        {
            relativeDifference = std::abs(newValues(i) - oldValues(i))/std::abs(newValues(i));
        }
        else
        {
            // both values zero, difference also zero
            continue;
        }

        if (relativeDifference > tolerance)
        {
            return false;
        }
    }

    return true;
}

double utils::convertBurialToStandardDefinition(
        const double burial,
        const double innerDiameter,
        const arma::rowvec& width)
{
    return burial + innerDiameter/2.0 + arma::sum(width);
}

double utils::convertBurialToStandardDefinition(
        const double burial,
        const double innerDiameter,
        const arma::subview_row<arma::mat::elem_type>& width)
{
    return burial + innerDiameter/2.0 + arma::sum(width);
}

arma::vec utils::convertBurialToStandardDefinition(
        const arma::vec& burial,
        const arma::vec& innerDiameter,
        const arma::colvec& width)
{
    return burial + innerDiameter/2.0 + arma::sum(width);
}

arma::vec utils::convertBurialToStandardDefinition(
        const arma::vec& burial,
        const arma::vec& innerDiameter,
        const arma::mat& width)
{
    // burial depth is defined as "from top of soil to center of pipe" in our code,
    // but given as distance from top of pipe to top of soil in thesis and in SIM
    // this means that we have to add half a diameter + the width of all layers to convert
    arma::vec convertedBurial = arma::zeros<arma::vec>(burial.n_elem);
    for (arma::uword i = 0; i < burial.n_elem; i++)
    {
        convertedBurial(i) = convertBurialToStandardDefinition(burial(i), innerDiameter(i), width.row(i));
    }

    return convertedBurial;
}

arma::vec utils::smooth(const arma::vec& x, const double smoothness)
{
    return 0.5 + 0.5*tanh((x - 0.5)/smoothness); // use x - 0.5 to get transition point to be between x = 0 and x = 1.0
}

arma::vec utils::smoothTransition(
        const arma::vec& firstPart,
        const arma::vec& lastPart,
        double smoothness,
        const double timeStep)
{
    // smooths transition from firstPart to lastPart by using a tanh function
    // idea from here: http://www.j-raedler.de/2010/10/smooth-transition-between-functions-with-tanh/

    // the smoothing is time-independent (or time-depending, depending on how you look at it),
    // meaning that changing the time step but using the same smoothness should result in the same smoothed curve,
    // if you also scale the x-axis correctly (meaning: use time on the x-axis, not number of time steps)

    // a smoothness of 1.0 will give a transition that uses approx. 10 minutes
    // a smoothness of 5.0 will use around 30 minutes

    arma::uword m = firstPart.n_elem;
    arma::uword n = lastPart.n_elem;

    // non-dimensionalizing smoothness
    // default time step of 60 gives no scaling
    smoothness *= 60.0/timeStep; // 60 is standard timestep

    arma::vec f = arma::join_cols(
                firstPart,
                arma::zeros<arma::vec>(n) + firstPart(m-1) // repeat last element of firstPart
                );
    arma::vec g = arma::join_cols(
                arma::zeros<arma::vec>(m) + lastPart(0), // repeat first element of lastPart
                lastPart
                );

    arma::vec x = arma::linspace(-int(m)+1, n, m+n);  // !!! cast to int before subtracting, since m is unsigned int!!!
    arma::vec s = smooth(x, smoothness);

    arma::vec h = (1-s) % f + s % g;

    return h;
}

arma::vec utils::createSmoothTransient(
        const double initialValue,
        const double finalValue,
        const arma::uword n,
        const double smoothness,
        const arma::uword dt)
{
    // transient will be centered in interval [0, n]
    if (n % 2 != 0)
    {
        std::cout << "n (" << n << ") should be divisible by 2!" << std::endl;
    }
    arma::vec a = arma::zeros<arma::vec>(n/2) + initialValue;
    arma::vec b = arma::zeros<arma::vec>(n/2) + finalValue;
    return smoothTransition(a, b, smoothness, dt);
}

arma::vec utils::findLogSpacedConcentricShellWidths(
        const double innerRadius,
        const double outerRadius,
        const arma::uword nShells)
{
    arma::vec logSpacedLayerRadii = arma::logspace(log10(innerRadius), log10(outerRadius), nShells+1);
    arma::vec logSpacedLayerWidths = arma::zeros<arma::vec>(nShells);
    for (arma::uword i = 0; i < nShells; i++)
    {
        logSpacedLayerWidths(i) = logSpacedLayerRadii(i+1) - logSpacedLayerRadii(i);
    }

    return logSpacedLayerWidths;
}

arma::cube utils::loadCubeFromFile(const std::string& filename, const arma::file_type& fileType)
{
    arma::cube c;
    if (!c.load(filename, fileType))
    {
        throw std::runtime_error("Couldn't load file \"" + filename + "\"");
    }

    return c;
}

arma::mat utils::loadMatFromFile(const std::string& filename, const arma::file_type& fileType)
{
    arma::mat m;
    if (!m.load(filename, fileType))
    {
        throw std::runtime_error("Couldn't load file \"" + filename + "\"");
    }

    return m;
}

arma::vec utils::loadVecFromFile(const std::string& filename, const arma::file_type& fileType)
{
    arma::vec v;
    if (!v.load(filename, fileType))
    {
        throw std::runtime_error("Couldn't load file \"" + filename + "\"");
    }

    return v;
}
