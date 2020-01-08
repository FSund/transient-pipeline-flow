#pragma once

#include <cmath>

/*!
 * This namespace contains a lot of useful constants.
 */
namespace constants
{
    //! The ratio of a circle's circumference to its diameter.
    static const double pi = 4.0*std::atan(1.0);

    //! The zero point of the Kelvin scale [degrees C].
    static const double kelvin = 273.15;

    //! Standard atmoshperic pressure [Pa].
    static const double standardAtmosphere = 101325;

    //! Convert Kelvin to Rankine [Ra/K].
    static const double kelvinToRankine = 9.0/5.0;

    //! Standard temperature [K].
    //! From Gassco.no (https://www.gassco.no/en/our-activities/what-is-natural-gas/glossary/).
    static const double standardTemperature = kelvin + 15;

    //! Standard pressure [Pa].
    //! From Gassco.no (https://www.gassco.no/en/our-activities/what-is-natural-gas/glossary/).
    static const double standardPressure = 101325;

    //! Universal gas constant [J/(K mol)].
    static const double gasConstant = 8.3144598;            // [J/mol K] == [m3 Pa /mol K]

    //! Universal gas constant [BTU/(R lbmol)]
    static const double gasConstantBTU = 1.98588;

    //! Convert Pascal to psf (pound force per square foot) [psf/Pa].
    static const double pascalToPoundForcePerSquareFoot = 0.3048*0.3048/4.4482216152605; // 4.448.. is lbf to Newton, 0.3048 is foot to meter

    //! Convert foot pound-force to Joule
    static const double footPoundForce = 1.355817948; // joule
    //! Convert slugs to kilogram
    static const double slug = 14.5939; // kg
    //! Convert degrees Rankine to degrees Kelvin
    static const double rankine = 5.0/9.0; // Kelvin
    //! Convert ft.lbf/slug.R to J/(kg K) (the unit of specific heat capacity)
    static const double footPoundForcePerSlugRankine = footPoundForce/(slug*rankine); // J/(kg K)

    //! The molar mass of air [g/mol].
    //! Used for calculating specific gravity.
    static const double molarMassOfAir = 28.967; // [g/mol] reference for calculating specific gravity
}
