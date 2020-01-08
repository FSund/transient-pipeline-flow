#include "equationofstate/equationofstatebase.hpp"

#include <cmath>
#include <limits>

#include "constants.hpp"

using arma::vec;

EquationOfStateBase::~EquationOfStateBase()
{}

EquationOfStateBase::EquationOfStateBase(const vec& composition)
{
    this->setComposition(composition);
}

vec EquationOfStateBase::evaluate(const double pressure, const double temperature, const vec& composition)
{
    this->setComposition(composition, false);
    return this->evaluate(pressure, temperature);
}

double EquationOfStateBase::calculateStandardDensity() const
{
    // [kg/Sm3] density of ideal gas at standard conditions
    const double standardIdealDensity =
            m_molarMassOfMixture*constants::standardPressure
            /(constants::gasConstant*constants::standardTemperature) // [g/Sm3]
            /1000.0; // convert to [kg/Sm3]

    const double compressibilityAtStandardConditions =
            this->calculateCompressibility(
                constants::standardPressure,
                constants::standardTemperature);

    // [kg/Sm3] density of real gas mixture
    const double standardDensity =
            standardIdealDensity/compressibilityAtStandardConditions;

    return standardDensity;
}

//double EquationOfStateBase::calculateCompressibility(const double pressure, const double temperature) const
//{
//    return this->calculateCompressibility(pressure, temperature);
//}

// TODO: reintroduce this function somewhere, but just as a way to convert from
// standard flow to mass flow
// maybe as part of physics.hpp?
//double EquationOfStateBase::findStandardDensity() const
//{
//    double standardIdealDensity = m_molarMassOfMixture*constants::standardPressure/(constants::gasConstant*constants::standardTemperature); // [g/Sm3] density of ideal gas at standard conditions
//    standardIdealDensity /= 1000.0; // convert to kg/Sm3

//    double compressibilityAtStandardConditions =
//            this->calculateCompressibility(
//                constants::standardPressure,
//                constants::standardTemperature);

//    double standardDensity = standardIdealDensity/compressibilityAtStandardConditions; // [kg/Sm3] density of real gas mixture

//    return standardDensity;
//}

bool EquationOfStateBase::setComposition(const vec& composition, const bool force)
{
    if (composition.n_elem != 10)
        throw std::invalid_argument("number of components not equal 10");

    const auto limit = 10*std::numeric_limits<double>::epsilon();
    if (force
            // check if composition change above limit
            || arma::sum(arma::abs(composition - m_composition)) > limit)
    {
        m_composition = composition;
        m_molarMassOfMixture = sum(composition % m_molarMass);

        return true; // composition changed
    }
    else
    {
        return false; // composition not changed
    }
}

double EquationOfStateBase::getMolarMassOfMixture() const
{
    return m_molarMassOfMixture;
}

const vec& EquationOfStateBase::getComposition() const
{
    return m_composition;
}
