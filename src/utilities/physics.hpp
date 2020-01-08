#pragma once

#include <armadillo>

class Composition;

namespace utils
{
    double calculateHeatCapacityConstantVolumeJFH(
            const double pressure);

    double calculateHeatCapacityConstantVolumeTGNet(
            const double molarMass,
            const double pressure,
            const double temperature);

    double calculateHeatCapacityConstantPressureJFH(
            const double molarMass,
            const double pressure,
            const double temperature);

    double calculateHeatCapacityConstantPressureLangelandsvik(
            const double molarMass,
            const double pressure,
            const double temperature);

    double calculateHeatCapacityConstantPressureTGNet(
            const double molarMassOfMixture,
            const double pressure,
            const double temperature);

    double calculateHeatCapacityConstantPressureKIO(
            const Composition& comp,
            const double pressure,
            const double temperature);

    double calculateIsobaricHeatCapacityJKH(
            const Composition& comp,
            const double pressure,
            const double temperature,
            const double Z = 0);

    arma::vec calculateViscosity(
            const arma::vec& molarMass,
            const arma::vec& temperature,
            const arma::vec& density);

    arma::vec calculateReynoldsNumber(
            const arma::vec& massFlow,
            const arma::vec& diameter,
            const arma::vec& viscosity);

    arma::vec calculateColebrookWhiteFrictionFactor(
            const arma::vec& sandGrainEquivalentRoughness,
            const arma::vec& diameter,
            const arma::vec& reynoldsNumber);

    double calculateColebrookWhiteFrictionFactor(
            const double sandGrainEquivalentRoughness,
            const double diameter,
            const double reynoldsNumber);

    double calculateHaalandFrictionFactor(
            const double sandGrainEquivalentRoughness,
            const double diameter,
            const double reynoldsNumber);

    namespace details
    {
        double KIOidealGasCP(
                const double specificGravity,
                const double temperature);

        double KIOdimensionlessResidualCP(
                const double reducedPressure,
                const double reducedTemperature);

        double JKHdimensionlessCP(
                const double molarMass,
                const double H2S,
                const double CO2,
                const double N2,
                const double H2,
                const double H2O,
                const double pressure,
                const double temperature,
                const double compressibility);

        double JKHidealGasCP(
                const Composition& comp,
                const double specificGravity,
                const double temperature);

        double JKHidealGasCP(
                const double specificGravity,
                const double H2S,
                const double CO2,
                const double N2,
                const double H2,
                const double H2O,
                const double temperature);

        double calculateHeatCapacityConstantPressureKIO(
                const Composition& comp,
                const double H2S,
                const double pressure,
                const double temperature);

        double colebrookWhiteFrictionFactor(
                const double sandGrainEquivalentRoughness,
                const double diameter,
                const double reynoldsNumber);

        double colebrookWhite(
                const double f,
                const double sandGrainEquivalentRoughness,
                const double diameter,
                const double reynoldsNumber);

        double colebrookWhiteDerivative(
                const double f,
                const double sandGrainEquivalentRoughness,
                const double diameter,
                const double reynoldsNumber);
    }
}
