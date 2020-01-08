#pragma once

#include <armadillo>

class AmbientFluid;

/*!
 * Namespace for some utilities.
 */
namespace utils {

    /*!
     * \brief Calculate thermal conductivity of natural gas at given pressure.
     *
     * This is some kind of empirical relation or fitting, with an unknown source.
     * It should probably be replaced with something else. See for example
     * <a href="https://doi.org/10.1016/j.jngse.2014.04.005"><i>A simple correlation to estimate natural gas thermal conductivity</i> (Azad Jarrahiana and Ehsan Heidaryan, Journal of Natural Gas Science and Engineering, Volume 18, May 2014)</a>.
     *
     * \param pressure Gas pressure [Pa].
     * \return Thermal conductivity [W/(m K)]
     */
    double calcGasThermalConductivity(const double pressure);

    /*!
     * \brief Calculate outer film coefficient for a given outer diameter and
     * AmbientFluid.
     *
     * This is just a wrapper around utils::calcOuterWallFilmCoefficient(const double, const double, const double, const double, const double, const double)
     *
     * \param diameter Outer diameter [m]
     * \param fluid AmbientFluid describing the fluid
     * \return Outer film coefficient [W/(m2 K)]
     */
    double calcOuterWallFilmCoefficient(
            const double diameter,
            const AmbientFluid& fluid
            );

    /*!
     * \brief Calculate the outer film coefficient for external flow normal to a
     * circular cylinder.
     *
     * Uses eq. 7.52 from Fundamentals of heat and mass transfer (7th Ed, 2011) (Bergman, Lavine, Incropera, DeWitt).
     *
     * \param diameter Outer diameter [m]
     * \param heatCapacityConstantPressure Fluid heat capacity (\f$c_p\f$) [J/(kg K)]
     * \param viscosity Fluid dynamic viscosity [Pa s] = [kg/m*s]
     * \param thermalConductivity Fluid thermal conductivity [W/(m K)]
     * \param density Fluid density [kg/m3]
     * \param velocity Fluid velocity [m/s]
     * \return Outer film coefficient [W/(m2 K)]
     */
    double calcOuterWallFilmCoefficient(
            const double diameter,
            const double heatCapacityConstantPressure = 4200, // [J/kg K]
            const double viscosity = 1.05/1000.0, // [Pa s] = [kg/m*s]
            const double thermalConductivity = 0.57, // [W/m K]
            const double density = 1020, // [kg/m3]
            const double velocity = 0.1 // [m/s]
            );

    /*!
     * \brief Calculate inner wall film coefficient for flow inside a cylinder.
     *
     * This uses the Dittus-Boelter equation at Reynolds numbers abve 1e4,
     * eq. 8.55 Fundamentals of heat and mass transfer (7th Ed, 2011) (Bergman, Lavine, Incropera, DeWitt)
     * at Reynolds number between 1e4 and 4e4, and returns 0 below 4e4.
     *
     * \param diameter Inner diameter [m]
     * \param fluidPressure Fluid pressure [Pa]
     * \param fluidReynoldsNumber Fluid Reynolds number [-]
     * \param fluidHeatCapacityConstantPressure Fluid heat capacity (\f$c_p\f$) [J/(kg K)]
     * \param fluidViscosity Fluid dynamic viscosity [Pa s] = [kg/m*s]
     * \return Inner film coefficient [W/(m2 K)]
     */
    double calcInnerWallFilmCoefficient(
            const double diameter,
            const double fluidPressure,
            const double fluidReynoldsNumber,
            const double fluidHeatCapacityConstantPressure,
            const double fluidViscosity);

    /*!
     * \brief Calculate equivalent burial layer thickness.
     *
     * For a pipeline buried in a medium at a given depth, this function
     * calculates the thickness of an equivalent cylinder shell of the same
     * medium around the pipeline which gives the same heat transfer between
     * the fluid in the pipeline and the ambient medium. The calculation is
     * based on equations in the documentation of the OLGA simulation software.
     *
     * Using the thickness from this function ensures the same results are
     * achieved with steady state and unsteady heat transfer models.
     *
     * \see utils::calcEquivalentBurialLayerRadius()
     *
     * \param innerDiameter Pipeline inner diameter [m]
     * \param wallThickness Pipeline wall thickness [m]
     * \param burialDepth Distance from top of pipe to top of burial medium [m]
     * \param burialMediumConductivity Thermal conductivity of burial medium [W/(m K)]
     * \return
     */
    double calcEquivalentBurialLayerWidth(
            const double innerDiameter,
            const double wallThickness,
            const double burialDepth,
            const double burialMediumConductivity = 2.0);

    /*!
     * \brief calcEquivalentBurialLayerRadius
     *
     * For a pipeline buried in a medium at a given depth, this function
     * calculates the outer radius of an equivalent cylinder shell of the same
     * medium around the pipeline which gives the same heat transfer between
     * the fluid in the pipeline and the ambient medium. The calculation is
     * based on equations in the documentation of the OLGA simulation software.
     *
     * \param innerDiameter Pipeline inner diameter [m]
     * \param wallThickness Pipeline wall thickness [m]
     * \param burialDepth Distance from top of pipe to top of burial medium [m]
     * \param burialMediumConductivity Thermal conductivity of burial medium [W/(m K)]
     * \return Outer radius of equivalent cylinder shell of burial medium [m]
     */
    double calcEquivalentBurialLayerRadius(
            const double innerDiameter,
            const double wallThickness,
            const double burialDepth,
            const double burialMediumConductivity = 2.0);

    /*!
     * \brief Calculate logarithmically (log10) spaced cylinder shell widths.
     *
     * This function calculates the widths of nShells cylinder shells,
     * logarithmically (log10) spaced between innerRadius and outerRadius.
     *
     * \param innerRadius Inner radius [m]
     * \param outerRadius Outer radius [m]
     * \param nShells Number of cylinder shells
     * \return arma::vec of shell widths [m]
     */
    arma::vec calcLogSpacedShellWidths(
            const double innerRadius,
            const double outerRadius,
            const arma::uword nShells = 10);

    /*!
     * \brief Calculate the widths of equivalent burial cylinder shells.
     *
     * This is just a wrapper around utils::calcEquivalentBurialLayerRadius()
     * and utils::calcLogSpacedShellWidths() that divides the shell into
     * several logarithmically spaced shells.
     *
     * \see utils::calcEquivalentBurialLayerRadius()
     * \see utils::calcLogSpacedShellWidths()
     *
     * \param innerDiameter Inner diameter [m]
     * \param wallThickness Wall thickness [m]
     * \param burialDepth Distance from top of pipe to top of burial medium [m]
     * \param burialMediumConductivity Thermal conductivity of burial medium [W/(m K)]
     * \param nShells Number of shells
     * \return Widths of logarithmically spaced cylinder shells.
     */
    arma::vec calcEquivalentBurialLayerWidths(
            const double innerDiameter,
            const double wallThickness,
            const double burialDepth,
            const double burialMediumConductivity = 2.0,
            const arma::uword nShells = 10);
}
