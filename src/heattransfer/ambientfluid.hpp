#pragma once

#include "heattransfer/material.hpp"

/*!
 * \brief The AmbientFluid class is a simple container class that defines the
 * conductivity, density, heat capacity, dynamic viscosity and velocity of a
 * fluid used by HeatTransfer.
 */
class AmbientFluid : public Material
{
public:
    /*!
     * \brief Basic constructor that requires all material properties.
     * \param velocity Fluid velocity [m/s]
     * \param viscosity Fluid dynamic viscosity [Pa s] = [kg/m*s]
     * \param conductivity Thermal conductivity [W/(m K)]
     * \param density Density [kg/m3]
     * \param heatCapacity Heat capacity at constant pressure (\f$c_p\f$) [J/(kg K)]
     */
    constexpr AmbientFluid(
            const double velocity,
            const double viscosity,
            const double conductivity,
            const double density,
            const double heatCapacity):
        Material(conductivity, density, heatCapacity),
        m_velocity(velocity),
        m_viscosity(viscosity)
    {}

    /*!
     * \brief Construct from Material, velocity and viscosity.
     * \param velocity Fluid velocity [m/s]
     * \param viscosity Fluid dynamic viscosity [Pa s] = [kg/m*s]
     * \param material Fluid Material
     */
    constexpr AmbientFluid(
            const double velocity,
            const double viscosity,
            const Material& material):
        Material(material),
        m_velocity(velocity),
        m_viscosity(viscosity)
    {}

    // read-only
    double velocity() const { return m_velocity; } //!< Velocity getter
    double viscosity() const { return m_viscosity; } //!< Viscosity getter

    // setters
//    double& viscosity() { return m_viscosity; }
//    double& velocity() { return m_velocity; }
//    double& conductivity() { return m_conductivity; }
//    double& density() { return m_density; }
//    double& heatCapacity() { return m_heatCapacity; }

    // initialized below
    static const AmbientFluid seawater; //!< Predefined seawater
    static const AmbientFluid air; //!< Predefined air

private:
    double m_velocity; //!< Fluid velocity [m/s]
    double m_viscosity; //!< Fluid dynamic viscosity [Pa s] = [kg/m*s]
};

constexpr inline AmbientFluid AmbientFluid::seawater = AmbientFluid(0.1, 1.05/1000.0, Material::seawater);
constexpr inline AmbientFluid AmbientFluid::air      = AmbientFluid(0.1, 15.11e-6,    Material::air); // 15.11e-6 // air viscosity at 20 C and 1 atm [Pa s] = [kg/m*s]
