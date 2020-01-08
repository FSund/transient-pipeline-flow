#pragma once

#include "heattransfer/material.hpp"

/*!
 * \brief The BurialMedium class is a simple container class that defines the
 * conductivity, density and heat capacity of a burial medium. Mostly used by
 * HeatTransfer.
 *
 * BurialMedium is a subclass of Material, but does not extend it in any way.
 */
class BurialMedium : public Material
{
public:
    /*!
     * \brief Basic constructor that requires all material properties.
     * \param conductivity Thermal conductivity [W/(m K)]
     * \param density Density [kg/m3]
     * \param heatCapacity Heat capacity at constant pressure (\f$c_p\f$) [J/(kg K)]
     */
    constexpr BurialMedium(
            const double conductivity,
            const double density,
            const double heatCapacity):
        Material(conductivity, density, heatCapacity)
    {}

    /*!
     * \brief Construct from Material (explicit).
     * \param material Material to copy.
     */
    constexpr explicit BurialMedium(
            const Material& material):
        Material(material)
    {}

    // setters
//    double& conductivity() { return m_conductivity; }
//    double& density() { return m_density; }
//    double& heatCapacity() { return m_heatCapacity; }

    // initialized below
    static const BurialMedium soil; //!< Predefined soil
};

constexpr inline BurialMedium BurialMedium::soil = BurialMedium(Material::soil);
