#pragma once

/*!
 * \brief The Material class is a simple container class that defines the
 * conductivity, density and heat capacity of a material.
 */
class Material
{
public:
    /*!
     * \brief Basic constructor that requires all material properties.
     * \param conductivity Thermal conductivity [W/(m K)]
     * \param density Density [kg/m3]
     * \param heatCapacity Heat capacity at constant pressure (\f$c_p\f$) [J/(kg K)]
     */
    constexpr Material(
            const double conductivity,
            const double density,
            const double heatCapacity):
        m_conductivity(conductivity),
        m_density(density),
        m_heatCapacity(heatCapacity)
    {}

    // read-only
    double conductivity() const { return m_conductivity; } //!< Conductivity getter
    double density() const { return m_density; } //!< Density getter
    double heatCapacity() const { return m_heatCapacity; } //!< Heat capacity getter

    // initialized below
    static const Material concrete; //!< Predefined concrete
    static const Material steel; //!< Predefined steel
    static const Material coating; //!< Predefined pipe coating
    static const Material soil; //!< Predefined soil
    static const Material seawater; //!< Predefined seawater
    static const Material air; //!< Predefined air

protected:
    double m_conductivity; //!< Thermal conductivity [W/(m K)]
    double m_density; //!< Density [kg/m3]
    double m_heatCapacity; //!< Heat capacity at constant pressure (\f$c_p\f$) [J/(kg K)]
};

constexpr inline Material Material::concrete = Material( 2.9,    3400,   650);
constexpr inline Material Material::steel    = Material(50,      7800,   590);
constexpr inline Material Material::coating  = Material( 0.74,   1300,  1900);
constexpr inline Material Material::soil     = Material( 2,      2000,  1000); // water saturated sand

constexpr inline Material Material::seawater = Material( 0.571,  1020,  4187);
constexpr inline Material Material::air      = Material( 0.0257, 1.225, 1012);
