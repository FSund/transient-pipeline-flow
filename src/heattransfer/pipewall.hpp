#pragma once

#include <vector>
#include <armadillo>

#include "material.hpp"

/*!
 * \brief The PipeWall class is a container class that defines the thickness
 * and Material properties of each layer a pipe consists of.
 */
class PipeWall
{
public:
    class Layer;

    /*!
     * \brief Construct PipeWall with a given number of layers.
     * \param nWallLayers Number of layers.
     */
    explicit PipeWall(const arma::uword nWallLayers):
        m_layers(std::vector<Layer>(nWallLayers))
    {}

    /*!
     * \brief Construct PipeWall from a vector of PipeWall::Layer.
     * \param wallLayers The layers of the wall.
     */
    explicit PipeWall(const std::vector<Layer>& wallLayers):
        m_layers(wallLayers)
    {}

    //! Get layer with index i. Returns const reference to member.
    const Layer& layer(const arma::uword i) const
    {
        return m_layers.at(i);
    }

    //! Get layer with index i. Returns (non-const) reference to member.
    Layer& layer(const arma::uword i)
    {
        return m_layers.at(i);
    }

    //! Get all pipe wall layers. Returns const reference to member.
    const std::vector<Layer>& layers() const { return m_layers; }

    //! Get number of layers in the pipe wall.
    arma::uword size() const { return m_layers.size(); }

    // initialized in .cpp
    static const PipeWall defaultPipeWall; //!< Predefined pipe wall.

private:
    //! m_layers Vector of PipeWall::Layer that the wall consists of.
    std::vector<Layer> m_layers;
};

/*!
 * \brief The PipeWall::Layer class is a simple container class that defines the
 * thickness and all other material properties of a single layer of pipe wall.
 */
class PipeWall::Layer : public Material
{
public:
    /*!
     * \brief Default constructor, sets all properties to -1.
     */
    constexpr Layer():
        Material(-1, -1, -1),
        m_thickness(-1)
    {}

    /*!
     * \brief Basic constructor that requires all material properties.
     * \param thickness Layer thickness [m]
     * \param conductivity Thermal conductivity [W/(m K)]
     * \param density Density [kg/m3]
     * \param heatCapacity Heat capacity at constant pressure (\f$c_p\f$) [J/(kg K)]
     */
    constexpr Layer(
            const double thickness,
            const double conductivity,
            const double density,
            const double heatCapacity):
        Material(conductivity, density, heatCapacity),
        m_thickness(thickness)
    {}

    /*!
     * \brief Construct from thickness and Material instance.
     * \param thickness Layer thickness [m]
     * \param material Pipe wall material.
     */
    constexpr Layer(
            const double thickness,
            const Material& material):
        Material(material),
        m_thickness(thickness)
    {}

    double thickness() const { return m_thickness; } //!< Thickness getter
    double conductivity() const { return m_conductivity; } //!< Conductivity getter
    double density() const { return m_density; } //!< Density getter
    double heatCapacity() const { return m_heatCapacity; } //!< Heat capacity getter

    // setters
    double& thickness() { return m_thickness; } //!< Thickness setter
    double& conductivity() { return m_conductivity; } //!< Conductivity setter
    double& density() { return m_density; } //!< Density setter
    double& heatCapacity() { return m_heatCapacity; } //!< Heat capacity setter

private:
    double m_thickness; //!< Layer wall thickness [m]
};
