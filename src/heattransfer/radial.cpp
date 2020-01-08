#include "heattransfer/radial.hpp"

#include "heattransfer/utils.hpp"
#include "heattransfer/pipewall.hpp"
#include "utilities/utilities.hpp"
#include "constants.hpp"

using arma::vec;
using arma::uvec;
using arma::uword;
using arma::zeros;

RadialHeatTransfer::RadialHeatTransfer(
        const double diameter,
        const PipeWall& pipeWall,
        const double burialDepth,
        const BurialMedium& burialMedium,
        const AmbientFluid& ambientFluid):
    m_diameter(diameter),
    m_burialDepth(burialDepth),
    m_burialMedium(burialMedium),
    m_ambientFluid(ambientFluid)
{
    if (pipeWall.size() <= 0)
        throw std::invalid_argument("PipeWall contains no elements");

    // transfer pipe wall properties directly
    // TODO: consider discretizing thick wall layers?
    m_width = zeros<vec>(pipeWall.size());
    m_conductivity = zeros<vec>(pipeWall.size());
    m_density = zeros<vec>(pipeWall.size());
    m_heatCapacity = zeros<vec>(pipeWall.size());
    for (uword i = 0; i < pipeWall.size(); i++)
    {
        m_width(i) = pipeWall.layer(i).thickness();
        m_conductivity(i) = pipeWall.layer(i).conductivity();
        m_density(i) = pipeWall.layer(i).density();
        m_heatCapacity(i) = pipeWall.layer(i).heatCapacity();
    }
    m_isBurialLayer = zeros<uvec>(pipeWall.size());

    // set up burial layers
    if (burialDepth > 0) // distance from top of pipe to top of burial medium
    {
        // calculate equivalent soil layers
        vec burialMediumWidths = utils::calcEquivalentBurialLayerWidths(
                    m_diameter, arma::sum(m_width), m_burialDepth,
                    burialMedium.conductivity());

        uword n = burialMediumWidths.n_elem;

        m_isBurialLayer = arma::join_vert(m_isBurialLayer, zeros<uvec>(n) + true);
        m_width =         arma::join_vert(m_width,         burialMediumWidths);
        m_conductivity =  arma::join_vert(m_conductivity,  zeros<vec>(n) + burialMedium.conductivity());
        m_density =       arma::join_vert(m_density,       zeros<vec>(n) + burialMedium.density());
        m_heatCapacity =  arma::join_vert(m_heatCapacity,  zeros<vec>(n) + burialMedium.heatCapacity());
    }

    // calculate inner and outer radii
    m_ri = zeros<vec>(size());
    m_ro = zeros<vec>(size());
    m_ri.fill(m_diameter/2.0);
    for (uword j = 1; j < size(); j++)
    {
        m_ri(j) += sum(m_width(arma::span(0, j-1)));
    }
    for (uword j = 0; j < size(); j++)
    {
        m_ro(j) = m_ri(j) + m_width(j);
    }

    // cross-sectional area of each shell
    m_crossSection = constants::pi*(pow(m_ro, 2.0) - pow(m_ri, 2.0)); // Cross section area of each shell == Ai
}

double RadialHeatTransfer::calculateOuterFilmCoefficient() const
{
    const double outerRadius = m_ro.tail(1)(0);
    const double outerDiameter = 2.0*outerRadius;
    const double ho = utils::calcOuterWallFilmCoefficient(outerDiameter, m_ambientFluid);

    return ho;
}
