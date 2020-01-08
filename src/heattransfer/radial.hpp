#pragma once

#include <armadillo>

#include "heattransfer/heattransferbase.hpp"
#include "heattransfer/burialmedium.hpp"
#include "heattransfer/ambientfluid.hpp"

class Pipeline;
class PipeWall;

/*!
 * \brief Base class for heat transfer calculation with 1d radial models.
 *
 * This contains the description of the burial medium, burial depth, ambient
 * fluid etc., as well as the 1d radial discretization of the surroundings of
 * the pipeline (in m_width, m_conductivity, etc.).
 */
class RadialHeatTransfer : public HeatTransferBase
{
public:
    /*!
     * \brief Constructor that sets up the discretization of the pipe surroundings.
     * \param diameter Pipe inner diameter [m]
     * \param pipeWall PipeWall instance
     * \param burialDepth Distance from top of pipe to top of burial medium [m]
     * \param burialMedium BurialMedium instance
     * \param ambientFluid AmbientFluid instance
     */
    RadialHeatTransfer(
            const double diameter,
            const PipeWall& pipeWall,
            const double burialDepth,
            const BurialMedium& burialMedium,
            const AmbientFluid& ambientFluid);

    /*!
     * \brief Calculate the outer film coefficient. Basically a wrapper around
     * utils::calcOuterWallFilmCoefficient(const double, const AmbientFluid&)
     * using the proper arguments.
     *
     * \return Outer film coefficient [W/(m2 K)]
     */
    double calculateOuterFilmCoefficient() const;

    //! The number of discretization elements.
    virtual arma::uword size() const { return m_width.n_elem; }

    /*!
     * \brief Make instance of HeatTransferState from heat flux. Override.
     *
     * Makes HeatTransferState instance with the correct temperature property
     * for 1d radial heat transfer models.
     *
     * \param heatFlux Heat flux [W/m2]
     * \return HeatTransferState instance with temperature.
     */
    virtual HeatTransferState makeState(const double heatFlux) const override
    {
        return HeatTransferState(heatFlux, arma::zeros<arma::vec>(size()));
    }

    /*!
     * \brief Make instance of HeatTransferState from heat flux. Override.
     *
     * Makes HeatTransferState instance with the correct temperature property
     * for 1d radial heat transfer models. Performs a linear interpolation
     * between gas temperature and ambient temperature. This could probably be
     * improved upon, since there is likely an analytical solution.
     *
     * \param heatFlux Heat flux [W/m2]
     * \param ambientTemperature Ambient temperature [K]
     * \param gasTemperature Gas temperature [K]
     * \return HeatTransferState instance with temperature.
     */
    virtual HeatTransferState makeState(
            const double heatFlux,
            const double gasTemperature,
            const double ambientTemperature) const override
    {
        arma::vec shellTemperature = arma::linspace(
                    gasTemperature, ambientTemperature, size());

        return HeatTransferState(heatFlux, shellTemperature);
    }

protected:
    double m_diameter; //!< Pipe inner diameter [m]

    double m_burialDepth; //!< Distance from top of pipe to top of burial medium [m]

    //! Description of the medium the pipeline is buried in (or on top
    //! of/suspended above depending on burial depth).
    BurialMedium m_burialMedium;
    //! Description of the fluid surrounding the pipeline.
    AmbientFluid m_ambientFluid;

    arma::vec m_width;        //!< Width of each discretization shell [m]
    arma::vec m_conductivity; //!< Thermal conductivity of each discretization shell [W/(m K)]
    arma::vec m_density;      //!< Density of each discretization shell [kg/m3]
    arma::vec m_heatCapacity; //!< Heat capacity of each discretization shell (\f$c_p\f$) [J/(kg K)]

    //! Bool that determines whether the discretization layers are due to burial
    //! in the burial medium or not.
    arma::uvec m_isBurialLayer;

    arma::vec m_crossSection; //!< Area/cross-section of each shell [m2]
    arma::vec m_ri;           //!< Inner radius [m]
    arma::vec m_ro;           //!< Outer radius [m]
};
