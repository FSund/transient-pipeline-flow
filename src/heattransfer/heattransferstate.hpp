#pragma once

#include <optional>
#include <armadillo>

/*!
 * \brief Container for the state of a HeatTransferBase instance. All
 * HeatTransferBase subclasses operate on instances of State, and return State.
 *
 * This uses std::optional for the temperature member, since this is only
 * used by UnsteadyHeatTransferBase.
 */
struct HeatTransferState
{
    /*!
     * \brief Construct from only heat flux (skip optional temperature vector).
     * \param heatFlux Heat flux between gas and surroundings [W/m2]
     */
    explicit HeatTransferState(const double heatFlux = 0);

    /*!
     * \brief Construct from heat flux and temperature vector, both with default values.
     * \param heatFlux Heat flux between gas and surroundings [W/m2]
     * \param temperature Pipe wall and burial medium temperatures when using 1d radial model [K]
     */
    explicit HeatTransferState(const double heatFlux, const arma::vec& temperature);

    /*!
     * \brief Temperature getter. This will throw an error if
     * State::m_temperature is not set, since this is optional. Check
     * State::hasTemperature() first if you want to be sure.
     *
     * \return Temperature [K]
     */
    const arma::vec& temperature() const { return m_temperature.value(); }

    /*!
     * \brief State::m_temperature is optional, so this getter returns true
     * if State::m_temperature has been set properly.
     *
     * \return True if State::m_temperature is set (optional).
     */
    bool hasTemperature() { return m_temperature.has_value(); }

    //! Heat flux getter [W/m2]
    double heatFlux() const { return m_heatFlux; }

    //! Set the temperature. Needs a setter since temperature is std::optional.
    void setTemperature(const arma::vec& temperature);

private:
    double m_heatFlux; //!< Heat flux [W/m2]
    std::optional<arma::vec> m_temperature; //!< Pipe wall temperature (optional) [K]
};
