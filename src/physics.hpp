#pragma once

#include <string>
#include <memory>
#include <armadillo>

#include "heattransfer/heattransfer.hpp"
#include "config.hpp"

class Pipeline;
class EquationOfState;

/*!
 * \brief The Physics class combines EquationOfState and HeatTransfer to
 * calculate the new state of a pipeline given a new pressure, temperature
 * and flow.
 *
 * This class operates on a Pipeline instance, and returns a new Pipeline
 * instance with updated properties (heat flow, Reynolds number etc.)
 */
class Physics
{
public:
    //! Declared to avoid the inline compiler-generated default destructor.
    virtual ~Physics();

    /*!
     * \brief Construct from Config and Pipeline.
     * \param state Pipeline instance
     * \param config Config instance
     */
    Physics(
            const Pipeline& state,
            const Config& config);

    /*!
     * \brief Construct from Pipeline and strings to select the equation of
     * state and heat transfer type to use.
     * \param state Pipeline description
     * \param eos Equation of state type ("BWRS", "GERG04" or "IdealGas")
     * \param heat Heat transfer type ("SteadyState", "Unsteady", "FixedUValue" or "FixedQValue")
     */
    explicit Physics(
            const Pipeline& state,
            const std::string& eos = "BWRS",
            const std::string& heat = "SteadyState");

    /*!
     * \brief Updates all derived properties, but does not evaluate the heat
     * transfer.
     * \param state Pipeline state to update
     */
    void updateDerivedProperties(
            Pipeline& state) const;

    /*!
     * \brief Initialize heat transfer of a Pipeline.
     * \param state Pipeline state to initializing heat transfer of
     */
    void initializeHeatTransferState(Pipeline& state) const;

    /*!
     * \brief Thermalize the heat transfer of a Pipeline.
     * \param state Pipeline state to thermalize
     */
    void thermalizeHeatTransfer(Pipeline& state) const;

    //! Get the size (number of grid points)
    arma::uword size() const { return m_heat->size(); }

    //! Get (const ref) EquationOfState
    const EquationOfState& equationOfState() const { return *m_eos; }

    //! Get (const ref) HeatTransfer
    const HeatTransfer& heatTransfer() const { return *m_heat; }

private:
    std::unique_ptr<EquationOfState> m_eos; //!< Equation of state
    std::unique_ptr<HeatTransfer> m_heat; //!< Heat transfer
};
