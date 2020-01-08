#pragma once

#include "heattransfer/heattransferbase.hpp"

/*!
 * \brief Implementation of HeatTransferBase for fixed heat flux.
 */
class FixedQValue: public HeatTransferBase
{
public:
    /*!
     * \brief Construct with heat flux.
     * \param q Wanted heat flux [W/m2]
     */
    explicit FixedQValue(const double q):
        m_qValue(q)
    {}

    //! Set the Q value
    void setQValue(const double q) const
    {
        m_qValue = q;
    }

    /*!
     * \brief Evaluate override.
     * \return Constant heat flux q [W/m2]
     */
    virtual HeatTransferState evaluate(
            const HeatTransferState& /*current*/,
            const double /*timeStep*/,
            const double /*ambientTemperature*/,
            const double /*gasPressure*/,
            const double /*gasTemperature*/,
            const double /*gasReynoldsNumber*/,
            const double /*gasHeatCapacityConstantPressure*/,
            const double /*gasViscosity*/) const override
    {
        return evaluateInternal();
    }

    /*!
     * \brief Internal method, exposed for unit testing.
     * \return Constant heat flux [W/m2]
     */
    HeatTransferState evaluateInternal() const
    {
        return HeatTransferState(m_qValue);
    }

private:
    mutable double m_qValue; //!< Heat flux q [W/m2]
};
