#pragma once

#include "heattransfer/heattransferbase.hpp"

/*!
 * \brief Implementation of HeatTransferBase for fixed U-value (total heat transfer
 * coefficient).
 */
class FixedUValue: public HeatTransferBase
{
public:
    /*!
     * \brief Construct with U-value.
     * \param U Total heat transfer coefficient [W/(m2 K)]
     */
    explicit FixedUValue(const double U):
        m_uValue(U)
    {}

    //! Set the U value
    void setUValue(const double U) const
    {
        m_uValue = U;
    }

    /*!
     * \brief Evaluate override.
     * \param ambientTemperature Ambient temperature [K]
     * \param gasTemperature Gas temperature [K]
     * \return Heat flux [W/m2]
     */
    virtual HeatTransferState evaluate(
            const HeatTransferState& /*current*/,
            const double /*timeStep*/,
            const double ambientTemperature,
            const double /*gasPressure*/,
            const double gasTemperature,
            const double /*gasReynoldsNumber*/,
            const double /*gasHeatCapacityConstantPressure*/,
            const double /*gasViscosity*/) const override
    {
        return evaluateInternal(gasTemperature, ambientTemperature);
    }

    /*!
     * \brief Internal method, exposed for unit testing.
     *
     * Calculates the heat flux q for the given ambient temperature and gas
     * temperature, from U-value m_uValue.
     *
     * \param gasTemperature Gas temperature [K]
     * \param ambientTemperature Ambient temperature [K]
     * \return Heat flux q [W/m2]
     */
    HeatTransferState evaluateInternal(const double gasTemperature, const double ambientTemperature) const
    {
        const double q = m_uValue*(gasTemperature - ambientTemperature);
        return HeatTransferState(q);
    }

private:
    mutable double m_uValue; //!< Total heat transfer coefficient [W/(m2 K)]
};
