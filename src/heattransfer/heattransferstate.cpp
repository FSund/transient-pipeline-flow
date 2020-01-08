#include "heattransferstate.hpp"

HeatTransferState::HeatTransferState(const double heatFlux):
    m_heatFlux(heatFlux)
{}

HeatTransferState::HeatTransferState(const double heatFlux, const arma::vec& temperature):
    m_heatFlux(heatFlux),
    m_temperature(temperature)
{}

void HeatTransferState::setTemperature(const arma::vec& temperature)
{
    m_temperature = temperature;
}
