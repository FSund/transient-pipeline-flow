#include "heattransfer/heattransferbase.hpp"

HeatTransferBase::~HeatTransferBase()
{}

HeatTransferState HeatTransferBase::makeState(const double heatFlux) const
{
    return HeatTransferState(heatFlux);
}

HeatTransferState HeatTransferBase::makeState(
        const double heatFlux,
        const double /*gasTemperature*/,
        const double /*ambientTemperature*/) const
{
    return makeState(heatFlux);
}
