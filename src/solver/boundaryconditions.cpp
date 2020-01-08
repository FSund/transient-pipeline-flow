#include "solver/boundaryconditions.hpp"

#include "pipeline.hpp"

using arma::uword;
typedef BoundaryConditions::SingleCondition SingleCondition;

BoundaryConditions::BoundaryConditions():
    BoundaryConditions(
        arma::zeros<arma::mat>(3, 2),
        Composition(),
        Composition()
    )
{}

BoundaryConditions::BoundaryConditions(
        const arma::mat& bcMat,
        const Composition& inletComposition,
        const Composition& outletComposition):
    m_inletFlow(bcMat(0, 0), true),
    m_inletPressure(bcMat(1, 0), false),
    m_inletTemperature(bcMat(2, 0), true),
    m_outletFlow(bcMat(0, 1), false),
    m_outletPressure(bcMat(1, 1), true),
    m_outletTemperature(bcMat(2, 1), false),
    m_inletComposition(inletComposition),
    m_outletComposition(outletComposition)
{
    if (bcMat.n_cols != 2)
        throw std::runtime_error("invalid number of columns (should be exactly 2)");
    if (bcMat.n_rows != 3)
        throw std::runtime_error("invalid number of rows (should be exactly 3)");
}

BoundaryConditions::BoundaryConditions(
        const double inletFlow,
        const double outletFlow,
        const double inletPressure,
        const double outletPressure,
        const double inletTemperature,
        const double outletTemperature,
        const Composition& inletComposition,
        const Composition& outletComposition):
    m_inletFlow(inletFlow, true),
    m_inletPressure(inletPressure, false),
    m_inletTemperature(inletTemperature, true),
    m_outletFlow(outletFlow, false),
    m_outletPressure(outletPressure, true),
    m_outletTemperature(outletTemperature, false),
    m_inletComposition(inletComposition),
    m_outletComposition(outletComposition)
{}

BoundaryConditions::BoundaryConditions(
        const SingleCondition inletFlow,
        const SingleCondition outletFlow,
        const SingleCondition inletPressure,
        const SingleCondition outletPressure,
        const SingleCondition inletTemperature,
        const SingleCondition outletTemperature,
        const Composition& inletComposition,
        const Composition& outletComposition):
    m_inletFlow(inletFlow),
    m_inletPressure(inletPressure),
    m_inletTemperature(inletTemperature),
    m_outletFlow(outletFlow),
    m_outletPressure(outletPressure),
    m_outletTemperature(outletTemperature),
    m_inletComposition(inletComposition),
    m_outletComposition(outletComposition)
{}

BoundaryConditions::BoundaryConditions(
        const Pipeline& state,
        const std::vector<std::string>& boundarySettings):
    BoundaryConditions(
        state.flow()(0),
        state.flow().tail(1)(0),
        state.pressure()(0),
        state.pressure().tail(1)(0),
        state.temperature()(0),
        state.temperature().tail(1)(0),
        state.composition().front(),
        state.composition().back()
        )
{
    setBoundarySettings(boundarySettings);
}

void BoundaryConditions::setBoundarySettings(const std::vector<std::string>& strings)
{
    if (strings.size() != 3)
        throw std::runtime_error("invalid number of strings (should be exactly 3)");

    setBoundarySettings(0, strings.at(0));
    setBoundarySettings(1, strings.at(1));
    setBoundarySettings(2, strings.at(2));
}

arma::uword BoundaryConditions::nActiveBoundaryConditions() const
{
    return uword(m_inletFlow.isActive()) + uword(m_outletFlow.isActive())
            + uword(m_inletPressure.isActive()) + uword(m_outletPressure.isActive())
            + uword(m_inletTemperature.isActive()) + uword(m_outletTemperature.isActive());
}

std::ostream& operator <<(std::ostream& out, const BoundaryConditions& c)
{
    out << "Flow:        " << c.inletFlow().value()        << ", " << c.outletFlow().value() << std::endl;
    out << "Pressure:    " << c.inletPressure().value()    << ", " << c.outletPressure().value() << std::endl;
    out << "Temperature: " << c.inletTemperature().value() << ", " << c.outletTemperature().value() << std::endl;
    out << "Composition: " << std::endl;
    out << arma::join_horiz(arma::vec(c.m_inletComposition), arma::vec(c.m_outletComposition));

    return out;
}
