#pragma once

#include <memory>
#include <vector>

class HeatTransferBase;
class HeatTransferState;
class Pipeline;

/*!
 * \brief The EquationOfState class is a wrapper around HeatTransferBase that
 * has one heat transfer instance per grid point, and wraps the
 * HeatTransferBase::evaluate() function.
 */
class HeatTransfer
{
public:
    //! Declared to avoid the inline compiler-generated default destructor.
    virtual ~HeatTransfer();

    /*!
     * \brief Construct from pipeline and string to determine the type of
     * heat transfer.
     * \param pipeline Pipeline instance
     * \param type Type of heat transfer ("BWRS", "GERG04" or "IdealGas")
     */
    HeatTransfer(const Pipeline& pipeline, const std::string& type);

    /*!
     * \brief This is a wrapper around HeatTransferBase::evaluate() that
     * calls that function for each grid point. This modifes the pipeline
     * argument.
     * \param state Current heat transfer state
     * \param timeStep Time step [s]
     * \param pipeline Pipeline instance
     */
    void evaluate(const std::vector<HeatTransferState>& state, const double timeStep, Pipeline& pipeline) const;

    //! std::vector-like at(i) getter
    const HeatTransferBase& at(std::size_t pos) const { return *m_heat->at(pos); }

    //! std::vector-like size() operator
    auto size() const { return m_heat->size(); }

private:
    //! Vector of HeatTransferBase instances, one for each grid point.
    std::unique_ptr<std::vector<std::unique_ptr<HeatTransferBase>>> m_heat;
};
