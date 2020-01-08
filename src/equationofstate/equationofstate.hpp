#pragma once

#include <vector>
#include <memory>
#include <armadillo>

class Pipeline;
class EquationOfStateBase;

/*!
 * \brief The EquationOfState class is a wrapper around EquationOfStateBase that
 * has one equation of state instance per grid point, and wraps the
 * EquationOfStateBase::evaluate() function.
 */
class EquationOfState
{
public:
    //! Declared to avoid the inline compiler-generated default destructor.
    ~EquationOfState();

    /*!
     * \brief Construct from pipeline and string to determine the type of
     * equation of state.
     * \param pipeline Pipeline instance
     * \param type Type of equation of state ("BWRS", "GERG04" or "IdealGas")
     */
    EquationOfState(const Pipeline& pipeline, const std::string& type);

    /*!
     * \brief This is a wrapper around EquationOfStateBase::evaluate() that
     * calls that function for each grid point.
     * \param state Pipeline instance
     * \return Matrix containing the output from each EquationOfStateBase
     * instance in each row. This way it's easy to update the relevant Pipeline
     * members, which are stored as column vectors.
     */
    arma::mat evaluate(const Pipeline& state);

    //! std::vector-like at(i) getter
    const EquationOfStateBase& at(std::size_t pos) const { return *m_eos->at(pos); }

    //! std::vector-like size() operator
    auto size() const { return m_eos->size(); }

private:
    //! Vector of EquationOfStateBase instances, one for each grid point.
    std::unique_ptr<std::vector<std::unique_ptr<EquationOfStateBase>>> m_eos;
};
