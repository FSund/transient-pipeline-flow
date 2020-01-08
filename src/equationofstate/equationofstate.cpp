#include "equationofstate.hpp"

#include "pipeline.hpp"
#include "equationofstate/bwrs.hpp"
#include "equationofstate/gerg04.hpp"
#include "equationofstate/idealgas.hpp"
#include "equationofstate/dummygas.hpp"

using std::cout;
using std::endl;
using std::vector;
using std::make_unique;
using std::unique_ptr;
using arma::uword;
using arma::mat;
using arma::vec;

EquationOfState::~EquationOfState()
{}

EquationOfState::EquationOfState(
        const Pipeline& state,
        const std::string& eos):
    m_eos(make_unique<vector<unique_ptr<EquationOfStateBase>>>())
{
    for (uword i = 0; i < state.size(); i++)
    {
        if (eos == "BWRS")
        {
            m_eos->push_back(make_unique<BWRS>(state.composition().at(i)));
        }
        else if (eos == "GERG04")
        {
            m_eos->push_back(make_unique<GERG04>(state.composition().at(i)));
        }
        else if (eos == "IdealGas")
        {
            m_eos->push_back(make_unique<IdealGas>(state.composition().at(i)));
        }
        else if (eos == "DummyGas")
        {
            m_eos->push_back(make_unique<DummyGas>(state.composition().at(i)));
        }
        else
        {
            std::string what = "invalid EOS type \"" + eos + "\"";
            throw std::invalid_argument(what);
        }
    }
}

arma::mat EquationOfState::evaluate(const Pipeline& state)
{
    if (state.size() != m_eos->size())
    {
        throw std::invalid_argument("state.size() != m_eos->size()");
    }

    mat output(state.size(), 7); // memory is not initialized
    for (uword i = 0; i < m_eos->size(); i++)
    {
        // using the simpler evaluate (without composition argument) can be
        // faster, but since we don't know if composition of Pipeline has
        // changed since previous iteration we just use this version to be sure
        // we get the correct answer
        const vec out = m_eos->at(i)->evaluate(state.pressure()(i), state.temperature()(i), state.composition().at(i));
        output.row(i).cols(0, 5) = out.t();

        output(i, 6) = m_eos->at(i)->getMolarMassOfMixture();
    }

    return output;
}
