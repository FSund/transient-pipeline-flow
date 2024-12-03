#include "timeseries.hpp"

#include "solver/boundaryconditions.hpp"

typedef BoundaryConditions::SingleCondition SingleCondition;

TimeSeries::TimeSeries(const arma::uword size, const arma::uword dt):
    m_timestamps(arma::linspace<arma::uvec>(0, dt*(size - 1), size)),
    m_inletFlow(arma::vec(size), false),
    m_inletPressure(arma::vec(size), false),
    m_inletTemperature(arma::vec(size), false),
    m_outletFlow(arma::vec(size), false),
    m_outletPressure(arma::vec(size), false),
    m_outletTemperature(arma::vec(size), false),
    m_inletComposition(std::vector<Composition>(size, Composition::defaultComposition)),
    m_outletComposition(std::vector<Composition>(size, Composition::defaultComposition))
{}

TimeSeries::TimeSeries(const arma::uvec& timestamps):
    TimeSeries(timestamps, std::vector<BoundaryConditions>(timestamps.n_elem))
{}

TimeSeries::TimeSeries(
        const Pipeline& pipeline,
        const arma::uword size,
        const arma::uword dt,
        const std::vector<std::string>& boundarySettings
        ):
    TimeSeries(
        dt,
        std::vector<BoundaryConditions>(size, BoundaryConditions(pipeline))
    )
{
    this->setBoundarySettings(boundarySettings);
}

TimeSeries::TimeSeries(
        const std::string& filename,
        const std::vector<std::string>& boundarySettings):
    TimeSeries(filename, 0, 0, boundarySettings)
{}


TimeSeries::TimeSeries(
        const std::string& filename,
        arma::uword lastRow,
        const std::vector<std::string>& boundarySettings):
    TimeSeries(filename, 0, lastRow, boundarySettings)
{}

void TimeSeries::loadFromMatrix(
        const arma::mat& bc)
{
    bool hasComposition;
    if (bc.n_cols == 1 + 6)
    {
        hasComposition = false;
    }
    else if (bc.n_cols == 1 + 6 + 20)
    {
        hasComposition = true;
    }
    else
    {
        throw std::runtime_error("invalid number of columns");
    }

    m_timestamps = arma::conv_to<arma::uvec>::from(arma::round(bc.col(0)));
    if (hasComposition)
    {
        for (std::size_t i = 0; i < bc.n_rows; i++)
        {
            m_inletComposition.push_back(Composition(bc.row(i).cols(4, 4 + 9).t()).normalize());
            m_outletComposition.push_back(Composition(bc.row(i).cols(17, 17 + 9).t()).normalize());
        }
        m_inletFlow.set(bc.col(1), true);
        m_inletPressure.set(bc.col(2), false);
        m_inletTemperature.set(bc.col(3), true);
        m_outletFlow.set(bc.col(14), false);
        m_outletPressure.set(bc.col(15), true);
        m_outletTemperature.set(bc.col(16), false);
    }
    else
    {
        m_inletComposition = std::vector<Composition>(bc.n_rows, Composition::defaultComposition);
        m_outletComposition = std::vector<Composition>(bc.n_rows, Composition::defaultComposition);
        m_inletFlow.set(bc.col(1), true);
        m_inletPressure.set(bc.col(2), false);
        m_inletTemperature.set(bc.col(3), true);
        m_outletFlow.set(bc.col(4), false);
        m_outletPressure.set(bc.col(5), true);
        m_outletTemperature.set(bc.col(6), false);
    }
}

TimeSeries::TimeSeries(
        const std::string& filename,
        arma::uword firstRow,
        arma::uword lastRow,
        const std::vector<std::string>& boundarySettings)
{   
    arma::mat bc;
    if (!bc.load(filename, arma::csv_ascii))
    {
        throw std::runtime_error("could not read file \"" + filename + "\"");
    }

    if (lastRow == 0)
    {
        lastRow = bc.n_rows - 1;
    }

    if (lastRow <= firstRow)
    {
        throw std::invalid_argument("lastRow <= firstRow");
    }

    if (lastRow >= bc.n_rows)
    {
        throw std::invalid_argument("lastRow more than the number of rows");
    }

    // discard unwanted rows
    bc = bc.rows(firstRow, lastRow);

    loadFromMatrix(bc);
    setBoundarySettings(boundarySettings);
}

TimeSeries::TimeSeries(
        const arma::mat& bc,
        const std::vector<std::string>& boundarySettings)
{
    loadFromMatrix(bc);
    setBoundarySettings(boundarySettings);
}

TimeSeries::TimeSeries(
        const arma::uvec& timestamps,
        const std::vector<BoundaryConditions>& boundaryConditions):
    m_timestamps(timestamps),
    m_inletFlow(arma::vec(boundaryConditions.size()), true),
    m_inletPressure(arma::vec(boundaryConditions.size()), false),
    m_inletTemperature(arma::vec(boundaryConditions.size()), true),
    m_outletFlow(arma::vec(boundaryConditions.size()), false),
    m_outletPressure(arma::vec(boundaryConditions.size()), true),
    m_outletTemperature(arma::vec(boundaryConditions.size()), false)
{
    if (timestamps.size() != boundaryConditions.size())
    {
        throw std::invalid_argument("lengths do not match");
    }

    m_inletComposition = std::vector<Composition>();
    m_outletComposition = std::vector<Composition>();
    for (std::size_t i = 0; i < boundaryConditions.size(); i++)
    {
        m_inletFlow.vec()(i)         = boundaryConditions.at(i).inletFlow();
        m_inletPressure.vec()(i)     = boundaryConditions.at(i).inletPressure();
        m_inletTemperature.vec()(i)  = boundaryConditions.at(i).inletTemperature();

        m_outletFlow.vec()(i)        = boundaryConditions.at(i).outletFlow();
        m_outletPressure.vec()(i)    = boundaryConditions.at(i).outletPressure();
        m_outletTemperature.vec()(i) = boundaryConditions.at(i).outletTemperature();

        m_inletComposition.push_back(boundaryConditions.at(i).inletComposition());
        m_outletComposition.push_back(boundaryConditions.at(i).outletComposition());
    }
}

TimeSeries::TimeSeries(
        const arma::uword dt,
        const std::vector<BoundaryConditions>& boundaryConditions):
    TimeSeries(
        arma::linspace<arma::uvec>(0, dt*(boundaryConditions.size() - 1), boundaryConditions.size()),
        boundaryConditions
        )
{}

void TimeSeries::setBoundarySettings(const std::vector<std::string> strings)
{
    if (strings.size() != 3)
        throw std::runtime_error("invalid number of strings (should be exactly 3)");

    for (arma::uword i = 0; i < strings.size(); i++)
    {
        // lambda
        Series* inlet = [this](const arma::uword i){
            if (i == 0)
                return &m_inletFlow;
            else if (i == 1)
                return &m_inletPressure;
            else if (i == 2)
                return &m_inletTemperature;
            else
                throw std::runtime_error("i not in (0, 1, 2)");
        }(i); // call lambda immediately

        // lambda
        Series* outlet = [this](const arma::uword i){
            if (i == 0)
                return &m_outletFlow;
            else if (i == 1)
                return &m_outletPressure;
            else if (i == 2)
                return &m_outletTemperature;
            else
                throw std::runtime_error("i not in (0, 1, 2)");
        }(i); // call lambda immediately

        if (strings.at(i) == "none")
        {
            inlet->setActive(false);
            outlet->setActive(false);
        }
        else if (strings.at(i) == "inlet")
        {
            inlet->setActive(true);
            outlet->setActive(false);
        }
        else if (strings.at(i) == "outlet")
        {
            inlet->setActive(false);
            outlet->setActive(true);
        }
        else if (strings.at(i) == "both")
        {
            inlet->setActive(true);
            outlet->setActive(true);
        }
        else
        {
            throw std::runtime_error("invalid setting \"" + strings.at(i) + "\"");
        }
    }
}

void TimeSeries::save(const std::string& filename)
{
    // make matrix
    arma::mat data(size(), 1 + 6 + 10 + 10);
    data.col(0) = arma::conv_to<arma::vec>::from(m_timestamps);

    // inlet
    data.col(1) = m_inletFlow.vec();
    data.col(2) = m_inletPressure.vec();
    data.col(3) = m_inletTemperature.vec();

    // outlet
    data.col(14) = m_outletFlow.vec();
    data.col(15) = m_outletPressure.vec();
    data.col(16) = m_outletTemperature.vec();

    for (size_t i = 0; i < size(); i++)
    {
        // could probably skip this loop and assign rowvec via span instead
        for (size_t j = 0; j < 10; j++) // loop over components
        {
            data(i, 4 + j) = m_inletComposition.at(i)(j);
            data(i, 17 + j) = m_outletComposition.at(i)(j);
        }
    }

    data.save(filename, arma::csv_ascii);
}

BoundaryConditionsStamped TimeSeries::at(std::size_t i) const
{
    BoundaryConditions bc(
        SingleCondition(m_inletFlow(i), m_inletFlow.isActive()),
        SingleCondition(m_outletFlow(i), m_outletFlow.isActive()),
        SingleCondition(m_inletPressure(i), m_inletPressure.isActive()),
        SingleCondition(m_outletPressure(i), m_outletPressure.isActive()),
        SingleCondition(m_inletTemperature(i), m_inletTemperature.isActive()),
        SingleCondition(m_outletTemperature(i), m_outletTemperature.isActive()),
        m_inletComposition.at(i),
        m_outletComposition.at(i)
    );

    return BoundaryConditionsStamped(m_timestamps(i), bc);
}

TimeSeries::operator std::vector<BoundaryConditionsStamped>() const
{
    // user-defined conversion
    std::vector<BoundaryConditionsStamped> timeSteps;
    for (std::size_t i = 0; i < this->size(); i++)
    {
        timeSteps.push_back(this->at(i));
    }

    return timeSteps;
}
