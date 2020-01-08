#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <armadillo>

#include "composition.hpp"

class Pipeline;
class BoundaryConditions;
class TimeStep;

/*!
 * \brief The TimeSeries class is used to store the boundary conditions for
 * several time steps.
 *
 * The boundary conditions are stored in a arma::mat for convenience, so we can
 * set the individual components like inlet flow via
 *
 *     TimeSeries ts(100);
 *     ts.inletFlow() = arma::linspace<arma::vec>(0, 250, 100);
 *     ts.outletPressure().fill(1e6);
 *     ts.inletTemperature() = arma::zeros<vec>(100) + 273.15;
 *
 * There is a user-defined conversion to vector<BoundaryConditions> which is
 * used when simulating each time step.
 */
class TimeSeries
{
public:
    /*!
     * \brief Construct time series with a given number of steps, and a fixed
     * time step.
     * \param nSteps Number of steps
     * \param dt Time step
     */
    TimeSeries(const arma::uword nSteps, const arma::uword dt = 60);

    /*!
     * \brief Construct from timestamps.
     * \param timestamps Timestamps
     */
    TimeSeries(const arma::uvec& timestamps);

    /*!
     * \brief Construct time series from Pipeline, copying flow, pressure,
     * temperature and composition from Pipeline.
     * \param pipeline Pipeline to copy boundary conditions from
     * \param nSteps Number of steps
     * \param dt Time step
     * \param boundarySettings Boundary settings to use
     */
    TimeSeries(
            const Pipeline& pipeline,
            const arma::uword nSteps,
            const arma::uword dt = 60,
            const std::vector<std::string>& boundarySettings = {"inlet", "outlet", "inlet"});

    /*!
     * \brief Wrapper for TimeSeries(const std::string&, arma::uword, arma::uword, const std::vector<std::string>&)
     * which loads the whole file.
     *
     * \see TimeSeries(const std::string&, arma::uword, arma::uword, const std::vector<std::string>>&)
     *
     * \param filename Path to CSV-file
     * \param boundarySettings Boundary settings to use
     */
    TimeSeries(
            const std::string& filename,
            const std::vector<std::string>& boundarySettings = {"inlet", "outlet", "inlet"});

    /*!
     * \brief Wrapper for TimeSeries(const std::string&, arma::uword, arma::uword, const std::vector<std::string>&),
     * with the option to skip the end of the csv-file.
     *
     * \see TimeSeries(const std::string&, arma::uword, arma::uword, const std::vector<std::string>>&)
     *
     * \param filename Path to CSV-file
     * \param lastRow Last row of CSV-file to load (zero-indexed)
     * \param boundarySettings Boundary settings to use
     */
    TimeSeries(
            const std::string& filename,
            arma::uword lastRow,
            const std::vector<std::string>& boundarySettings = {"inlet", "outlet", "inlet"});

    /*!
     * \brief Load boundary conditions from a CSV.
     *
     * The CSV has to be in one of two formats; either 7 columns
     *  - timestamp [s]
     *  - inlet flow [kg/s]
     *  - inlet pressure [Pa]
     *  - inlet temperature [K]
     *  - outlet flow [kg/s]
     *  - outlet pressure[Pa]
     *  - outlet temperature [K]
     *
     * or 27 columns
     *  - timestamp [s]
     *  - inlet flow [kg/s]
     *  - inlet pressure [Pa]
     *  - inlet temperature [K]
     *  - inlet composition (10 cols) [fraction]
     *  - outlet flow [kg/s]
     *  - outlet pressure [Pa]
     *  - outlet temperature [K]
     *  - outlet composition (10 cols) [fraction]
     *
     * where inlet and outlet composition have the columns
     *  - C1, C2, C3, iC4, nC4, iC5, nC5, C6, N2, CO2
     *
     * \param filename Path to CSV-file
     * \param firstRow First row of CSV-file to load (zero-indexed)
     * \param lastRow Last row of CSV-file to load (zero-indexed)
     * \param boundarySettings Boundary settings to use
     */
    TimeSeries(
            const std::string& filename,
            arma::uword firstRow,
            arma::uword lastRow,
            const std::vector<std::string>& boundarySettings = {"inlet", "outlet", "inlet"});

    /*!
     * \brief Construct from matrix.
     *
     * The matrix has to be in one of two formats; either 7 columns
     *  - timestamp [s]
     *  - inlet flow [kg/s]
     *  - inlet pressure [Pa]
     *  - inlet temperature [K]
     *  - outlet flow [kg/s]
     *  - outlet pressure[Pa]
     *  - outlet temperature [K]
     *
     * or 27 columns
     *  - timestamp [s]
     *  - inlet flow [kg/s]
     *  - inlet pressure [Pa]
     *  - inlet temperature [K]
     *  - inlet composition (10 cols) [fraction]
     *  - outlet flow [kg/s]
     *  - outlet pressure [Pa]
     *  - outlet temperature [K]
     *  - outlet composition (10 cols) [fraction]
     *
     * where inlet and outlet composition have the columns
     *  - C1, C2, C3, iC4, nC4, iC5, nC5, C6, N2, CO2
     *
     * \param data Matrix with boundary conditions
     * \param boundarySettings Boundary settings to use
     */
    TimeSeries(
            const arma::mat& data,
            const std::vector<std::string>& boundarySettings = {"inlet", "outlet", "inlet"});

    /*!
     * \brief Construct from timestamps and boundary conditions.
     * \param timestamps Vector of timestamps [s]
     * \param boundaryConditions Vector of boundary conditions
     */
    TimeSeries(
            const arma::uvec& timestamps,
            const std::vector<BoundaryConditions>& boundaryConditions);

    /*!
     * \brief Construct with constant time step and vector of boundary conditions.
     * \param dt Time step [s]
     * \param boundaryConditions BoundaryConditions instance
     */
    TimeSeries(
            const arma::uword dt,
            const std::vector<BoundaryConditions>& boundaryConditions);

    /*!
     * \brief Set boundary settings via vector of string, or brace-init-list
     * like `{"inlet", "outlet", "inlet"}`
     * \param settings Vector of strings or brace-init-list, like `{"inlet", "outlet", "inlet"}`.
     *
     * Valid options are `"inlet"`, `"outlet"`, `"both"` and `"none"`.
     */
    void setBoundarySettings(const std::vector<std::string> settings);

    void save(const std::string& filename);

    // getters
    //! Get (const ref) timestamps.
    const auto& timestamps() const { return m_timestamps; }
    //! Get (const ref) inlet composition.
    const std::vector<Composition>& inletComposition() const { return m_inletComposition; }
    //! Get (const ref) outlet composition
    const std::vector<Composition>& outletComposition() const { return m_outletComposition; }

    const auto& inletFlow()         const { return m_inletFlow; } //!< Get (const ref) inlet flow.
    const auto& inletPressure()     const { return m_inletPressure; } //!< Get (const ref) inlet pressure.
    const auto& inletTemperature()  const { return m_inletTemperature; } //!< Get (const ref) inlet temperature.
    const auto& outletFlow()        const { return m_outletFlow; } //!< Get (const ref) outlet flow.
    const auto& outletPressure()    const { return m_outletPressure; } //!< Get (const ref) outlet pressure.
    const auto& outletTemperature() const { return m_outletTemperature; } //!< Get (const ref) outlet temperature.

    // getters/setters
    //! Get (ref) timestamps.
    auto& timestamps() { return m_timestamps; }
    //! Get (ref) inlet composition.
    std::vector<Composition>& inletComposition() { return m_inletComposition; }
    //! Get (ref) outlet composition
    std::vector<Composition>& outletComposition() { return m_outletComposition; }

    auto& inletFlow()         { return m_inletFlow; } //!< Get (ref) inlet flow.
    auto& inletPressure()     { return m_inletPressure; } //!< Get (ref) inlet pressure.
    auto& inletTemperature()  { return m_inletTemperature; } //!< Get (ref) inlet temperature.
    auto& outletFlow()        { return m_outletFlow; } //!< Get (ref) outlet flow.
    auto& outletPressure()    { return m_outletPressure; } //!< Get (ref) outlet pressure.
    auto& outletTemperature() { return m_outletTemperature; } //!< Get (ref) outlet temperature.

    //! User-defined conversion to vector of TimeStep. Convenient when we want
    //! to loop over all time steps, e.g. in Simulator::simulate().
    operator std::vector<TimeStep>() const;

    //! std::vector-like at(i) getter
    TimeStep at(std::size_t pos) const;

    //! Get size (number of grid points).
    arma::uword size() const { return m_timestamps.n_rows; }

    //! Contains a single property, and information regarding whether it is
    //! an active boundary condition or not
    class Series
    {
    public:
        friend class TimeSeries;

        //!
        //! \brief Default construct. Not active by default.
        //! \param active If it is an active boundary condition
        //!
        Series(const bool active=false):
            m_isActive(active)
        {}

        //!
        //! \brief Construct from arma::vec. Active by default.
        //! Allows implicit conversion from arma::vec to Series.
        //! \param value Value
        //!
        Series(const arma::vec& value):
            m_isActive(true),
            m_vec(value)
        {}

        //!
        //! \brief Construct from arma::vec and bool active.
        //! \param value Value
        //! \param active If it is an active boundary condition
        //!
        Series(const arma::vec& value, const bool active):
            m_isActive(active),
            m_vec(value)
        {}

        //! Set new value. Makes it an active boundary condition.
        void set(const arma::vec& value)
        {
            m_isActive = true;
            m_vec = value;
        }

        //! Set new value and active status.
        void set(const arma::vec& value, const bool active)
        {
            m_isActive = active;
            m_vec = value;
        }

        //! Armadillo-like fill. Fills member vector with value, and sets it active.
        void fill(const double value)
        {
            m_isActive = true;
            m_vec.fill(value);
        }

        //! Set the active/inactive property
        void setActive(const bool active)
        {
            m_isActive = active;
        }

        // copy assignment operator
        //!
        //! \brief Copy assignment operator from arma::vec. Sets it active.
        //!
        //! Allows for example
        //!
        //!     Series s;
        //!     s = arma::zeros<vec>(10) + 100;
        //!
        //! \param value New value
        //! \return Series instance
        //!
        Series& operator= (const arma::vec& value)
        {
            m_vec = value;
            m_isActive = true;

            return *this;
        }

        //! Get element by index
        double operator()(const arma::uword i) const { return m_vec(i); }

        //! Get the size of the member arma::vec
        arma::uword size() const { return m_vec.n_elem; }

        //! If the Series is an active boundary condition
        bool isActive() const { return m_isActive; }

        //! If the Series is an active boundary condition
        bool isActive() { return m_isActive; }

        //! Get (const ref) member arma::vec m_vec
        const arma::vec& vec() const { return m_vec; }

    private:
        //! Get (non-const ref) member arma::vec m_vec (private)
        arma::vec& vec() { return m_vec; }

        //! If the Series is an active boundary condition
        bool m_isActive = false;

        //! The value of the Series.
        arma::vec m_vec; // not initialized
    };

private:
    arma::uvec m_timestamps; //!< Timestamps [s]

    Series m_inletFlow; //!< Inlet flow [kg/s]
    Series m_inletPressure; //!< Inlet pressure [Pa]
    Series m_inletTemperature; //!< Inlet temperature [K]
    Series m_outletFlow; //!< Outlet flow [kg/s]
    Series m_outletPressure; //!< Outlet pressure [Pa]
    Series m_outletTemperature; //!< Outlet temperature [K]

    std::vector<Composition> m_inletComposition; //!< Inlet composition [fraction]
    std::vector<Composition> m_outletComposition; //!< Outlet composition [fraction]

    /*!
     * \brief Load boundary conditions from matrix (private)
     * \param bc Arma matrix with boundary conditions
     */
    void loadFromMatrix(const arma::mat& bc);
};
