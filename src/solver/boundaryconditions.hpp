#pragma once

#include <ostream>
#include <vector>
#include <string>
#include <armadillo>

#include "composition.hpp"

class Pipeline;

/*!
 * \brief The BoundaryConditions class is a container for the boundary
 * conditions at the inlet and outlet of a pipeline, for a single time step.
 *
 * The class uses a
 * <a href="http://arma.sourceforge.net/docs.html#adv_constructors_mat">armadillo fixed size matrix</a>
 * to guarantee that the boundary conditions have the correct format
 * (3 properties at the inlet and the outlet)
 * as well as the Composition class, which uses the same trick to ensure the
 * correct number of components (10).
 *
 * This class also has some convenience getters like inletFlow() etc., and a
 * std::ostream& operator << overload to allow nice printing of a
 * BoundaryConditions object.
 */
class BoundaryConditions
{
public:
    /*!
     * \brief BoundaryConditions with no properties set.
     */
    BoundaryConditions();

    /*!
     * \brief Constructor with defaults for inlet and outlet composition.
     *
     * This is marked explicit, to not allow implicit conversion between
     * arma::mat and BoundaryConditions.
     *
     * \param boundaryConditionMat Boundary conditions in form of a matrix
     * \param inletComposition Inlet composition
     * \param outletComposition Outlet composition
     */
    explicit BoundaryConditions(
            const arma::mat& boundaryConditionMat,
            const Composition& inletComposition = Composition::defaultComposition,
            const Composition& outletComposition = Composition::defaultComposition);

    /*!
     * \brief Construct from arma::vec for each property at inlet and outlet,
     * and composition.
     * \param inletFlow Inlet flow [kg/s]
     * \param outletFlow Outlet flow [kg/s]
     * \param inletPressure Inlet pressure [Pa]
     * \param outletPressure Outlet pressure [Pa]
     * \param inletTemperature Inlet temperature [K]
     * \param outletTemperature Outlet temperature [K]
     * \param inletComposition Inlet composition
     * \param outletComposition Outlet composition
     */
    BoundaryConditions(
            const double inletFlow,
            const double outletFlow,
            const double inletPressure,
            const double outletPressure,
            const double inletTemperature,
            const double outletTemperature,
            const Composition& inletComposition = Composition::defaultComposition,
            const Composition& outletComposition = Composition::defaultComposition);

    /*!
     * \brief Construct from a Pipeline object.
     *
     * \param pipeline Pipeline instance to get boundary conditions from
     * \param boundarySettings Boundary settings to apply. See
     * setBoundarySettings() examples for valid input.
     */
    explicit BoundaryConditions(
            const Pipeline& pipeline,
            const std::vector<std::string>& boundarySettings = {"inlet", "outlet", "inlet"});

    /*!
     * \brief Set boundary settings via vector of string, or brace-init-list
     * like `{"inlet", "outlet", "inlet"}`
     * \param strings Vector of strings or brace-init-list, like `{"inlet", "outlet", "inlet"}`.
     *
     * Valid options are `"inlet"`, `"outlet"`, `"both"` and `"none"`.
     */
    void setBoundarySettings(const std::vector<std::string>& strings);

    //! Returns the number of active boundary conditions.
    arma::uword nActiveBoundaryConditions() const;

    // getters
    inline const auto& inletFlow() const { return m_inletFlow; } //!< Get inlet flow [kg/s]
    inline const auto& outletFlow() const { return m_outletFlow; } //!< Get outlet flow [kg/s]

    inline const auto& inletPressure() const { return m_inletPressure; } //!< Get inlet pressure [Pa]
    inline const auto& outletPressure() const { return m_outletPressure; } //!< Get outlet pressure [Pa]

    inline const auto& inletTemperature() const { return m_inletTemperature; } //!< Get inlet temperature [K]
    inline const auto& outletTemperature() const { return m_outletTemperature; } //!< Get outlet temperature [k]

    inline const Composition& inletComposition() const { return m_inletComposition; } //!< Get inlet Composition
    inline const Composition& outletComposition() const { return m_outletComposition; } //!< Get outlet Composition

    // setters
    inline auto& inletFlow() { return m_inletFlow; } //!< Get inlet flow [kg/s]
    inline auto& outletFlow() { return m_outletFlow; } //!< Get outlet flow [kg/s]

    inline auto& inletPressure() { return m_inletPressure; } //!< Get inlet pressure [Pa]
    inline auto& outletPressure() { return m_outletPressure; } //!< Get outlet pressure [Pa]

    inline auto& inletTemperature() { return m_inletTemperature; } //!< Get inlet temperature [K]
    inline auto& outletTemperature() { return m_outletTemperature; } //!< Get outlet temperature [k]

    inline Composition& inletComposition() { return m_inletComposition; } //!< Inlet Composition setter
    inline Composition& outletComposition() { return m_outletComposition; } //!< Outlet Composition setter

    //! Get (const ref) inlet property from index
    //! 0 == flow
    //! 1 == pressure
    //! 2 == temperature
    inline const auto& inlet(const arma::uword i) const
    {
        if (i == 0)
            return inletFlow();
        else if (i == 1)
            return inletPressure();
        else if (i == 2)
            return inletTemperature();
        else
            throw std::runtime_error("invalid index i");
    }

    //! Get (const ref) outlet property from index
    //! 0 == flow
    //! 1 == pressure
    //! 2 == temperature
    inline const auto& outlet(const arma::uword i) const
    {
        if (i == 0)
            return outletFlow();
        else if (i == 1)
            return outletPressure();
        else if (i == 2)
            return outletTemperature();
        else
            throw std::runtime_error("invalid index i");
    }

    /*!
     * \brief Overload the ostream::operator<< so we can do "cout << bc" and get
     * nice print of the boundary conditions.
     */
    friend std::ostream& operator << (std::ostream &out, const BoundaryConditions& c);

    //!
    //! \brief The SingleCondition class is used to store a single boundary
    //! condition (like flow, pressure or temperature), and also contains
    //! information on whether the condition is "active", meaning if it is
    //! meant to be a constrain when solving the governing equations.
    //!
    class SingleCondition
    {
    public:
        friend class BoundaryConditions;

        //!
        //! \brief Construct from value and isActive bool
        //! \param value Value at boundary
        //! \param active If the condition is active
        //!
        SingleCondition(const double value, const bool active):
            m_value(value),
            m_isActive(active)
        {}

        //! User-defined conversion to double, which allows implicit conversion.
        inline operator double() const { return m_value; }

        //! User-defined conversion to double, which allows implicit conversion.
        inline operator double() { return m_value; }

        //! Get the value
        double value() const { return m_value; }

        //! Get the value
        double value() { return m_value; }

        //! If the condition is active or not
        bool isActive() const { return m_isActive; }

        //! If the condition is active or not
        bool isActive() { return m_isActive; }

    private:
        //! Set the active state of the condition
        void setActive(const bool active) { m_isActive = active; }

        double m_value; //!< The value at the boundary
        bool m_isActive; //!< If the condition is active or not
    };

    /*!
     * \brief Construct BoundaryConditions::SingleCondition for each property at
     * inlet and outlet, and composition.
     * \param inletFlow Inlet flow [kg/s]
     * \param outletFlow Outlet flow [kg/s]
     * \param inletPressure Inlet pressure [Pa]
     * \param outletPressure Outlet pressure [Pa]
     * \param inletTemperature Inlet temperature [K]
     * \param outletTemperature Outlet temperature [K]
     * \param inletComposition Inlet composition
     * \param outletComposition Outlet composition
     */
    BoundaryConditions(
            const SingleCondition inletFlow,
            const SingleCondition outletFlow,
            const SingleCondition inletPressure,
            const SingleCondition outletPressure,
            const SingleCondition inletTemperature,
            const SingleCondition outletTemperature,
            const Composition& inletComposition = Composition::defaultComposition,
            const Composition& outletComposition = Composition::defaultComposition);

private:
    SingleCondition m_inletFlow; //!< Inlet flow [kg/s]
    SingleCondition m_inletPressure; //!< Inlet pressure [Pa]
    SingleCondition m_inletTemperature; //!< Inlet temperature [K]
    SingleCondition m_outletFlow; //!< Outlet flow [kg/s]
    SingleCondition m_outletPressure; //!< Outlet presure [Pa]
    SingleCondition m_outletTemperature; //!< Outlet temperature [K]

    Composition m_inletComposition; //!< Inlet Composition [fraction]
    Composition m_outletComposition; //!< Outlet Composition [fraction]

    //! Set a single boundary setting (private)
    void setBoundarySettings(const arma::uword i, const std::string& setting)
    {
        if (setting == "none")
        {
            this->inlet(i).setActive(false);
            this->outlet(i).setActive(false);
        }
        else if (setting == "inlet")
        {
            this->inlet(i).setActive(true);
            this->outlet(i).setActive(false);
        }
        else if (setting == "outlet")
        {
            this->inlet(i).setActive(false);
            this->outlet(i).setActive(true);
        }
        else if (setting == "both")
        {
            this->inlet(i).setActive(true);
            this->outlet(i).setActive(true);
        }
        else
        {
            throw std::runtime_error("invalid setting \"" + setting + "\"");
        }
    }

    //! Get (non-const ref) inlet property from index (private)
    //! 0 == flow
    //! 1 == pressure
    //! 2 == temperature
    inline SingleCondition& inlet(const arma::uword i)
    {
        if (i == 0)
            return inletFlow();
        else if (i == 1)
            return inletPressure();
        else if (i == 2)
            return inletTemperature();
        else
            throw std::runtime_error("invalid index i");
    }

    //! Get (non-const ref) outlet property from index (private)
    //! 0 == flow
    //! 1 == pressure
    //! 2 == temperature
    inline SingleCondition& outlet(const arma::uword i)
    {
        if (i == 0)
            return outletFlow();
        else if (i == 1)
            return outletPressure();
        else if (i == 2)
            return outletTemperature();
        else
            throw std::runtime_error("invalid index i");
    }
};

/*!
 * \brief The BoundaryConditionsStamped class is a subclass of BoundaryConditions with an extra
 * member to store the timestamp of the boundary conditions.
 */
class BoundaryConditionsStamped : public BoundaryConditions
{
public:
    /*!
     * \brief Construct from timestamp and BoundaryConditions instance.
     * \param timestamp Timestamp [s]
     * \param boundaryConditions BoundaryConditions instance
     */
    BoundaryConditionsStamped(
        const arma::uword timestamp,
        const BoundaryConditions &boundaryConditions) : BoundaryConditions(boundaryConditions),
                                                        m_timestamp(timestamp)
    {}

    //! Get (copy of) timestamp
    arma::uword timestamp() const { return m_timestamp; }

private:
    arma::uword m_timestamp; //!< Timestamp [s]
};

