#pragma once

#include <armadillo>

#include <vector>

#include "composition.hpp"
#include "heattransfer/pipewall.hpp"
#include "heattransfer/burialmedium.hpp"
#include "heattransfer/ambientfluid.hpp"
#include "heattransfer/heattransferbase.hpp"
#include "advection/batchtrackingstate.hpp"

class BoundaryConditions;

class Pipeline
{
public:
    ////////////////////////////////////////////////////////////////////////////////
    class State
    {
    public:
        friend class Pipeline;

        State(
                const arma::vec& gridPoints,
                const arma::vec& pressure,
                const arma::vec& temperature,
                const arma::vec& flow,
                const std::vector<Composition>& composition);

        State(
                const arma::vec& gridPoints,
                const arma::vec& pressure,
                const arma::vec& temperature,
                const arma::vec& flow,
                const Composition& composition = Composition::defaultComposition);

        // getters
        // main/governing properties
        const arma::vec& flow() const                           { return m_flow; }
        const arma::vec& pressure() const                       { return m_pressure; }
        const arma::vec& temperature() const                    { return m_temperature; }
        const std::vector<Composition>& composition() const     { return m_composition; }

        // derived properties
        const arma::vec& heatCapacityConstantVolume() const     { return m_heatCapacityConstantVolume; }
        const arma::vec& heatCapacityConstantPressure() const   { return m_heatCapacityConstantPressure; }
        const arma::vec& density() const                        { return m_density; }
        const arma::vec& viscosity() const                      { return m_viscosity; }
        const arma::vec& specificGasConstant() const            { return m_specificGasConstant; }
        const arma::vec& molarMass() const                      { return m_molarMass; }

        const arma::vec& compressibilityFactor() const          { return m_compressibilityFactor; }
        const arma::vec& dZdtAtConstantPressure() const         { return m_temperatureDerivativeConstantPressure; }
        const arma::vec& dZdpAtConstantTemperature() const      { return m_pressureDerivativeConstantTemperature; }
        const arma::vec& dZdtAtConstantDensity() const          { return m_temperatureDerivativeConstantDensity; }

        const arma::vec& velocity() const                       { return m_velocity; }
        const arma::vec& frictionFactor() const                 { return m_frictionFactor; }
        const arma::vec& reynoldsNumber() const                 { return m_reynoldsNumber; }

        const arma::vec& ambientTemperature() const             { return m_ambientTemperature; }
        const arma::vec& heatFlow() const                       { return m_heatFlow; }

//        const std::vector<HeatTransferState>& heatTransferState() const { return m_heatTransferState; }
//        const BatchTrackingState& batchTrackingState() const  { return m_batchTrackingState; }
//        bool batchTrackingIsInitialized() const                 { return m_batchTrackingIsInitialized; }

//        // setters
//        // main/governing properties
//        arma::vec& flow()                           { return m_flow; }
//        arma::vec& pressure()                       { return m_pressure; }
//        arma::vec& temperature()                    { return m_temperature; }
//        std::vector<Composition>& composition();

//        // derived properties
//        arma::vec& heatCapacityConstantVolume()     { return m_heatCapacityConstantVolume; }
//        arma::vec& heatCapacityConstantPressure()   { return m_heatCapacityConstantPressure; }
//        arma::vec& density()                        { return m_density; }
//        arma::vec& viscosity()                      { return m_viscosity; }
//        arma::vec& specificGasConstant()            { return m_specificGasConstant; }
//        arma::vec& molarMass()                      { return m_molarMass; }

//        arma::vec& compressibilityFactor()          { return m_compressibilityFactor; }
//        arma::vec& dZdtAtConstantPressure()         { return m_temperatureDerivativeConstantPressure; }
//        arma::vec& dZdpAtConstantTemperature()      { return m_pressureDerivativeConstantTemperature; }
//        arma::vec& dZdtAtConstantDensity()          { return m_temperatureDerivativeConstantDensity; }

//        arma::vec& velocity()                       { return m_velocity; }
//        arma::vec& frictionFactor()                 { return m_frictionFactor; }
//        arma::vec& reynoldsNumber()                 { return m_reynoldsNumber; }

//        arma::vec& ambientTemperature()             { return m_ambientTemperature; }
//        arma::vec& heatFlow()                       { return m_heatFlow; }

//        std::vector<HeatTransferState>& heatTransferState() { return m_heatTransferState; }
//        BatchTrackingState& batchTrackingState()  { return m_batchTrackingState; }
//        bool& batchTrackingIsInitialized()          { return m_batchTrackingIsInitialized; }

        /*!
         * \brief For pretty printing.
         *
         * So we can do
         *
         *     std::cout << Pipeline::State;
         *
         * and get nice printing.
         */
        friend std::ostream& operator << (std::ostream &out, const State& c);

    protected:
        State(const arma::vec& gridPoints);

        // input/governing properties
        arma::vec m_flow; //!< Flow [kg/s]
        arma::vec m_pressure; //!< Gas pressure [Pa]
        arma::vec m_temperature; //!< Gas temperature [K]
        std::vector<Composition> m_composition; //!< Gas composition (fractions)

        // output/derived properties
        arma::vec m_heatCapacityConstantVolume; //!< Gas heat capacity at constant volume \f$c_v\f$ [J/(kg K)]
        arma::vec m_heatCapacityConstantPressure; //!< Gas heat capacity at constant pressure \f$c_p\f$ [J/(kg K)]
        arma::vec m_density; //!< Gas density [kg/m3]
        arma::vec m_viscosity; //!< Gas dynamic viscosity [Pa s] = [kg/m*s]
        arma::vec m_specificGasConstant; //!< Specific gas constant of the gas [J/(kg K)]
        arma::vec m_molarMass; //!< Molar mass of the gas [g/mol]

        arma::vec m_compressibilityFactor; //!< Compressibility factor [-]
        //! Partial derivative of the compressibility factor with respect to
        //! temperature, at constant pressure,
        //! \f$\frac{\partial Z}{\partial T}|_p\f$ [-]
        arma::vec m_temperatureDerivativeConstantPressure;
        //! Partial derivative of the compressibility factor with respect to
        //! pressure, at constant temperature,
        //! \f$\frac{\partial Z}{\partial p}|_T\f$ [-]
        arma::vec m_pressureDerivativeConstantTemperature;
        //! Partial derivative of the compressibility factor with respect to
        //! temperature, at constant density,
        //! \f$\frac{\partial Z}{\partial T}|_\rho\f$ [-]
        arma::vec m_temperatureDerivativeConstantDensity;

        arma::vec m_velocity; //!< Gas velocity [m/s]
        arma::vec m_frictionFactor; //!< Friction factor [-]
        arma::vec m_reynoldsNumber; //!< Reynolds number [-]

        arma::vec m_ambientTemperature; //!< Ambient temperature [K]
        arma::vec m_heatFlow; //!< Heat flow q [W/m2]

        std::vector<HeatTransferState> m_heatTransferState; //!< Heat transfer state
        bool m_heatTransferIsInitialized; //!< If heat transfer is initialized

        BatchTrackingState m_batchTrackingState; //!< Batch tracking state
        bool m_batchTrackingIsInitialized; //!< If batch tracking is initialized
    };
    ////////////////////////////////////////////////////////////////////////////////

    //!
    //! \brief Construct from number of grid points and pipeline length.
    //! \param nGridPoints Number of grid points
    //! \param length Total length [m]
    //!
    Pipeline(const arma::uword nGridPoints = 100, const double length = 100e3);

    arma::uword size() const                                { return m_gridPoints.n_elem; }

    //! Get current inlet and outlet state as BoundaryConditions instance.
    BoundaryConditions getBoundaryConditions() const;

    //! Enable batch tracking. Also initializes the batch tracking state.
    void enableBatchTracking();

    //! Initialize the batch tracking state.
    void initializeBatchTracking();

    //! Set the length of the pipeline. This updates the grid points and the
    //! batch tracking state.
    void setLength(const double length);

    arma::uword timestamp() const { return m_timestamp; }
    arma::uword& timestamp() { return m_timestamp; }

    const Pipeline::State& state() const { return m_state; }

    // getters
    // constants
    double length() const                                   { return m_length; }
    const arma::vec& gridPoints() const                     { return m_gridPoints; }
    const arma::vec& diameter() const                       { return m_diameter; }
    const arma::vec& height() const                         { return m_height; }
    const arma::vec& elevation() const                      { return m_height; }
    const arma::vec& roughness() const                      { return m_roughness; }

    const arma::vec& burialDepth() const                    { return m_burialDepth; }
    const std::vector<PipeWall>& pipeWall() const           { return m_pipeWall; }
    const std::vector<BurialMedium>& burialMedium() const   { return m_burialMedium; }
    const std::vector<AmbientFluid>& ambientFluid() const   { return m_ambientFluid; }

    // main/governing properties
    const arma::vec& flow() const                           { return prop().m_flow; }
    const arma::vec& pressure() const                       { return prop().m_pressure; }
    const arma::vec& temperature() const                    { return prop().m_temperature; }
    const std::vector<Composition>& composition() const     { return prop().m_composition; }

    // derived properties
    const arma::vec& heatCapacityConstantVolume() const     { return prop().m_heatCapacityConstantVolume; }
    const arma::vec& heatCapacityConstantPressure() const   { return prop().m_heatCapacityConstantPressure; }
    const arma::vec& density() const                        { return prop().m_density; }
    const arma::vec& viscosity() const                      { return prop().m_viscosity; }
    const arma::vec& specificGasConstant() const            { return prop().m_specificGasConstant; }
    const arma::vec& molarMass() const                      { return prop().m_molarMass; }

    const arma::vec& compressibilityFactor() const          { return prop().m_compressibilityFactor; }
    const arma::vec& dZdtAtConstantPressure() const         { return prop().m_temperatureDerivativeConstantPressure; }
    const arma::vec& dZdpAtConstantTemperature() const      { return prop().m_pressureDerivativeConstantTemperature; }
    const arma::vec& dZdtAtConstantDensity() const          { return prop().m_temperatureDerivativeConstantDensity; }

    const arma::vec& velocity() const                       { return prop().m_velocity; }
    const arma::vec& frictionFactor() const                 { return prop().m_frictionFactor; }
    const arma::vec& reynoldsNumber() const                 { return prop().m_reynoldsNumber; }

    const arma::vec& ambientTemperature() const             { return prop().m_ambientTemperature; }
    const arma::vec& heatFlow() const                       { return prop().m_heatFlow; }

    const std::vector<HeatTransferState>& heatTransferState() const { return prop().m_heatTransferState; }
    bool heatTransferIsInitialized() const                  { return prop().m_heatTransferIsInitialized; }

    const BatchTrackingState& batchTrackingState() const    { return prop().m_batchTrackingState; }
    bool batchTrackingIsInitialized() const                 { return prop().m_batchTrackingIsInitialized; }

    // settings
    bool constantComposition() const                        { return m_constantComposition; }

    // helpers
    //! Get inlet composition
    const arma::vec& inletComposition() const;
    //! Get outlet composition
    const arma::vec& outletComposition() const;

    // setters
    // constants
    arma::vec& gridPoints()                     { return m_gridPoints; }
    arma::vec& diameter()                       { return m_diameter; }
    arma::vec& height()                         { return m_height; }
    arma::vec& elevation()                      { return m_height; }
    arma::vec& roughness()                      { return m_roughness; }

    arma::vec& burialDepth()                    { return m_burialDepth; }
    std::vector<PipeWall>& pipeWall()           { return m_pipeWall; }
    std::vector<BurialMedium>& burialMedium()   { return m_burialMedium; }
    std::vector<AmbientFluid>& ambientFluid()   { return m_ambientFluid; }

    // main/governing properties
    arma::vec& flow()                           { return prop().m_flow; }
    arma::vec& pressure()                       { return prop().m_pressure; }
    arma::vec& temperature()                    { return prop().m_temperature; }

    // derived properties
    arma::vec& heatCapacityConstantVolume()     { return prop().m_heatCapacityConstantVolume; }
    arma::vec& heatCapacityConstantPressure()   { return prop().m_heatCapacityConstantPressure; }
    arma::vec& density()                        { return prop().m_density; }
    arma::vec& viscosity()                      { return prop().m_viscosity; }
    arma::vec& specificGasConstant()            { return prop().m_specificGasConstant; }
    arma::vec& molarMass()                      { return prop().m_molarMass; }

    arma::vec& compressibilityFactor()          { return prop().m_compressibilityFactor; }
    arma::vec& dZdtAtConstantPressure()         { return prop().m_temperatureDerivativeConstantPressure; }
    arma::vec& dZdpAtConstantTemperature()      { return prop().m_pressureDerivativeConstantTemperature; }
    arma::vec& dZdtAtConstantDensity()          { return prop().m_temperatureDerivativeConstantDensity; }

    arma::vec& velocity()                       { return prop().m_velocity; }
    arma::vec& frictionFactor()                 { return prop().m_frictionFactor; }
    arma::vec& reynoldsNumber()                 { return prop().m_reynoldsNumber; }

    arma::vec& ambientTemperature()             { return prop().m_ambientTemperature; }
    arma::vec& heatFlow()                       { return prop().m_heatFlow; }

    std::vector<HeatTransferState>& heatTransferState() { return prop().m_heatTransferState; }
    bool& heatTransferIsInitialized()           { return prop().m_heatTransferIsInitialized; }

    BatchTrackingState& batchTrackingState()    { return prop().m_batchTrackingState; }
    bool& batchTrackingIsInitialized()          { return prop().m_batchTrackingIsInitialized; }

    // settings
    bool& constantComposition()                 { return m_constantComposition; }

    /*!
     * \brief Set the composition. This also updates the batch tracking state to
     * reflect the composition change.
     * \param composition New composition at each grid point
     */
    void updateComposition(const std::vector<Composition>& composition);

    /*!
     * \brief Set the composition of the whole pipeline. This also updates the
     * batch tracking state to reflect the composition change.
     * \param composition New composition (same for the whole pipeline)
     */
    void updateComposition(const Composition& composition);

    /*!
     * \brief Set the composition without updating batch tracking state. This is
     * meant for internal use in Solver etc.
     * \param composition New composition at each grid point
     */
    void setCompositionUnsafe(const std::vector<Composition>& composition);

private:
    // constants
    double m_length; //!< Total length [m]
    arma::vec m_gridPoints; //!< Grid points [m]
    arma::vec m_diameter; //!< Inner diameter [m]
    arma::vec m_height; //!< Height profile/elevation [m]
    arma::vec m_roughness; //!< Sand grain equivalent roughness [m]

    arma::vec m_burialDepth; //!< Distance from burial medium to top of pipe [m]
    std::vector<PipeWall> m_pipeWall; //!< Pipe wall description
    std::vector<BurialMedium> m_burialMedium; //!< Burial medium description
    std::vector<AmbientFluid> m_ambientFluid; //!< Ambient fluid description

    // settings
    bool m_constantComposition; //!< If we simulate using constant composition or not

    // other
    State m_state; //!< Pipeline::State instance with all properties like flow etc.
    arma::uword m_timestamp = 0; //!< Current timestamp [s]

    // shorthand for internal use
    const State& prop() const { return m_state; }
    State& prop() { return m_state; }
};
