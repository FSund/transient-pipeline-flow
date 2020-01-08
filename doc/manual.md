# User guide {#user_guide}

\tableofcontents

## Typical usage

To set up a simulation of a pipeline, you typically start by creating a 
Pipeline instance which describes the pipeline.

    const double length = 10e3; // length
    const double N = 10; // number of grid points
    Pipeline pipeline(N, length);
    pipeline.roughness().fill(1e-6); // roughness in meters
    pipeline.ambientTemperature().fill(278); // 278 Kelvin
    pipeline.diameter().fill(0.7); // inner diameter of 70 cm
    /* pipeline wall, burial depth, burial medium, ambient fluid, etc. */

(See the members of Pipeline to see which other properties can be set.)

Then you have to initialize the flow, pressure and temperature of the pipeline. 
This should have some relation to the boundary condiditions that will be used 
during the simulation step later, or the program will have issues with
convergence.

    pipeline.flow().fill(0); // no initial flow
    pipeline.pressure().fill(10e6); // 100 bar
    pipeline.temperature().fill(280); // 280 Kelvin

You then create a Config instance, which has different options for the program

    Config config;
    config.equationOfState = "GERG04";
    config.outputPath = "./results/";
    config.samplingInterval = 5*60; // sample every 5 minutes
    /*...*/

Before we can simulate we need some boundary conditions

    TimeSeries bc(100, 60); // 100 steps of 60 seconds
    bc.inletFlow().fill(100); // 100 kg/s
    bc.outletPressure().fill(10e6); // 100 bar
    bc.inletTemperature().fill(280); // 280 Kelvin

Finally create a Simulator instance from the Pipeline and Config, and perform 
the simulation

    Simulator simulator(pipeline, config);
    simulator.simulate(bc);

The results are then stored to the folder in Config::outputPath, and the final 
state can also be inspected via

    std::cout << simulator.state().flow() << std::endl;
    std::cout << simulator.state().pressure() << std::endl;
    std::cout << simulator.state().temperature() << std::endl;
    /*...*/

See the files located in the `examples` for working examples of how to set up a simulation.

## Advanced options

### Sampler

Which results to save can be configured via the Sampler member of Simulator,
accessible via Simulator::sampler(). Flow, pressure, temperature and outlet 
composition are stored by default.

Add properties via Sampler::addPropertyToPrint, which takes in getters member 
functions of Pipeline, as follows

    simulator.sampler().addPropertyToPrint(&Pipeline::reynoldsNumber);
    simulator.sampler().addPropertyToPrint(&Pipeline::density);
    /*...*/

See the documentation of Pipeline to see which getters/properties are
available.

The convenience functions Pipeline::inletComposition() and 
Pipeline::outletComposition() have been defined, and can also be stored as 
follows

    simulator.sampler().addPropertyToPrint(&Pipeline::inletComposition);
    simulator.sampler().addPropertyToPrint(&Pipeline::outletComposition);

### Boundary conditions

Boundary conditions for simulating the pipeline are stored in TimeSeries 
instances. This consists of a series of timestamps, and flow, pressure, 
temperature and composition at the inlet and outlet of the pipeline at each 
time stamp, and also information regarding which property is a constraint at 
the inlet/outlet when solving the governing equations.

TimeSeries can be created in many different ways, but perhaps the most useful 
one is by loading a CSV-file, via 
TimeSeries::TimeSeries(const std::string&, const std::vector<std::string>&)

    TimeSeries bc("boundaryConditions.csv");

When using this constructor the boundary settings have to be specified, either
via an extra argument to the constructor

    TimeSeries bc("boundaryConditions.csv", {"inlet", "outlet", "intlet"});

or via a call to TimeSeries::setBoundarySettings() after construction

    bc.setBoundarySettings({"inlet", "outlet", "intlet"});

Check the documentation of 
TimeSeries::TimeSeries(const std::string&, arma::uword, arma::uword, const std::vector<std::string>&)
to see how the CSV-file should be structured and what units are expected.

Other useful examples is constructing an empty TimeSeries and setting all 
properties manually

    // create empty TimeSeries
    const int N = 100; // number of steps
    const double dt = 60; // 60 seconds
    TimeSeries bc(N, 60);

    // set properties
    bc.inletFlow() = arma::linspace(0, 100, N); // overloaded copy-assignment operator
    bc.outletPressure().fill(10e6);
    bc.inletTemperature().fill(280);

When we use this method the properties that are set via fill() or assignment 
are automatically marked to be used as boundary conditions. For example, after
calling

    bc.inletFlow().fill(100);

`bc.inletFlow().isActive()` will return `true`, and the inlet flow will be used
as a boundary condition when solving the governing equations.

The composition is stored as `std::vector<Composition>`, so it can be a bit 
more fiddly to set up, but can be done for example via

    bc.inletComposition() = std::vector<Composition>(N, Composition::defaultComposition);

### Equation of state
At the moment there are two different equations of state implemented in the 
application:
 - The <a href="https://en.wikipedia.org/wiki/Benedict%E2%80%93Webb%E2%80%93Rubin_equation#The_BWRS_equation_of_state">BWRS</a> (Benedict–Webb–Rubin-Starling) equation of state, in the class BWRS
 - The <a href="http://www.gerg.eu/public/uploads/files/publications/technical_monographs/tm15_04.pdf">GERG-2004</a> Wide-Range Equation of State for Natural Gases and Other Mixtures, in the class GERG04

This is selected via Config::equationOfState.

    config.equationOfState = "BWRS";

or

    config.equationOfState = "GERG04";

### Heat transfer

Four different heat transfer models are implemented:
 - Unsteady 1d radial heat transfer, in the class UnsteadyHeatTransfer
 - Steady state heat transfer, in the class SteadyStateHeatTransfer
 - Fixed Q value (steady state), in the class FixedQValue
 - Fixed U value (steady state), in the class FixedUValue

This is selected via Config::heatTransfer, for example

    config.equationOfState = "SteadyState";

The properties governing heat transfer can be controlled via the Pipeline 
instance passed to the Simulator constructor. This is mainly the inner 
diameter, the PipeWall, burial depth, BurialMedium, and AmbientFluid members of 
Pipeline.

Heat transfer is implemented with one HeatTransferBase instance at each grid 
point.  So in theory it is possible to have different heat transfer models at 
each grid point.  But this is not exposed via any interface at the moment (but 
should not be difficult to implement).

**NB:** At the moment there are no good ways of controlling the Q-value and U-value of 
FixedQValue and FixedUValue. When a config is passed to the Simulator 
constructor, the config is passed on to the 
Physics::Physics(const Pipeline&, const Config&) constructor, and 
Config::heatTransfer is passed on to finally reach HeatTransfer::makeSingle 
(not documented), which just sets Q and U to zero. It *can* be set as follows, 
but the interface is not very user-friendly, and it makes use of mutable member 
variables, which is kind of an anti-pattern

    dynamic_cast<const FixedQValue&>(sim.physics().heatTransfer().at(0)).setQValue(q);

(`FixedQValue&` needs to be `const` as the `.physics()` getter returns 
`const Physics&`, and that is why mutable members are used -- and it is
preferred to not introduce non-const getters for Simulator::m_physics etc.)

### Equation solver
The Solver class takes care of solving the governing equations.
This has several configuration options, which are documented below.

#### Convergence and relaxation
The equation solver uses an iterative approach, which checks for convergence
after each iteration. The tolerance used when checking if the system has
converged is controlled via Config::tolerances

    config.tolerances = {0.001, 0.001, 0.001}; // {flow, pressure, temperature}

The convergence check itself is performed in Solver::differencesWithinTolerance.

The tolerance can be configured as either a *relative* or an *absolute* limit,
via Config::toleranceType as follows

    config.toleranceType = "relative";

or

    config.toleranceType = "absolute";

"relative" is the most commonly used option. The maximum number of iterations
that are performed before returning is controlled via Config::maxIterations

    config.maxIterations = 200;

Relaxation factors are also implemented to allow for easier convergence.
The default values are 1.0 for flow and pressure, and 2/3 for temperature.

    config.relaxationFactors = {1, 1, 2/3.0}; // {flow, pressure, temperature}

There is also an option to force a given number of iterations, skipping the
convergence check. This is enabled via 

    config.bruteForce = true;

and the number of iterations is controlled via

    config.maxIterations = 10;

#### Energy equation
Two different types of energy equations are implemented;
the internal energy form, and the enthalpy form (see *Form of energy equation 
in gas-pipeline simulations* (Filip Sund and Tor Ytrehus, 2018) for more info).
This can be selected via the Config::discretizer parameter.

    std::string discretizer = "InternalEnergy";

or

    std::string discretizer = "Enthalpy";
