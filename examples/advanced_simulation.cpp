#include <iostream>

#include "transflow.hpp"

// Example using most of the features of Pipeline and Simulator
int main()
{
    // set up pipeline
    const arma::uword N = 10;
    Pipeline pipeline(N, 100e3); // 10 grid points, 100 km pipeline
    pipeline.roughness().fill(1e-5); // add some friction
    pipeline.diameter().fill(0.8);
//    pipeline.height() = arma::linspace(0, 100, N);
    pipeline.burialDepth().fill(2.0);

    // initialize pipeline
    pipeline.flow().fill(0); // start with no flow
    pipeline.pressure().fill(10e6); // 100 bar
    pipeline.temperature().fill(280);

    //// heat transfer stuff
    // custom burial medium with conductivity of 4.0, density of 2500 kg/m3,
    // and heat capacity of 700  [J/(kg K)]
    BurialMedium burialMedium(4.0, 2500, 700);
    pipeline.burialMedium() = std::vector<BurialMedium>(N, burialMedium);

    // custom ambient fluid -- sea water with velocity 0.2 m/s and increased viscosity
    AmbientFluid ambientFluid(0.2, 1.1/1000.0, Material::seawater);
    pipeline.ambientFluid() = std::vector<AmbientFluid>(N, ambientFluid);

    // custom pipe wall
    PipeWall pipeWall = PipeWall::defaultPipeWall;
    pipeWall.layer(2).thickness() = 0.12; // custom concrete thickness
    pipeline.pipeWall() = std::vector<PipeWall>(N, pipeWall);

    // increase ambient temperature
    pipeline.ambientTemperature() += 5;

    // set up simulator
    Config config;
    config.outputPath = "./output/";
    config.equationOfState = "GERG04";
    config.heatTransfer = "Unsteady";
    Simulator sim(pipeline, config); // 10 grid points

    // simulate with slowly increasing flow from 0 to 100 kg/s in 20 minutes
    const arma::uword dt = 60; // 1 minute time step
    arma::uword nSteps = 20;
    TimeSeries bc(nSteps, dt); // 20 steps of 60 seconds
    bc.inletFlow() = arma::linspace(0, 100, bc.size());
    bc.outletPressure().fill(10e6);
    bc.inletTemperature().fill(pipeline.temperature()(0));
    sim.simulate(bc);

    std::cout << "After ramp-up" << std::endl;
    std::cout << sim.state().flow() << std::endl;
    std::cout << sim.state().pressure() << std::endl;
    std::cout << sim.state().temperature() << std::endl;

    // simulate with constant boundary conditions for 12*12 steps of
    // 5*60 seconds (12 hours total)
    bc = TimeSeries(5*60, std::vector<BoundaryConditions>(12*12, bc.at(bc.size() - 1)));
    bc.timestamps() += sim.pipeline().timestamp() + dt;
    sim.simulate(bc);

    std::cout << "After 12 hours" << std::endl;
    std::cout << sim.state().flow() << std::endl;
    std::cout << sim.state().pressure() << std::endl;
    std::cout << sim.state().temperature() << std::endl;

    // simulate zero inlet flow and some outlet flow
    nSteps = 5*60;
    bc = TimeSeries(sim.pipeline(), nSteps, dt);
    bc.timestamps() += sim.pipeline().timestamp() + dt;
    bc.inletFlow().fill(0); // close off inlet
    bc.outletFlow().fill(5); // slowly drain via outlet
    bc.setBoundarySettings({"both", "outlet", "inlet"});

    sim.simulate(bc);

    std::cout << "After zero inlet flow and some outlet flow" << std::endl;
    std::cout << sim.state().flow() << std::endl;
    std::cout << sim.state().pressure() << std::endl;
    std::cout << sim.state().temperature() << std::endl;
}
