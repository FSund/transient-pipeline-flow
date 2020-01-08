#include <iostream>

#include "transflow.hpp"

int main()
{
    const arma::uword N = 10; // 10 grid points
    const double length = 10e3; // 10 km
    Pipeline pipeline(N, length);
    pipeline.roughness().fill(1e-6); // add some friction
    pipeline.ambientTemperature().fill(273.15 + 4); // 4 Celsius

    // initialize pipeline
    pipeline.pressure().fill(10e6); // 100 bara
    pipeline.flow().fill(100); // 100 kg/s
    pipeline.temperature().fill(280); // 280 Kelvin

    Config config;
    config.outputPath = "./output/";
    config.equationOfState = "BWRS";

    Simulator sim(pipeline, config);

    const arma::uword nSteps = 100;
    const arma::uword dt = 60;
    TimeSeries bc(nSteps, dt); // 100 steps of 60 seconds
    bc.inletFlow().fill(100); // 100 kg/s
    bc.outletPressure().fill(10e6); // 100 bara
    bc.inletTemperature().fill(280); // 280 Kelvin

    arma::vec nIterations = sim.simulate(bc);

    std::cout << sim.state().flow() << std::endl;
    std::cout << sim.state().pressure() << std::endl;
    std::cout << sim.state().temperature() << std::endl;
    std::cout << nIterations << std::endl;
}
