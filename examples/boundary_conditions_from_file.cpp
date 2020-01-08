#include <iostream>

#include "transflow.hpp"

int main()
{
    // load boundary conditions from file
    const std::string path = std::string(TRANSFLOW_RESOURCE_PATH) + "/examples/bc-with_composition.csv";
    const arma::uword lastRow = 1000; // last row to load (zero-indexed)
    TimeSeries bc(path, lastRow);

    std::cout << "Number of time steps: " << bc.size() << std::endl;
    std::cout << "Timestamps from " << bc.timestamps()(0) << " s to " << bc.timestamps().tail(1)(0) << " s" << std::endl;
    std::cout << "Total time: " << double(bc.timestamps()(1) - bc.timestamps()(0))/60.0 << " minutes" << std::endl;

    const arma::uword N = 800;
    Pipeline pipeline(N, 800e3); // 100 grid points, 800 km pipeline
    pipeline.roughness().fill(1e-5); // add some friction
    pipeline.constantComposition() = false;
    pipeline.ambientTemperature().fill(constants::kelvin + 5);

    // initialize pipeline
    pipeline.flow().fill(bc.inletFlow()(0));
    pipeline.pressure() = arma::linspace(bc.inletPressure()(0), bc.outletPressure()(0), N);
    pipeline.temperature().fill(bc.inletTemperature()(0));

    Config config;
    config.outputPath = "./output/";

    Simulator sim(pipeline, config);
    sim.simulate(bc);

    std::cout << "Inlet flow:         " << sim.state().flow()(0) << std::endl;
    std::cout << "Outlet flow:        " << sim.state().flow().tail(1)(0) << std::endl;
    std::cout << "Inlet pressure:     " << sim.state().pressure()(0) << std::endl;
    std::cout << "Outlet pressure:    " << sim.state().pressure().tail(1)(0) << std::endl;
    std::cout << "Inlet temperature:  " << sim.state().temperature()(0) << std::endl;
    std::cout << "Outlet temperature: " << sim.state().temperature().tail(1)(0) << std::endl;
    std::cout << "Inlet composition:  " << sim.state().composition().front().vec().t();
    std::cout << "Outlet composition: " << sim.state().composition().back().vec().t();
}
