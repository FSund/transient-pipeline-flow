#include <iostream>
#include <string>

#include "transflow.hpp"

using std::cout;
using std::endl;

int EP2();

int main()
{
    return EP2();
}

Pipeline makeEP2()
{
    arma::mat pipedata;
    if (!pipedata.load("D:/Simulations/EP2_heat_capacity/pipedata_with_height.csv"))
        throw std::runtime_error("File not loaded");

    const arma::uword N = 658; // number of grid points
    const double length = arma::sum(pipedata.col(8))*1000;
    cout << "length: " << length << endl;
    Pipeline pipeline(N, length);

    const arma::vec points = pipeline.gridPoints();
    const arma::vec loc = pipedata.col(0)*1000;
    pipeline.height() = utils::LinearInterpolator::getValuesAtPoints(loc, pipedata.col(2), points, 1);
    pipeline.burialDepth() = utils::LinearInterpolator::getValuesAtPoints(loc, pipedata.col(9), points, 1);

    pipeline.diameter().fill(1.016);
    pipeline.roughness().fill(0.485e-6);

    for (auto& bm : pipeline.burialMedium())
    {
        const double density = pipedata(0, 10);
        const double conductivity = pipedata(0, 11);
        const double heatCapacity = pipedata(0, 12);
        bm = BurialMedium(conductivity, density, heatCapacity);
    }

    for (auto& fluid : pipeline.ambientFluid())
    {
        const double velocity = pipedata(0, 13);
        const double viscosity = pipedata(0, 14);
        fluid = AmbientFluid(velocity, viscosity, Material::seawater);
    }

    // set up pipeline wall
    const arma::vec layer1thickness = utils::LinearInterpolator::getValuesAtPoints(loc, pipedata.col(1)/1e3, points, 1);
    const arma::vec layer2thickness = utils::LinearInterpolator::getValuesAtPoints(loc, pipedata.col(16)/1e3, points, 1);
    const arma::vec layer3thickness = utils::LinearInterpolator::getValuesAtPoints(loc, pipedata.col(20)/1e3, points, 1);
    const arma::vec layer1density = utils::LinearInterpolator::getValuesAtPoints(loc, pipedata.col(3), points, 1);
    const arma::vec layer2density = utils::LinearInterpolator::getValuesAtPoints(loc, pipedata.col(17), points, 1);
    const arma::vec layer3density = utils::LinearInterpolator::getValuesAtPoints(loc, pipedata.col(21), points, 1);
    const arma::vec layer1conductivity = utils::LinearInterpolator::getValuesAtPoints(loc, pipedata.col(4), points, 1);
    const arma::vec layer2conductivity = utils::LinearInterpolator::getValuesAtPoints(loc, pipedata.col(18), points, 1);
    const arma::vec layer3conductivity = utils::LinearInterpolator::getValuesAtPoints(loc, pipedata.col(22), points, 1);
    const arma::vec layer1heatCapacity = utils::LinearInterpolator::getValuesAtPoints(loc, pipedata.col(5), points, 1);
    const arma::vec layer2heatCapacity = utils::LinearInterpolator::getValuesAtPoints(loc, pipedata.col(19), points, 1);
    const arma::vec layer3heatCapacity = utils::LinearInterpolator::getValuesAtPoints(loc, pipedata.col(23), points, 1);
    for (size_t i = 0; i < pipeline.size(); i++)
    {
        std::vector<PipeWall::Layer> layers;
        layers.push_back(PipeWall::Layer(layer1thickness(i), layer1conductivity(i), layer1density(i), layer1heatCapacity(i)));
        layers.push_back(PipeWall::Layer(layer2thickness(i), layer2conductivity(i), layer2density(i), layer2heatCapacity(i)));
        layers.push_back(PipeWall::Layer(layer3thickness(i), layer3conductivity(i), layer3density(i), layer3heatCapacity(i)));

        pipeline.pipeWall().at(i) = PipeWall(layers);
    }

    return pipeline;
}

int EP2()
{
    Pipeline pipeline = makeEP2();

    TimeSeries bc(std::string("D:/Simulations/franpipe_oneyear/bc.csv"), 10500, 0);
    bc.setBoundarySettings({"inlet", "outlet", "inlet"});

    cout << bc.at(0) << endl;

    // initialize pipeline
    pipeline.flow().fill(bc.at(0).inletFlow());
    pipeline.pressure() = arma::linspace(bc.at(0).inletPressure(), bc.at(0).outletPressure(), pipeline.size());
    pipeline.temperature().fill(bc.at(0).inletTemperature());

    Config config;

    // make thermalized state first
    {
        arma::uvec timestamps = arma::linspace<arma::uvec>(0, 2*24*60*60, 24*12+1);
        TimeSeries therm(timestamps, std::vector<BoundaryConditions>(timestamps.size(), bc.at(0)));
        Simulator sim(pipeline, config);
        sim.simulate(therm);
        pipeline = sim.pipeline();
        pipeline.timestamp() = 0;

        cout << "Done thermalizing" << endl;
    }

//    pipeline.timestamp() = 3147600;
//    config.MVA = true;
//    config.samplingInterval = 5*60;
//    config.outputPath = "D:/Simulations/franpipe_oneyear/results/online_mva_p2_const/";

//    Simulator sim(pipeline, config);
//    sim.sampler().setIndicesToSample({0, sim.pipeline().size() - 1});

//    sim.simulate(bc);

    return 0;
}
