#include "debug.hpp"

#include <functional>
#include <filesystem>

#include "sampler.hpp"
#include "pipeline.hpp"

using namespace arma;
using namespace std;

typedef const arma::vec& (Pipeline::*PipelineMemFun)() const; // CONST!

TEST_CASE("Sampler")
{
    auto dir = std::filesystem::current_path() / "output";
    Sampler sampler(dir);

    // check that adding properties creates files in the correct place
    // first remove any existing file
    std::filesystem::remove(dir / "heatFlow.csv");
    CHECK(!std::filesystem::exists(dir / "heatFlow.csv"));
    // then add the property
    sampler.addPropertyToPrint(&Pipeline::heatFlow);
    CHECK(std::filesystem::exists(dir / "heatFlow.csv"));

    // check that sampling works as expected
    Pipeline pipeline(5);
    pipeline.heatFlow() = {1, 2, 3, 4, 10};
    // sample to file
    sampler.sample(pipeline);
    // load file and compare contents
    mat data;
    data.load((dir / "heatFlow.csv").string(), arma::csv_ascii);
    vec heatFlow = data.row(0).t();
    CHECK(equal(heatFlow, vec({0, 1, 2, 3, 4, 10}))); // first zero is timestamp
}
TEST_CASE("interval")
{
    // check that Sampler only samples at the given intervals
    auto dir = std::filesystem::current_path() / "output";
    Sampler sampler(dir, 240);

    Pipeline pipeline(5);
    pipeline.timestamp() = 1000;
    CHECK(sampler.sample(pipeline));
    CHECK_FALSE(sampler.sample(pipeline)); // repeated sample at same timestamp
    CHECK_FALSE(sampler.sample(pipeline)); // repeated sample at same timestamp

    pipeline.timestamp() += 239;
    CHECK_FALSE(sampler.sample(pipeline));
    CHECK_FALSE(sampler.sample(pipeline));

    pipeline.timestamp() += 1;
    CHECK(sampler.sample(pipeline));
    CHECK_FALSE(sampler.sample(pipeline));
}

TEST_CASE("some functions")
{
    SUBCASE("getSampleLabel")
    {
        CHECK(Sampler::getSampleLabel(&Pipeline::flow) == "flow");
        CHECK(Sampler::getSampleLabel(&Pipeline::pressure) == "pressure");
        CHECK(Sampler::getSampleLabel(&Pipeline::temperature) == "temperature");
        CHECK(Sampler::getSampleLabel(&Pipeline::heatCapacityConstantVolume) == "heatCapacityConstantVolume");
        CHECK(Sampler::getSampleLabel(&Pipeline::heatCapacityConstantPressure) == "heatCapacityConstantPressure");
        CHECK(Sampler::getSampleLabel(&Pipeline::density) == "density");
        CHECK(Sampler::getSampleLabel(&Pipeline::viscosity) == "viscosity");
        CHECK(Sampler::getSampleLabel(&Pipeline::specificGasConstant) == "specificGasConstant");
        CHECK(Sampler::getSampleLabel(&Pipeline::molarMass) == "molarMass");
        CHECK(Sampler::getSampleLabel(&Pipeline::compressibilityFactor) == "compressibilityFactor");
        CHECK(Sampler::getSampleLabel(&Pipeline::dZdtAtConstantPressure) == "dZdtAtConstantPressure");
        CHECK(Sampler::getSampleLabel(&Pipeline::dZdpAtConstantTemperature) == "dZdpAtConstantTemperature");
        CHECK(Sampler::getSampleLabel(&Pipeline::dZdtAtConstantDensity) == "dZdtAtConstantDensity");
        CHECK(Sampler::getSampleLabel(&Pipeline::velocity) == "velocity");
        CHECK(Sampler::getSampleLabel(&Pipeline::frictionFactor) == "frictionFactor");
        CHECK(Sampler::getSampleLabel(&Pipeline::reynoldsNumber) == "reynoldsNumber");
        CHECK(Sampler::getSampleLabel(&Pipeline::ambientTemperature) == "ambientTemperature");
        CHECK(Sampler::getSampleLabel(&Pipeline::heatFlow) == "heatFlow");

        // inlet/outlet comp
        CHECK(Sampler::getSampleLabel(&Pipeline::inletComposition) == "inletComposition");
        CHECK(Sampler::getSampleLabel(&Pipeline::outletComposition) == "outletComposition");
    }

    SUBCASE("make output dir")
    {
        const std::string dirString = std::string(TRANSFLOW_RESOURCE_PATH) + "/test/sampler/existing_folder";
        CHECK_NOTHROW(Sampler::makeOutputDir(dirString));

        auto path = std::filesystem::path(std::filesystem::absolute(std::string(TRANSFLOW_RESOURCE_PATH) + "/test/sampler/existing_folder"));
        CHECK(Sampler::makeOutputDir(dirString) == path);

        CHECK_THROWS(Sampler::makeOutputDir(std::string(TRANSFLOW_RESOURCE_PATH) + "/test/sampler/file_without_extension"));
    }
}
