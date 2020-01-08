#include "sampler.hpp"

#include <functional> // std::invoke
#include <iomanip>
#include "pipeline.hpp"

Sampler::Sampler(
        const std::filesystem::path& path,
        const arma::uword interval,
        const bool append,
        const arma::uvec indicesToSample):
    m_outputDir(makeOutputDir(path)),
    m_printInterval(interval),
    m_indicesToSample(indicesToSample),
    m_append(append)
{
    addPropertyToPrint(&Pipeline::flow);
    addPropertyToPrint(&Pipeline::pressure);
    addPropertyToPrint(&Pipeline::temperature);
}

void Sampler::addPropertyToPrint(PropertyGetter samplingFunction)
{
    auto label = getSampleLabel(samplingFunction);
    addPropertyToPrint(samplingFunction, label);
}

void Sampler::addPropertyToPrint(PropertyGetter samplingFunction, const std::string& label)
{
    std::filesystem::path filePath = m_outputDir / (label + ".csv");

    auto mode = std::ofstream::out;
    if (m_append)
        mode = std::ofstream::app;

    std::ofstream file(filePath, mode);

    // check if we managed to open file
    if (file.is_open())
    {
        file << std::setiosflags(std::ios::right | std::ios::scientific);
        file.precision(8);

        m_outputFiles.push_back(std::move(file));
//        std::cout << "added output file \"" << filePath << "\"" << std::endl;
    }
    else
    {
        throw std::runtime_error("could not create file \"" + filePath.string() + "\"");
    }

    m_samplers.push_back(details::PropertyToSample(label, samplingFunction));
}

bool Sampler::sample(const Pipeline& pipeline, const bool force)
{
    if (pipeline.timestamp() - m_timeOfLastPrint < m_printInterval
            && pipeline.timestamp() > 0
            && !force)
    {
        return false; // don't print if too little time has passed
    }

    for (std::size_t i = 0; i < m_samplers.size(); i++)
    {
        const arma::vec& data = std::invoke(m_samplers.at(i).function, pipeline);

        auto& fout = m_outputFiles.at(i);

        // timestamp
        fout << std::setw(8) << pipeline.timestamp()
             << ",";

        if (m_samplers.at(i).function == &Pipeline::inletComposition
                || m_samplers.at(i).function == &Pipeline::outletComposition
                || m_indicesToSample.n_elem == 0) // all elements
        {
            // print elements manually to avoid arma pretty printing
            // first n-1 points
            for (std::size_t i = 0; i < data.n_elem - 1; i++)
            {
                fout << std::setw(16) << data(i)
                     << ",";
            }
            // final point
            fout << std::setw(16) << data.tail(1)(0);
            fout << std::endl;
        }
        else
        {
            // print elements manually to avoid arma pretty printing
            if (m_indicesToSample.n_elem > 1)
            {
                // first n-1 points
                for (arma::uword i : m_indicesToSample(arma::span(0, m_indicesToSample.n_elem - 2)))
                {
                    fout << std::setw(16) << data(i)
                         << ",";
                }
            }
            // final point
            fout << std::setw(16) << data(m_indicesToSample.back());
            fout << std::endl;
        }
    }

    m_timeOfLastPrint = pipeline.timestamp();

    return true;
}

std::filesystem::path Sampler::makeOutputDir(const std::filesystem::path& path)
{
    // convert to absolute form
    const auto outputDir = std::filesystem::absolute(path);

    // check if dir exists
    if (std::filesystem::exists(outputDir))
    {
        if (std::filesystem::is_directory(outputDir))
        {
//            std::cout << "using existing output directory \"" << outputDir.string() << "\"" << std::endl;
        }
        else
        {
            throw std::runtime_error("wanted output path \"" + outputDir.string() + "\" already exists, but isn't a directory");
        }
    }
    else
    {
        // create dir if it doesn't exist
        if (std::filesystem::create_directories(outputDir))
        {
//            std::cout << "created directory \"" << outputDir << "\"" << std::endl;
        }
        else
        {
            throw std::runtime_error("could not create directory \"" + outputDir.string() + "\"");
        }
    }

    return outputDir;
}

Sampler& Sampler::setIndicesToSample(const arma::uvec& indices)
{
    m_indicesToSample = indices;
    return *this;
}

std::string Sampler::getSampleLabel(PropertyGetter samplingFunction)
{
    if (samplingFunction == PropertyGetter(&Pipeline::flow))
        return "flow";
    else if (samplingFunction == PropertyGetter(&Pipeline::pressure))
        return "pressure";
    else if (samplingFunction == PropertyGetter(&Pipeline::temperature))
        return "temperature";
    else if (samplingFunction == PropertyGetter(&Pipeline::inletComposition))
        return "inletComposition";
    else if (samplingFunction == PropertyGetter(&Pipeline::outletComposition))
        return "outletComposition";
    else if (samplingFunction == PropertyGetter(&Pipeline::heatCapacityConstantVolume))
        return "heatCapacityConstantVolume";
    else if (samplingFunction == PropertyGetter(&Pipeline::heatCapacityConstantPressure))
        return "heatCapacityConstantPressure";
    else if (samplingFunction == PropertyGetter(&Pipeline::density))
        return "density";
    else if (samplingFunction == PropertyGetter(&Pipeline::viscosity))
        return "viscosity";
    else if (samplingFunction == PropertyGetter(&Pipeline::specificGasConstant))
        return "specificGasConstant";
    else if (samplingFunction == PropertyGetter(&Pipeline::molarMass))
        return "molarMass";
    else if (samplingFunction == PropertyGetter(&Pipeline::compressibilityFactor))
        return "compressibilityFactor";
    else if (samplingFunction == PropertyGetter(&Pipeline::dZdtAtConstantPressure))
        return "dZdtAtConstantPressure";
    else if (samplingFunction == PropertyGetter(&Pipeline::dZdpAtConstantTemperature))
        return "dZdpAtConstantTemperature";
    else if (samplingFunction == PropertyGetter(&Pipeline::dZdtAtConstantDensity))
        return "dZdtAtConstantDensity";
    else if (samplingFunction == PropertyGetter(&Pipeline::velocity))
        return "velocity";
    else if (samplingFunction == PropertyGetter(&Pipeline::frictionFactor))
        return "frictionFactor";
    else if (samplingFunction == PropertyGetter(&Pipeline::reynoldsNumber))
        return "reynoldsNumber";
    else if (samplingFunction == PropertyGetter(&Pipeline::ambientTemperature))
        return "ambientTemperature";
    else if (samplingFunction == PropertyGetter(&Pipeline::heatFlow))
        return "heatFlow";
    else
        throw std::runtime_error("unknown sampling function");
}
