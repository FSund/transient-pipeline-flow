#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <filesystem>
#include <armadillo>

class Pipeline;

/*!
 * \ingroup typedefs
 * Typedef for pointer to const member function (getter) in Pipeline.
 *
 * Makes it trivial to declare a member-function pointer
 *
 *     PropertyGetter p = &Pipeline::flow;
 *
 * and declaring functions that receive member-function pointers
 *
 *     void function(PropertyGetter p) { ... }
 *
 * and declaring functions that return member-function pointers
 *
 *     PropertyGetter function() { ... }
 *
 * From here: https://isocpp.org/wiki/faq/pointers-to-members#typedef-for-ptr-to-memfn
 */
typedef const arma::vec& (Pipeline::*PropertyGetter)() const;

namespace details
{

/*!
 * \brief The PropertyToSample struct contains information about a single
 * property from Pipeline that we are going to sample.
 *
 * PropertyGetter is a typedef for a member function pointer of Pipeline,
 * and the label is an identifier that is also used to determine the file name.
 */
struct PropertyToSample
{
    /*!
     * \brief Construct from label and a Pipeline member function pointer.
     * \param label Label
     * \param function Pipeline fember function (getter) pointer
     */
    PropertyToSample(const std::string& label, PropertyGetter function):
        label(label),
        function(function)
    {}

    std::string label; //!< Label
    PropertyGetter function; //!< Pipeline member function (getter) pointer
};

}; // end namespace details

/*!
 * \brief The Sampler class is used to sample selected Pipeline properties
 * during simulations.
 *
 * Adding new properties to sample is done via Sampler::addPropertyToPrint.
 *
 * Control of which properties to sample is implemented via pointers to member
 * functions/getters of Pipeline, which are stored in m_samplers, together
 * with a label.
 */
class Sampler
{
public:    
    /*!
     * \brief Construct given output directory and (optional) print interval.
     * \param path Where to store results
     * \param interval How often to sample (max) [s]
     * \param append If we should append to (true) or overwrite existing output files (false)
     */
    Sampler(
            const std::filesystem::path& path,
            const arma::uword interval = 60,
            const bool append = false,
            const arma::uvec indicesToSample = {});

    /*!
     * \brief Add property to save to file. This automatically determines the
     * label and filename of the property.
     * \param samplingFunction Pointer to Pipeline property getter
     */
    void addPropertyToPrint(PropertyGetter samplingFunction);

    /*!
     * \brief Add property to print with explicit label.
     * \param samplingFunction Pointer to Pipeline property getter
     * \param label Property label (filename will be \<label\>.csv)
     */
    void addPropertyToPrint(PropertyGetter samplingFunction, const std::string& label);

    /*!
     * \brief Sample current state. This is usually only called by Simulator,
     * but can also be called to force a sample.
     * \param pipeline Pipeline instance to sample
     * \param forceSample [bool]
     * \return true if a sample was recorded, else false
     */
    bool sample(const Pipeline& pipeline, const bool forceSample = false);

    //! Get (const ref) output directory.
    const std::filesystem::path& outputDir() const { return m_outputDir; }

    //! Automatically determines property label from Pipeline getter function
    //! pointer.
    static std::string getSampleLabel(PropertyGetter samplingFunction);

    //! Make output directory from string. Does some checks to see if directory
    //! exists, and creates it if it doesn't exists. Throws error if directory
    //! is a file, or if directory can't be created.
    static std::filesystem::path makeOutputDir(const std::filesystem::path& path);

    //! Set which indices to sample. Empty vector means that we sample all
    //! points.
    Sampler& setIndicesToSample(const arma::uvec& indices);

private:
    //! Output file streams, one for each property. Should be in the same order
    //! as m_samplers.
    std::vector<std::ofstream> m_outputFiles;
    //! details::PropertyToSample instances with the label and sampling function for
    //! each property
    std::vector<details::PropertyToSample> m_samplers;
    //! Output directory
    const std::filesystem::path m_outputDir;

    //! How often to print (max) [s]
    arma::uword m_printInterval = 60;
    //! Time of last print [s]
    arma::uword m_timeOfLastPrint = 0;
    //! Which grid points to sample. Empty vector means that we sample all
    //! points.
    arma::uvec m_indicesToSample {};

    //! If we append to (true) or overwrite (false) existing files
    bool m_append;
};
