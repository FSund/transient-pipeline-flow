#pragma once

#include <stdexcept>

namespace utils
{
    class physics_error: public std::range_error
    {
        using std::range_error::range_error; // inherit constructors
    };

    class temperature_range_error: public physics_error
    {
        using physics_error::physics_error; // inherit constructors
    };

    class pressure_range_error: public physics_error
    {
        using physics_error::physics_error; // inherit constructors
    };

    class no_convergence_error: public std::runtime_error
    {
        using runtime_error::runtime_error; // inherit constructors
        unsigned long long m_nIterations;

    public:
        no_convergence_error(const std::string& __arg, const unsigned long long nIterations):
            runtime_error(__arg), // inherited constructor
            m_nIterations(nIterations)
        {}

        unsigned long long nIterations() const { return m_nIterations; }
    };

    class linalg_error: public std::runtime_error
    {
        using runtime_error::runtime_error; // inherit constructors
    };

    class no_solution_found: public linalg_error
    {
        using linalg_error::linalg_error; // inherit constructors
    };
}
