#include "composition.hpp"

#include <limits>
#include <iomanip>

Composition::Composition():
    Composition(arma::zeros<arma::vec>(10))
{}

Composition::Composition(const arma::vec& composition):
    m_composition(composition)
{}

bool Composition::isNormalized() const
{
    // use numeric_limits<double>::epsilon()*m_composition.n_elem
    // since we have epsilon per component, and mat is Mat<double>
    return std::abs((1.0 - arma::sum(m_composition))) <
            std::numeric_limits<double>::epsilon()*m_composition.n_elem;
}

Composition& Composition::normalize()
{
    double sum = arma::sum(m_composition);
    if (sum <= 0)
        throw std::runtime_error("sum of components <= 0");
    m_composition /= sum;
    return *this;
}

Composition Composition::normalized() const
{
     return Composition(m_composition).normalize();
}

std::ostream& operator <<(std::ostream& out, const Composition& c)
{
    // store flags
    std::ios_base::fmtflags flags(out.flags());

    const int width = 10;
    out.setf(std::ios::left);
    out << std::setw(width) << "C1";
    out << std::setw(width) << "C2";
    out << std::setw(width) << "C3";
    out << std::setw(width) << "iC4";
    out << std::setw(width) << "nC4";
    out << std::setw(width) << "iC5";
    out << std::setw(width) << "nC5";
    out << std::setw(width) << "C6+";
    out << std::setw(width) << "N2";
    out << std::setw(width) << "CO2";
    out << std::endl;

//    out.setf(std::ios::fixed);
    for (std::size_t i = 0; i < c.m_composition.size(); i++)
    {
        out << std::setprecision(7) << std::fixed << std::setw(width) << c.m_composition[i];
    }

    // restore flags
    out.flags(flags);
    return out;
}

const Composition Composition::defaultComposition = Composition(arma::vec({89.16, 7.3513, 0.5104, 0.0311, 0.0251, 0.0024, 0.0009, 0.0, 0.6980, 2.2208})/100);

bool operator==(const Composition& lhs, const Composition& rhs)
{
    return arma::all(
            arma::abs(lhs.m_composition - rhs.m_composition)
            < std::numeric_limits<double>::epsilon()
        );
}

