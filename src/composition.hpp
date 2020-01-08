#pragma once

#include <ostream>
#include <limits>
#include <armadillo>

/*!
 * \brief The Composition class is a simple container class for the composition
 * of natural gas. We use arma::vec::fixed<10> to fix the number of components
 * at 10.
 *
 * The composition is stored in order C1, C2, C3, iC4, nC4, iC5, nC5, C6, N2, CO2.
 *
 * The only member variable is the composition of the gas. All other functions
 * are for convenience, getters, normalization etc.
 */
class Composition
{
public:
    /*!
     * \brief Empty constructor. Returns composition with all fractions set to zero.
     */
    Composition();

    /*!
     * \brief Construct for arma::vec. Explicit to disallow implicit conversion from arma::vec.
     * \param composition
     */
    explicit Composition(const arma::vec& composition);

    //! Get fraction of component i.
    inline double operator()(const arma::uword i) const { return m_composition(i); }

    //! Get/set fraction of component i.
    inline double& operator()(const arma::uword i) { return m_composition(i); }

    // user-defined conversion (implicit conversion)
    // this allows Composition to be a drop-in replacement for the arma::vec
    // in most cases, which is both positive and negative
    /*!
     * \brief User-defined implicit conversion from Composition to arma::vec.
     *
     * This allows Composition to be a drop-in replacement for the arma::vec
     * in most cases.
     */
    inline operator arma::vec() const { return m_composition; }

    //! Get (const ref) member arma::vec
    inline const arma::vec& vec() const { return m_composition; }

    //! Get (copy of) member arma::vec
    inline arma::vec vec() { return m_composition; }

    /*!
     * \brief For pretty printing.
     *
     * So we can do
     *
     *     std::cout << composition;
     *
     * and get nice printing.
     */
    friend std::ostream& operator << (std::ostream &out, const Composition& c);

    //! Overloaded comparison operator
    friend bool operator==(const Composition& lhs, const Composition& rhs);

    //! Returns true if composition is normalized.
    bool isNormalized() const;

    //! Normalizes the composition and returns refeference.
    Composition& normalize();

    //! Returns normalized copy (self is const.).
    Composition normalized() const;

    //! Number of elements (fixed)
    const static arma::uword n_elem = 10;

    //! Default composition.
    static const Composition defaultComposition; // initialized in .cpp

    // getters
    inline double C1()  const { return m_composition[0]; } //!< Get C1 fraction
    inline double C2()  const { return m_composition[1]; } //!< Get C2 fraction
    inline double C3()  const { return m_composition[2]; } //!< Get C3 fraction
    inline double iC4() const { return m_composition[3]; } //!< Get iC4 fraction
    inline double nC4() const { return m_composition[4]; } //!< Get nC4 fraction
    inline double iC5() const { return m_composition[5]; } //!< Get iC5 fraction
    inline double nC5() const { return m_composition[6]; } //!< Get nC5 fraction
    inline double C6()  const { return m_composition[7]; } //!< Get C6+ fraction
    inline double N2()  const { return m_composition[8]; } //!< Get N2 fraction
    inline double CO2() const { return m_composition[9]; } //!< Get CO2 fraction

    // setters
    inline double& C1()  { return m_composition[0]; } //!< Get/set C1 fraction
    inline double& C2()  { return m_composition[1]; } //!< Get/set C2 fraction
    inline double& C3()  { return m_composition[2]; } //!< Get/set C3 fraction
    inline double& iC4() { return m_composition[3]; } //!< Get/set iC4 fraction
    inline double& nC4() { return m_composition[4]; } //!< Get/set nC4 fraction
    inline double& iC5() { return m_composition[5]; } //!< Get/set iC5 fraction
    inline double& nC5() { return m_composition[6]; } //!< Get/set nC5 fraction
    inline double& C6()  { return m_composition[7]; } //!< Get/set C6+ fraction
    inline double& N2()  { return m_composition[8]; } //!< Get/set N2 fraction
    inline double& CO2() { return m_composition[9]; } //!< Get/set CO2 fraction

private:
    /*!
     * \brief Store composition as fixed size vector, so the size can never be
     * changed. Will throw error if we try to construct Composition with
     * incorrect number of elements.
     *
     * Fractions are stored in order C1, C2, C3, iC4, nC4, iC5, nC5, C6, N2, CO2.
     */
    arma::vec::fixed<10> m_composition;
};
