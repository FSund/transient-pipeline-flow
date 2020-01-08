#include "utilities/numerics.hpp"

#include <iomanip>

using std::cout;
using std::endl;

using arma::vec;
using arma::zeros;
using arma::uword;
using arma::sword;

// Method for solving tridiagonal linear system of equations
// Source: p. 51 in  "Numerical Recipes"
// Solves for a vector x[1..n] the tridiagonal linear set given by equation (2.4.1).
/*
    Store tridiagonal matrix as 3 vectors instead of matrix

    | b1  c1   0  ...                 |   |  x1  |   |  r1  |
    | a2  b2  c2  ...                 |   |  x2  |   |  r2  |
    |             ...                 | * | ...  | = | ···  |
    |             ...  aN−1 bN−1 cN−1 |   | xN−1 |   | rN−1 |
    |             ...   0    aN   bN  |   |  xN  |   |  rN  |

    Note that a(0) and c(n) are undefined and are not referenced by the routine
*/
arma::vec utils::tridag(const vec &a, const vec &b, const vec &c, const vec &r, const uword n)
{
    // TODO: add some assertions about length of vectors here?
    // remove n argument, and use f.ex a.n_elem instead?

    if (b(0) == 0.0) cout << "Error 1 in tridag" << endl;
    // If this happens then you should rewrite your equations as a set of order N − 1, with u2 trivially eliminated.

    double bet;
    vec gam = zeros<vec>(n);
    vec u = zeros<vec>(n);

    bet = b(0);
    u(0) = r(0)/bet;

    for (uword j = (2-1); j <= (n-1); j++) // Decomposition and forward substitution.
    {
        gam(j) = c(j-1)/bet;
        bet = b(j) - a(j)*gam(j);

        if (bet == 0.0)	cout << "Error 2 in tridag" << endl; // Algorithm fails; see below.

        u(j) = (r(j) - a(j)*u(j-1))/bet;
    }

    for (sword j = (sword(n) - 2); j >= 0; j--) // Backsubstitution.
    {
        u(uword(j)) -= gam(uword(j) + 1)*u(uword(j) + 1);
    }

    return u;
}
