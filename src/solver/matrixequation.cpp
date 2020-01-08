#include "solver/matrixequation.hpp"

#include "utilities/errors.hpp"
#include "solver/boundaryconditions.hpp"

using arma::mat;
using arma::cube;
using arma::uword;
using arma::vec;
using arma::zeros;
using arma::sp_mat;
using arma::umat;
using std::endl;
using std::cout;

MatrixEquation::~MatrixEquation()
{}

mat MatrixEquation::solve(
        const uword nGridPoints,
        const uword nEquationsAndVariables,
        const BoundaryConditions& boundaryConditions) const
{
    // solves matrix equation Ax = b, and rearranges the output vec into a mat which includes boundaryconditions
    vec x;
    try
    {
        x = this->solveMatrixEquation();
    }
    catch (const std::runtime_error&) // spsolve only throws this, and only when it can't find a solution
    {
        if (m_coefficients.n_rows == m_coefficients.n_cols) // square matrix, can use sparse solver
        {
            arma::superlu_opts settings;
            settings.pivot_thresh = 0.5;
            settings.equilibrate = true;
            cout << "MatrixEquation::solve() could not find a solution. Trying a final trick before giving up";
            try
            {
                x = arma::spsolve(m_coefficients, m_constants, "superlu", settings); // sparse solver, uses superlu by default
            }
            catch (const std::runtime_error&) // spsolve only throws this, and only when it can't find a solution
            {
                cout << " - failed." << endl;
                throw utils::no_solution_found("MatrixEquation::solve(): arma::spsolve() could not find a solution.");
            }
            cout << " - success." << endl;
        }
        else // non-square matrix
        {
            throw utils::no_solution_found("MatrixEquation::solve(): could not find a solution.");
        }
    }
    catch (...) // in case dense solver solve() throws something we're not aware of
    {
        cout << "MatrixEquation::solve(): Unknown error when solving matrix equation, re-throwing exception." << endl;
        std::exception_ptr eptr = std::current_exception(); // capture;
        std::rethrow_exception(eptr);
    }

    return reshapeSolverOutput(x, boundaryConditions, nGridPoints, nEquationsAndVariables);
}

void MatrixEquation::fillCoefficientMatrixAndConstantsVector(
        const uword nGridPoints,
        const uword nEquationsAndVariables,
        const BoundaryConditions& boundaryConditions,
        const cube& term_i,
        const cube& term_ipp,
        const mat& boundaryTerms)
{
    // fills in the coefficient matrix A (m_coefficients) and the constants vector b (m_constants)
    const uword N = nEquationsAndVariables; // per grid point
    if (term_i.n_slices != N || term_ipp.n_slices != N || boundaryTerms.n_cols != N)
    {
        cout << "ERROR: MatrixEquation::fillMatrixAndVector(): something is wrong, wrong number of elements!" << endl;
        exit(EXIT_FAILURE);
    }

    // (i,j,k) == (row, col, slice) == (grid point, equation #, variable)
    // For cube and 3D field classes, access the element/object stored at the i-th row, j-th column and k-th slice.

    // in term_i and term_ipp
    // dimensions (nGridPoints - 1, nEquationsAndVariables, nEquationsAndVariables)
    // each slice contains the factors for one variable (x, y, z, ...)
    // each column contains the factors for one equation
    // each row contains the factors for one grid point (y_0, y_1, y_2, ..., y_N)

    // if we have more than 3 boundary conditions (for example if we have reverse flow at outlet we need outlet temperature as bc, so will give both inlet and outlet temperature)
    // we remove one variable from the system -> remove one row from unknowns x and one column from coefficients matrix A
    const uword nElements = nGridPoints - 1; // number of finite difference "elements"
    if (boundaryConditions.nActiveBoundaryConditions() < N)
        throw std::runtime_error("too few boundary conditions");
    const uword nExtraConditions = boundaryConditions.nActiveBoundaryConditions() - N; // number of boundary conditions exceeding N (the minimum amount)
    const uword nRows = N*nElements; // nEquations for each element
    const uword nCols = N*nElements - nExtraConditions; // nVariables for each element, but need to subtract for boundary conditions

//    m_coefficients = sp_mat(nRows, nCols); // don't construct here, use batch insertion later
    m_constants = zeros<vec>(nRows); // the constants/knowns in each equation for each element

    // use batch insertion
    // we need 2*nEquations*nVariables elements for each element I, but subtract nEquations for each boundary condition (removes one column from A)
    //  - the actual insertion is performed in the last line of this function
    uword nElementsToInsert = 2*N*N*nElements - (N + nExtraConditions)*N;

    // locations is a dense matrix of type umat, with a size of 2 x N, where N is the number of values to be inserted;
    // the location of the i-th element is specified by the contents of the i-th column of the locations matrix,
    // where the row is in locations(0,i), and the column is in locations(1,i)
    umat locations = zeros<umat>(2, nElementsToInsert);
    vec values = zeros<vec>(nElementsToInsert);
    uword c = 0; // counter of where we are in locations and values

    uword colOffset = 0; // only used to find starting point for col0

    // fill in column by column
    // TODO: write more documentation

    // do boundary terms first
    // boundaryTerms dimensions (nGridPoints - 1, nEquationsAndVariables)
    for (uword element = 0; element < nElements; element++)
    {
        for (uword eq = 0; eq < N; eq++)
        {
            m_constants(N*element + eq) = boundaryTerms(element, eq);
        }
    }

    // first grid point
    // each column contains coefficients of y_i^(n+1) from only the equations for the first element (I = 1)
    {
        uword element = 0;
        uword col = 0; // column in A matrix
        for (uword var = 0; var < N; var++) // variable number (flow, pressure, temperature)
        {
            if (boundaryConditions.inlet(var).isActive())
            {
                for (uword eq = 0; eq < N; eq++) // equation number (continuity, momentum, energy)
                {
                    // term y_i^(n+1) is KNOWN, so the whole terms goes into the constants/known vector b
                    m_constants(eq) -= term_i(element, eq, var)*boundaryConditions.inlet(var);
                }
                colOffset++; // have filled in one column, so increment this offset
            }
            else
            {
                for (uword eq = 0; eq < N; eq++) // equation number (continuity, momentum, energy)
                {
                    // term y_i^(n+1) is unknown, so the coefficient in front of that goes into the unknowns matrix A
                    uword row = eq;
                    locations(0, c) = row;
                    locations(1, c) = col;
                    values(c) = term_i(element, eq, var);
                    c++;
                }
                col++; // only incremented if not a boundary condition at the inlet
            }
        }
    }

    // central columns, which follow the pattern [m1, p1, T1, m2, p2, T2, ...]
    uword row0 = 0;
    uword col0 = N - colOffset; // start at correct column number

    for (uword element = 0; element < (nElements - 1); element++) // loop over elements I, skip last one since we use term_i(element+1) and term_ipp(lastOne) is handled by the final loop below
    {
        for (uword eq = 0; eq < N; eq++) // equation number (continuity, momentum, energy)
        {
            for (uword var = 0; var < N; var++) // variable number (flow, pressure, temperature)
            {
                uword row;
                uword col = col0 + var; // one column per iteration

                // coeffient of y_(i+1)^(n+1) from equation for element I --> y_1
                row = row0 + eq;
                locations(0, c) = row;
                locations(1, c) = col;
                values(c) = term_ipp(element, eq, var);
                c++;

                // coefficient of y_i^(n+1) from equation for element I+1 --> y_1
                row = row0 + (eq + N); // one block down
                locations(0, c) = row;
                locations(1, c) = col;
                values(c) = term_i(element+1, eq, var); // grid + 1 to get element I+1, term_i to get y_i
                c++;
            }
        }
        row0 += N;
        col0 += N;
    }

    // final grid point
    // each column contains coefficients of y_(i+1)^(n+1) from only the equations for the final element
    {
        uword element = nElements - 1; // final element
        uword col = col0;
        for (uword var = 0; var < N; var++) // variable number (flow, pressure, temperature)
        {
            if (boundaryConditions.outlet(var).isActive())
            {
                for (uword eq = 0; eq < N; eq++) // equation number (continuity, momentum, energy)
                {
                    uword row = row0 + eq;
                    m_constants(row) -= term_ipp(element, eq, var)*boundaryConditions.outlet(var); // y_(i+1)^(n+1) for element I
                }
                colOffset++; // unused here...
            }
            else
            {
                for (uword eq = 0; eq < N; eq++) // equation number (continuity, momentum, energy)
                {
                    // coeffient of y_(i+1)^(n+1) from equation for final element
                    uword row = row0 + eq;
                    locations(0, c) = row;
                    locations(1, c) = col;
                    values(c) = term_ipp(element, eq, var);
                    c++;
                }
                col++; // only incremented if not a boundary condition at the outlet
            }
        }
    }

    // options for sp_mat batch insertion constructor
    bool add_values = false; // when set to true, identical locations are allowed, and the values at identical locations are added
    bool sort_locations = true; // If sort_locations is set to false, the locations matrix is assumed to contain locations that are already sorted according to column-major ordering
    bool check_for_zeros = true; // If check_for_zeros is set to false, the values vector is assumed to contain no zero values
    uword n_rows = nRows;
    uword n_cols = nCols;

//    sp_mat test(add_values, locations, values, n_rows, n_cols, sort_locations, check_for_zeros);
//    cout << arma::accu(arma::abs(Amatrix - test)) << endl;

    // use batch insertion to populate Amatrix
    m_coefficients = sp_mat(add_values, locations, values, n_rows, n_cols, sort_locations, check_for_zeros); // use batch insertion constructor
}

mat MatrixEquation::reshapeSolverOutput(
        const vec& x,
        const BoundaryConditions& boundaryConditions,
        const uword nGridPoints,
        const uword nVariables) const
{
    // transforms output from arma::solve() (vec) to matrix, and adds boundary conditions where they belong
    mat output = zeros<mat>(nGridPoints, nVariables);

    // count number of inlet bc's
    uword nInletBCs = 0;
    for (uword i = 0; i < nVariables; i++)
    {
        if (boundaryConditions.inlet(i).isActive())
        {
            nInletBCs++;
        }
    }

    // do y_0 and inlet boundary conditions first
    uword i0 = 0; // index in x
    uword grid = 0; // grid point
    for (uword var = 0; var < nVariables; var++)
    {
        if (boundaryConditions.inlet(var).isActive())
        {
            output(grid, var) = boundaryConditions.inlet(var);
        }
        else
        {
            output(grid, var) = x(i0);
            i0++;
        }
    }
    // i0 should now have the correct index

    grid = 1;
    while (i0 < nVariables*(nGridPoints - 2)) // stop before outlet
    {
        for (uword var = 0; var < nVariables; var++)
        {
            output(grid, var) = x(i0 + var);
        }
        i0 += nVariables;
        grid++;
    }

    // finally, do y_N and outlet BCs
    grid = nGridPoints - 1;
    // i0 should be correct
    for (uword var = 0; var < nVariables; var++)
    {
        if (boundaryConditions.outlet(var).isActive())
        {
            output(grid, var) = boundaryConditions.outlet(var); // 1 selects outlet
        }
        else
        {
            output(grid, var) = x(i0);
            i0++;
        }
    }

    return output;
}

// private
vec MatrixEquation::solveMatrixEquation() const
{
    // solve matrix equation A*x = b
    vec x;
    if (m_coefficients.n_rows > m_coefficients.n_cols)
    {
        // over-determined, need to use the regular solver
        // WARNING: this is slower and increases memory usage, since we need to convert the sparse matrix to a dense matrix first
        // it's probably possible to avoid this by solving the extra equation(s) for m_i/p_i/T_i and inserting into a equation i+1/i-1, thereby eliminating some rows
        mat A(m_coefficients); // convert sparse matrix to dense matrix
        x = arma::solve(A, m_constants); // dense matrix solver
    }
    else if (m_coefficients.n_rows < m_coefficients.n_cols)
    {
        throw std::runtime_error("under-determined system, this is not implemented (probably caused by a user error)");
    }
    else
    {
        // critically determined system (square matrix)
        x = arma::spsolve(m_coefficients, m_constants, "superlu"); // sparse matrix solver, uses superlu by default
    }

    return x;
}
