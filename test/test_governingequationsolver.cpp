#include "debug.hpp"
#include "solver/governingequationsolver.hpp"
#include "solver/discretizer/enthalpy.hpp"
#include "solver/discretizer/internalenergy.hpp"
#include "solver/boundaryconditions.hpp"
#include "physics.hpp"
#include "pipeline.hpp"
#include <iomanip>

#include <memory>
#include <tuple>

using std::unique_ptr;
using std::vector;
using std::string;
using std::endl;
using std::cout;
using arma::uword;
using arma::vec;
using arma::zeros;
using arma::endr;
using arma::mat;

TEST_SUITE_BEGIN("GoverningEquationSolver");

TEST_CASE("constructors")
{
    SUBCASE("strings")
    {
        CHECK_NOTHROW(GoverningEquationSolver(10));
        CHECK_NOTHROW(GoverningEquationSolver<InternalEnergyDiscretizer>(10));
        CHECK_NOTHROW(GoverningEquationSolver<EnthalpyDiscretizer>(10));
    }

    SUBCASE("Discretizer")
    {
        CHECK_NOTHROW(GoverningEquationSolver(std::make_unique<InternalEnergyDiscretizer>(10)));
        CHECK_NOTHROW(GoverningEquationSolver(std::make_unique<EnthalpyDiscretizer>(10)));
    }
}

TEST_CASE("Uniform flow test")
{
    uword nGridPoints = 10;
    vector<unique_ptr<GoverningEquationSolverBase>> solvers;
    solvers.push_back(std::make_unique<GoverningEquationSolver<EnthalpyDiscretizer>>(nGridPoints));
    solvers.push_back(std::make_unique<GoverningEquationSolver<InternalEnergyDiscretizer>>(nGridPoints));

    for (const string m : {"inlet", "outlet", "none", "both"})
    {
        for (const string p : {"inlet", "outlet", "both", "none"})
        {
            if (m == "none" && p != "both") continue;
            if (p == "none" && m != "both") continue;
            for (auto& solver : solvers)
            {
                const vec pressure = zeros<vec>(nGridPoints) + 1e6;
                const vec temperature = zeros<vec>(nGridPoints) + 273.15 + 10;
                const vec massFlow = zeros<vec>(nGridPoints) + 100;

                Pipeline state(nGridPoints);
                state.constantComposition() = true;
                state.pressure() = pressure;
                state.temperature() = temperature;
                state.flow() = massFlow;
                state.diameter().fill(1.0);
                state.setLength(1e4);
                state.roughness().fill(0);
                state.height().fill(0);

                // update properties
                const arma::uword dt = 60;
                Physics physics(state, "BWRS", "SteadyState");
                physics.updateDerivedProperties(state);

                const BoundaryConditions bc(state, {m, p, "inlet"}); // takes inlet and outlet from state

                const mat output = solver->solve(dt, state, state, bc);

                vec tolerances;
                if (solver->isOverDetermined(bc))
                {
                    tolerances = {1e-7, 3e-9, 1e-8};
                }
                else
                {
                    tolerances = {1e-12, 1e-14, 3e-14};
                }
                for (uword i = 0; i < 3; i++) // loop over variables
                {
                    vec relDiff = output.col(i);
                    relDiff = relDiff - bc.inlet(i).value();
                    relDiff = arma::abs(relDiff/bc.inlet(i));
                    CHECK(arma::all(relDiff < tolerances(i)) == true);

                    // print smore more info if test failed
                    if (arma::all(relDiff < tolerances(i)) != true)
                    {
                        const auto propertyOutsideTolerance = [i]()
                        {
                            if (i == 0) return "flow";
                            else if (i == 1) return "pressure";
                            else return "temperature";
                        }();
//                        cout << "property = " << propertyOutsideTolerance << ", \ndiff = " << (diff/bc(i)).t() << ", bc = [" << m << ", " << p << "], n = " << solver.getBoundarySettings().nActiveBoundaryConditions() << ", tolerance = " << tolerances(i) << endl;
                        cout << "property = " << propertyOutsideTolerance << ", max diff = " << arma::max(arma::abs(relDiff)) << ", overdetermined: " << solver->isOverDetermined(bc) << endl;
                    }
                }
            }
        }
    }
}

// same as the above, only more grid points and longer time step
TEST_CASE("Uniform flow test 2")
{
    uword nGridPoints = 100;
    vector<unique_ptr<GoverningEquationSolverBase>> solvers;
    solvers.push_back(std::make_unique<GoverningEquationSolver<EnthalpyDiscretizer>>(nGridPoints));
    solvers.push_back(std::make_unique<GoverningEquationSolver<InternalEnergyDiscretizer>>(nGridPoints));

    for (const string m : {"inlet", "outlet", "none", "both"})
    {
        for (const string p : {"inlet", "outlet", "both", "none"})
        {
            if (m == "none" && p != "both") continue;
            if (p == "none" && m != "both") continue;
            for (auto& solver : solvers)
            {
                const vec pressure = zeros<vec>(nGridPoints) + 1e6;
                const vec temperature = zeros<vec>(nGridPoints) + 273.15 + 10;
                const vec massFlow = zeros<vec>(nGridPoints) + 100;

                // no source/sink terms (friction etc.)
                Pipeline currentState(nGridPoints);
                currentState.pressure() = pressure;
                currentState.temperature() = temperature;
                currentState.flow() = massFlow;
                currentState.diameter().fill(1.0);
                currentState.setLength(1e4);

                // update derived properties
                Physics physics(currentState);
                physics.updateDerivedProperties(currentState);
//                physics.initializeBatchTracking(currentState);


                Pipeline state(currentState);
                const BoundaryConditions bc(state, {m, p, "inlet"}); // takes inlet and outlet from state
                const arma::uword dt = 3000;
                const mat output = solver->solve(dt, state, state, bc);

                vec tolerances;
                if (solver->isOverDetermined(bc))
                {
                    tolerances = {1e-7, 3e-9, 1e-7};
                }
                else
                {
                    tolerances = {3e-11, 2e-13, 6e-13};
                }
                for (uword i = 0; i < 3; i++) // loop over variables
                {
                    vec relDiff = output.col(i);
                    relDiff = relDiff - bc.inlet(i).value();
                    relDiff = arma::abs(relDiff/bc.inlet(i));
                    CHECK(arma::all(relDiff < tolerances(i)) == true);

                    // print some more info if test failed
                    if (arma::all(relDiff < tolerances(i)) != true)
                    {
                        const auto propertyOutsideTolerance = [i]()
                        {
                            if (i == 0) return "flow";
                            else if (i == 1) return "pressure";
                            else return "temperature";
                        }();
//                        cout << "property = " << propertyOutsideTolerance << ", \ndiff = " << (diff/bc(i)).t() << ", bc = [" << m << ", " << p << "], n = " << solver.getBoundarySettings().nActiveBoundaryConditions() << ", tolerance = " << tolerances(i) << endl;
                        cout << "property = " << propertyOutsideTolerance << ", max diff = " << arma::max(arma::abs(relDiff)) << ", overdetermined: " << solver->isOverDetermined(bc) << endl;
                    }
                }
            }
        }
    }
}

TEST_SUITE_END();
