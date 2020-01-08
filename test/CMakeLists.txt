set(TEST_SOURCES
    debug.cpp
    test_eos.cpp
    test_eos_bwrs.cpp
    test_eos_gerg04.cpp
    test_eos_idealgas.cpp
    test_heattransfer.cpp
    test_heat_steadystate.cpp
    test_heat_unsteady.cpp
    test_heat_radial.cpp
    test_utils.cpp
    test_constants.cpp
    test_governingequationsolver.cpp
    test_matrix_equation.cpp
    test_discretizer_enthalpy.cpp
    test_discretizer_internalenergy.cpp
    test_batchtracking.cpp
    test_pipeline.cpp
    test_boundaryconditions.cpp
    test_composition.cpp
    test_physics.cpp
    test_simulator.cpp
    test_solver.cpp
    test_materials.cpp
    test_timeseries.cpp
    test_sampler.cpp
    test_phys_utils.cpp
    test_linearinterpolator.cpp
)

add_library(doctest INTERFACE)
target_include_directories(doctest INTERFACE 
    "${CMAKE_SOURCE_DIR}/third_party/doctest"
)

add_executable(test_runner test_runner.cpp ${TEST_SOURCES})
target_link_libraries(test_runner transflow doctest)

# disable warning about unknown pragmas
# since we have #pragma clang diagnostic ignored "-Woverloaded-shift-op-parentheses"
# we get a warning about unknown pragmas, which is turned into an error since we
# have warnings as errors enabled
target_compile_options(test_runner PRIVATE "-Wno-unknown-pragmas")

# don't print errors, since doctest catches the exceptions that are thrown
# alongside the prints, and the prints interfere with doctest's output
target_compile_definitions(test_runner PRIVATE ARMA_DONT_PRINT_ERRORS)

add_test(all_tests test_runner)
