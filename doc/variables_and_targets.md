# CMake variables and targets {#variables_and_targets}

\tableofcontents

This project provides some CMake variables for use during configuration and other explicitly specified build targets.

\section variables CMake Variables

| Variable           | Description                                                                                    | Possible Values                                 | Default Value |
|--------------------|------------------------------------------------------------------------------------------------|-------------------------------------------------|---------------|
| `CMAKE_BUILD_TYPE` | On a single configuration generator, this string determines the build type of the application. | `Debug`/`Release`/`RelWithDebInfo`/`MinSizeRel` | `Release`     |
| `BUILD_DOCS`       | Used to determine if documentation will be generated.                                          | `ON`/`OFF`                                      | `ON`          |
| `BUILD_TESTS`      | Used to determine if the test executable should be built.                                      | `ON`/`OFF`                                      | `ON`          |
| `BUILD_EXAMPLES`   | Used to determine of the examples should be built.                                             | `ON`/`OFF`                                      | `ON`          |

\section targets Build Targets

| Target           | Description                                                    |
|------------------|----------------------------------------------------------------|
| `[Nothing]`      | Build the application.                                         |
| `test`           | If tests were built, then run all tests.                       |
| `doc`            | If `GEN_DOCS=ON`, then generate documentation using `Doxygen`. |
