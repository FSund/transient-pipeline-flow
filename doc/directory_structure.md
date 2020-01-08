# Directory structure {#directory_structure}

Ignoring all extraneous files, this is the structure that this project uses: 

    project-name/
    ├── CMakeLists.txt
    ├── cmake
    ├── doc/
    │   ├── CMakeLists.txt
    │   ├── Doxyfile.in
    │   ├── directory_structure.md
    │   ├── main_page.md
    │   ├── manual.md
    │   └── ...
    ├── examples/
    │   ├── CMakeLists.txt
    │   ├── simple_simulation.cpp
    │   └── ...
    ├── res/
    │   ├── equationofstate
    │   ├── examples
    │   └── test
    ├── src/
    │   ├── config.hpp.in
    │   └── ...
    ├── test/
    │   ├── CMakeLists.txt
    │   ├── debug.hpp
    │   ├── debug.cpp
    │   ├── test_runner.cpp
    │   └── ...
    ├── third_party/
    │   └── doctest/
    │       └── doctest.h
    └── tools/
        └── profile_slider.py

This may all look very complex, but here is an explanation of all of the directories:

| Directory     | Purpose                                                                                                                        |
|---------------|--------------------------------------------------------------------------------------------------------------------------------|
| `cmake`       | Contains all CMake related configuration files.                                                                                |
| `doc`         | Contains Doxygen configuration files, which can be used to create documentation for the project using the CMake target `doc`.  |
| `examples`    | Contains examples that show how to use the application.                                                                        |
| `res`         | Contains different resources required to run the application.                                                                  |
| `src`         | Contains the source code for the application. Private level includes and implementations are also in this directory.           |
| `test`        | Contains test files to unit test the application.                                                                              |
| `third_party` | Contains CMake configuration files for third party dependencies or single header files for projects with single includes.      |
| `tools`       | Contains different useful tools, for example for analysis of results.                                                          |
