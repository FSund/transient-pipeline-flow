### Building the code {#compilation}

Prerequisites:

To compile the program the following libraries are required
 - [CMake](https://cmake.org/install/) (At least v3.1): For building the code.
 - [Armadillo](http://arma.sourceforge.net/) (tested with 9.600.6)
 - [SuperLU](https://github.com/xiaoyeli/superlu) (version 5.2 required by Armadillo)
 - [HDF5](https://www.hdfgroup.org/downloads/hdf5) (tested with 1.10.5)
 - [LAPACK](http://www.netlib.org/lapack/) (tested with 3.8.0)

Optionally, [doxygen](https://github.com/doxygen/doxygen) (tested with 1.8.16) is required to build the documentation.

Check out `../README.md` for instructions on how to set up a complete build system from the ground up on Windows, using MSYS2.

### Compilation

Compiling the code is as simple as

```
git clone <transflow repo>
cd transflow
mkdir build
cd build
cmake -G"MSYS Makefiles" ..
```

If all libraries are installed correctly, the `cmake` command should complete successfully, and compilation can be performed via

```
cmake --build .
```

### Running the Tests

In order to run the tests, the option `BUILD_TESTS` needs to be set to `ON` during configuration. Then, invoking the tests is as simple as running:

    ctest

from the build folder.

If you want some more output from the tests you can also run the test executable directly

    ./test/test_runner

### Documentation

In order to build the documentation, the option `GEN_DOCS` needs to be set to `ON` during configuration. After compiling with `cmake` the documentation can be built via

    doxygen ./doc/Doxyfile

Then the documentation should be accessible by viewing the file

    ./doc/html/index.html

in a web browser.