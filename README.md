# TransFlow

A transient pipeline flow simulation library.

View the [html documentation online](https://fsund.github.io/transient-pipeline-flow/) or check `docs/` for documentation.

# Windows/MSYS2 installation
 * Install [MSYS2](https://www.msys2.org/) (only supported by Vista and more recent Windows versions)

## Installing required packages
 * Open MSYS2 MSYS terminal

```
pacman -Syyu
```

 * Kill the window after the process has completed, and open a new MSYS2 MSYS terminal
 * Install the required packages:

```
pacman -Syu
pacman -S base-devel git
pacman -S mingw-w64-x86_64-toolchain
pacman -S mingw-w64-x86_64-gdb  mingw-w64-x86_64-gcc  mingw-w64-x86_64-cmake
pacman -S mingw-w64-x86_64-armadillo mingw-w64-x86_64-hdf5 mingw-w64-x86_64-openblas mingw-w64-x86_64-lapack mingw-w64-x86_64-arpack mingw-w64-x86_64-doxygen 
```

 * If you want to use Qt Creator to edit the code, you can install this here. Be aware that this takes up around 7 Gb.

```
pacman -S mingw-w64-x86_64-qt-creator
```

### SuperLU (required by Armadillo, for sparse matrix support)
 * Open `MSYS2 MinGW 64-bit` terminal
 * Navigate to somewhere we can download and compile SuperLU
 * Download, compile and install SuperLU

```
pacman -S unzip
wget https://github.com/xiaoyeli/superlu/archive/v5.2.0.zip
unzip v5.2.0.zip
cd superlu-5.2.0
mkdir build
cd build
cmake -G'MSYS Makefiles' -Denable_blaslib=OFF -DCMAKE_INSTALL_PREFIX="C:/msys64/mingw64/" ..
make install
```

## Compiling TransFlow
 * Open MSYS2 MinGW 64-bit terminal
 * Set up git for BitBucket if you haven't done so already: https://confluence.atlassian.com/bitbucket/set-up-an-ssh-key-728138079.html
 * Navigave to where you want to keep transflow
 * Get the most recent version from git

```
pacman -S unzip
git clone git@bitbucket.org:polytecresearch/transflow.git
cd transflow
mkdir build
cd build
```

 * Configure

```
cmake -G"MSYS Makefiles" ..
```

 * Check that output from the above command looks okay
 * Compile

```
cmake --build .
```

 * Run any of the compiled executables (in `build/examples/` or `build/test/`)