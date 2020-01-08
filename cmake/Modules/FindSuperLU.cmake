find_path(SUPERLU_INCLUDE_DIRS
  NAMES
  supermatrix.h
  PATHS
  $ENV{SUPERLUDIR}
  ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES
  superlu
  SRC
)

find_path(SUPERLU_INCLUDE_DIRS
  NAMES supermatrix.h
  PATH_SUFFIXES "superlu" "SuperLU" "include/superlu" "include" "SRC"
)

find_library(SUPERLU_LIBRARIES superlu PATHS $ENV{SUPERLUDIR} ${LIB_INSTALL_DIR} PATH_SUFFIXES lib)
  
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SUPERLU DEFAULT_MSG
                                  SUPERLU_INCLUDE_DIRS SUPERLU_LIBRARIES)

mark_as_advanced(SUPERLU_INCLUDE_DIRS SUPERLU_LIBRARIES)

if(SUPERLU_FOUND)
  set(SUPERLU_INCLUDE_DIRS ${SUPERLU_INCLUDE_DIR})
  set(SUPERLU_LIBRARIES    ${SUPERLU_LIBRARY})
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determining location of ${SUPERLU_WITH_VERSION} succeeded:\n"
    "Include directory: ${SUPERLU_INCLUDE_DIRS}\n"
    "Library directory: ${SUPERLU_LIBRARIES}\n\n")
endif()