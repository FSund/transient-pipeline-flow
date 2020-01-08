find_package(Doxygen)

if (DOXYGEN_FOUND)
    set(DOXYGEN_IN Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    option(DOXYGEN_QUIET "Quiet Doxygen" ON)
    option(DOXYGEN_WARNINGS "Show Doxygen warnings" OFF)
    
    set(DOXYGEN_QUIET_VAR NO)
    if (DOXYGEN_QUIET)
        set(DOXYGEN_QUIET_VAR YES)
    endif (DOXYGEN_QUIET)

    set(DOXYGEN_WARNINGS_VAR NO)
    if (DOXYGEN_WARNINGS)
        set(DOXYGEN_WARNINGS_VAR YES)
    endif (DOXYGEN_WARNINGS)

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    
    add_custom_target(doc
        ALL  # build by default/as part of ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating documentation with Doxygen"
        VERBATIM
    )

    message("DOXYGEN_EXECUTABLE: ${DOXYGEN_EXECUTABLE}")
    message("DOXYGEN_OUT: ${DOXYGEN_OUT}")
else (DOXYGEN_FOUND)
    message("Doxygen needs to be installed to generate documentation!")
endif (DOXYGEN_FOUND)
