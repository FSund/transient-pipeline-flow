add_custom_target(examples
    ALL  # build by default/as part of ALL
    )

# from here: https://stackoverflow.com/a/14307434/1850917
file(GLOB APP_SOURCES "*.cpp")
foreach(testsourcefile ${APP_SOURCES})
    get_filename_component(testname ${testsourcefile} NAME_WE)
    add_executable(${testname} ${testsourcefile})
    target_link_libraries(${testname} transflow)
    add_dependencies(examples ${testname})
endforeach(testsourcefile ${APP_SOURCES})
