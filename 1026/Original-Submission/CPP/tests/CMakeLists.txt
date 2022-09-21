function(create_test test source)

    set(name ${test}_${BLA_VEN})
    add_executable(${name} ${source}.cpp)
    if (GTEST_DEPENDENCY)
        add_dependencies(${name} ${GTEST_DEPENDENCY})
    endif ()

    target_link_libraries(${name} PRIVATE CALS_${BLA_VEN})
    target_link_libraries(${name} PRIVATE ${GTEST_LINK_FLAGS})

    add_test(NAME ${name} COMMAND ${name})
endfunction(create_test)

set(tests als cals)
foreach (test_dir ${tests})
    add_subdirectory(${test_dir})
endforeach ()
