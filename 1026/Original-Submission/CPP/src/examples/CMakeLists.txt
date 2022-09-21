add_executable(driver driver.cpp)
add_compilation_flags(driver "${opt_flags}")
target_link_libraries(driver PUBLIC ${CP_CALS})
