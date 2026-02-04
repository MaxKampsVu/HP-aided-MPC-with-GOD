set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/tests)

set(testbin "")

macro(add_testfile source_name)
    add_executable(${source_name}_test ${source_name}.cpp)
    target_link_libraries(${source_name}_test PUBLIC Boost::unit_test_framework Threads::Threads NTL GMP EMPTool protocols)
    # add_test(NAME ${source_name}_test COMMAND ${source_name}_test)
    list(APPEND testbin ${source_name}_test)
endmacro()

add_testfile(io)
add_testfile(rand)
add_testfile(utils)
add_testfile(mpc)

add_custom_target(tests)
add_dependencies(tests ${testbin} create_ot_data)
