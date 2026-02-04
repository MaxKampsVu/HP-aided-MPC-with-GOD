set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/benchmarks)

set(benchbin "")

macro(add_benchmark source_name)
    add_executable(${source_name} ${source_name}.cpp utils.cpp)
    target_link_libraries(${source_name} Boost::system Boost::program_options nlohmann_json::nlohmann_json Threads::Threads NTL GMP EMPTool protocols)
    list(APPEND benchbin ${source_name})
endmacro()

# benchmarks 
add_benchmark(wasat_mpc)
add_benchmark(wasat_offline)
add_benchmark(wasat_online)
add_benchmark(alhena_mpc)
add_benchmark(alhena_offline)
add_benchmark(alhena_online)
add_benchmark(wasat_ope)
add_benchmark(alhena_ope)


add_custom_target(benchmarks)
add_dependencies(benchmarks ${benchbin} create_ot_data)
