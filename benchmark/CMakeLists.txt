set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/benchmarks)

set(benchbin "")

macro(add_benchmark source_name)
    add_executable(${source_name} ${source_name}.cpp utils.cpp)
    target_link_libraries(${source_name} Boost::system Boost::program_options nlohmann_json::nlohmann_json Threads::Threads NTL GMP EMPTool AsyncAsterisk)
    list(APPEND benchbin ${source_name})
endmacro()

# benchmarks 
add_benchmark(hm_async_asterisk_offline)
add_benchmark(hm_async_asterisk_online)
add_benchmark(hm_async_asterisk_mpc)
add_benchmark(dm_async_asterisk_offline)
add_benchmark(dm_async_asterisk_online)
add_benchmark(dm_async_asterisk_mpc)
add_benchmark(dm_async_asterisk_god_mpc)
add_benchmark(dm_async_asterisk_god_offline)
add_benchmark(dm_async_asterisk_god_online)
add_benchmark(dm_sync_asterisk_god_mpc)
add_benchmark(dm_sync_asterisk_god_offline)
add_benchmark(dm_sync_asterisk_god_online)

add_benchmark(dm_async_asterisk_god_ope)
add_benchmark(dm_sync_asterisk_god_ope)

add_benchmark(dm_async_asterisk_ope)


add_custom_target(benchmarks)
add_dependencies(benchmarks ${benchbin} create_ot_data)
