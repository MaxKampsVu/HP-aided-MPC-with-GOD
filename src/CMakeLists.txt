add_library(protocols
    utils/helpers.cpp
    utils/circuit.cpp
    utils/types.cpp
    protocols/rand_gen_pool.cpp
    protocols/ot_provider_ha.cpp
    protocols/offline_phase.cpp
    protocols/online_phase.cpp)

target_include_directories(protocols PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/io ${CMAKE_CURRENT_SOURCE_DIR}/utils ${CMAKE_CURRENT_SOURCE_DIR}/protocols)
target_link_libraries(protocols PUBLIC Boost::system EMPTool NTL GMP)
