add_library(AsyncAsterisk
    utils/helpers.cpp
    utils/circuit.cpp
    utils/types.cpp
    asyncAsterisk/rand_gen_pool.cpp
    asyncAsterisk/ot_provider.cpp
    asyncAsterisk/ot_provider_ha.cpp
    asyncAsterisk/dm_offline_evaluator.cpp
    asyncAsterisk/dm_online_evaluator.cpp
    asyncAsterisk/hm_offline_evaluator.cpp
    asyncAsterisk/hm_online_evaluator.cpp
    asyncAsterisk/dmgod_offline_evaluator.cpp
    asyncAsterisk/dmgod_online_evaluator.cpp)

target_include_directories(AsyncAsterisk PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/io ${CMAKE_CURRENT_SOURCE_DIR}/utils ${CMAKE_CURRENT_SOURCE_DIR}/asyncAsterisk)
target_link_libraries(AsyncAsterisk PUBLIC Boost::system EMPTool NTL GMP)
