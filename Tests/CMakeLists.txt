#===================================
#
# Tests Project
#
#===================================
add_executable(Tests
    Tests_main.h Tests_main.cpp
    test_sand_hypoplasticity.cpp
    test_input_from_MPM.cpp
    )

target_include_directories(
    Tests PUBLIC
    ${MATERIAL_MODELS_INCLUDE_DIR}
    )

target_link_libraries(
    Tests PUBLIC
    MaterialModels
    )