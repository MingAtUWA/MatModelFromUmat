#===================================
#
# TestsF90 Project
#
#===================================
add_executable(
    TestsF90
    main.f90
    )

target_link_libraries(
    TestsF90 PUBLIC
    SandHypoplasticity
    )
