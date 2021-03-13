cmake_minimum_required(VERSION 3.15.0)

project(MatModelFromUmat)

enable_language(Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(FortranUtils/)

# umats and their wrapper in Fortran
add_subdirectory(SandHypoplasticity/)
add_subdirectory(SandHypoplasticityCavitation/)

add_subdirectory(MaterialModels/)

add_subdirectory(Tests/)

add_subdirectory(TestsF90/)
