#===================================
#
# MaterialModels Project
#
#===================================
add_library(MaterialModels STATIC
    LinkList.hpp
    MaterialModels_pcp.h MaterialModels_pcp.cpp
    MaterialModel.h
    MaterialModelInitializer.h MaterialModelInitializer.cpp
    SandHypoplasticityByUmat.h SandHypoplasticityByUmat.cpp
    SandHypoplasticityCavitationByUmat.h SandHypoplasticityCavitationByUmat.cpp
    )

set(MATERIAL_MODELS_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE STRING INTERNAL FORCE)

target_link_libraries(
    MaterialModels PRIVATE
    SandHypoplasticity
    )

if (MSVC)
    target_compile_definitions(MaterialModels PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_options(MaterialModels PRIVATE -YcMaterialModels_pcp.h)
endif()
