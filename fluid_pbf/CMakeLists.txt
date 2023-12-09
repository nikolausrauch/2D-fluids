set(SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/fluid_pbf.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/pbfsystem.cpp"
    )
    
    
set(HDR
        "${CMAKE_CURRENT_SOURCE_DIR}/pbfsystem.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/kernel.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/scene.h"
    )

add_executable( fluid_pbf ${SRC} ${HDR} )
target_link_libraries( fluid_pbf PRIVATE viewer_static utility )

target_compile_features( fluid_pbf PUBLIC cxx_std_17 )
set_target_properties( fluid_pbf PROPERTIES CXX_EXTENSIONS OFF )

if(OpenMP_CXX_FOUND AND COMPILE_WITH_OPENMP)
    target_link_libraries(fluid_pbf PUBLIC OpenMP::OpenMP_CXX)
endif()
