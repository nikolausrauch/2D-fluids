set(SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/fluid_pcisph.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/pcisph.cpp"
    )
    
    
set(HDR
        "${CMAKE_CURRENT_SOURCE_DIR}/pcisph.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/kernel.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/scene.h"
    )

add_executable( fluid_pcisph ${SRC} ${HDR} )
target_link_libraries( fluid_pcisph PRIVATE viewer_static utility )

target_compile_features( fluid_pcisph PUBLIC cxx_std_17 )
set_target_properties( fluid_pcisph PROPERTIES CXX_EXTENSIONS OFF )

if(OpenMP_CXX_FOUND AND COMPILE_WITH_OPENMP)
    target_link_libraries(fluid_pcisph PUBLIC OpenMP::OpenMP_CXX)
endif()

