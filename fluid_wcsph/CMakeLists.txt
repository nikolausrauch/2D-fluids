set(SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/fluid_wcsph.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/wcsph.cpp"
    )
    
    
set(HDR
        "${CMAKE_CURRENT_SOURCE_DIR}/wcsph.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/kernel.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/scene.h"
    )

add_executable( fluid_wcsph ${SRC} ${HDR} )
target_link_libraries( fluid_wcsph PRIVATE viewer_static utility )

target_compile_features( fluid_wcsph PUBLIC cxx_std_17 )
set_target_properties( fluid_wcsph PROPERTIES CXX_EXTENSIONS OFF )

if(OpenMP_CXX_FOUND AND COMPILE_WITH_OPENMP)
    target_link_libraries(fluid_wcsph PUBLIC OpenMP::OpenMP_CXX)
endif()

