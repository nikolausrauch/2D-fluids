set(SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/fluid_iisph.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/iisph.cpp"
    )
    
    
set(HDR
        "${CMAKE_CURRENT_SOURCE_DIR}/iisph.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/kernel.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/scene.h"
    )

add_executable( fluid_iisph ${SRC} ${HDR} )
target_link_libraries( fluid_iisph PRIVATE viewer_static utility )

target_compile_features( fluid_iisph PUBLIC cxx_std_17 )
set_target_properties( fluid_iisph PROPERTIES CXX_EXTENSIONS OFF )

if(OpenMP_CXX_FOUND AND COMPILE_WITH_OPENMP)
    target_link_libraries(fluid_iisph PUBLIC OpenMP::OpenMP_CXX)
endif()

