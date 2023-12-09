set(SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/perfmonitor.cpp"
    )
    
    
set(HDR
        "${CMAKE_CURRENT_SOURCE_DIR}/perfmonitor.h"
    )

add_library( utility ${SRC} ${HDR} )
target_compile_features( utility PUBLIC cxx_std_17 )
set_target_properties( utility PROPERTIES CXX_EXTENSIONS OFF )

target_include_directories( utility PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    )


if(OpenMP_CXX_FOUND AND COMPILE_WITH_OPENMP)
    target_link_libraries(utility PUBLIC OpenMP::OpenMP_CXX)
endif()
