set(SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/examplescenes.cpp"
    )
    
    
set(HDR
        "${CMAKE_CURRENT_SOURCE_DIR}/examplescenes.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/rendering.h"
    )

add_library( utility ${SRC} ${HDR} )
target_link_libraries( utility PUBLIC glm_static viewer_static fluid_sim )
target_compile_features( utility PUBLIC cxx_std_17 )
set_target_properties( utility PROPERTIES CXX_EXTENSIONS OFF )

target_include_directories( utility PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    )
