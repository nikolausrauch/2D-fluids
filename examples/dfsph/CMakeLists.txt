set(SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/fluid_dfsph.cpp"
    )
    
    
set(HDR
        ""
    )

add_executable( fluid_dfsph ${SRC} ${HDR} )
target_link_libraries( fluid_dfsph PRIVATE viewer_static fluid_sim utility )

target_compile_features( fluid_dfsph PUBLIC cxx_std_17 )
set_target_properties( fluid_dfsph PROPERTIES CXX_EXTENSIONS OFF )
