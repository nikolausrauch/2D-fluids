set(SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/fluid_pcisph.cpp"
    )
    
    
set(HDR
        ""
    )

add_executable( fluid_pcisph ${SRC} ${HDR} )
target_link_libraries( fluid_pcisph PRIVATE viewer_static fluid_sim utility )

target_compile_features( fluid_pcisph PUBLIC cxx_std_17 )
set_target_properties( fluid_pcisph PROPERTIES CXX_EXTENSIONS OFF )
