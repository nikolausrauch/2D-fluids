set(SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/fluid_iisph.cpp"
    )
    
    
set(HDR
        ""
    )

add_executable( fluid_iisph ${SRC} ${HDR} )
target_link_libraries( fluid_iisph PRIVATE viewer_static fluid_sim utility )

target_compile_features( fluid_iisph PUBLIC cxx_std_17 )
set_target_properties( fluid_iisph PROPERTIES CXX_EXTENSIONS OFF )
