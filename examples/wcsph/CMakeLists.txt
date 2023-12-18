set(SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/fluid_wcsph.cpp"
    )
    
    
set(HDR
        ""
    )

add_executable( fluid_wcsph ${SRC} ${HDR} )
target_link_libraries( fluid_wcsph PRIVATE viewer_static fluid_sim utility )

target_compile_features( fluid_wcsph PUBLIC cxx_std_17 )
set_target_properties( fluid_wcsph PROPERTIES CXX_EXTENSIONS OFF )
