set(SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/fluid_pbf.cpp"
    )
    
    
set(HDR
        ""
    )

add_executable( fluid_pbf ${SRC} ${HDR} )
target_link_libraries( fluid_pbf PRIVATE viewer_static fluid_sim utility )

target_compile_features( fluid_pbf PUBLIC cxx_std_17 )
set_target_properties( fluid_pbf PROPERTIES CXX_EXTENSIONS OFF )
