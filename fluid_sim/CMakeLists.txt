set(SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/simulation.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/wcsph.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/pcisph.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/iisph.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/dfsph.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/pbf.cpp"

        "${CMAKE_CURRENT_SOURCE_DIR}/utility/perfmonitor.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/utility/particlehelper.cpp"
    )
    
    
set(HDR
        "${CMAKE_CURRENT_SOURCE_DIR}/simulation.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/wcsph.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/pcisph.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/iisph.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/dfsph.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/pbf.h"

        "${CMAKE_CURRENT_SOURCE_DIR}/sphkernel.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/utility/perfmonitor.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/utility/hashgrid.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/utility/particlehelper.h"
    )

add_library( fluid_sim ${SRC} ${HDR} )
target_link_libraries( fluid_sim PUBLIC glm_static )
target_compile_features( fluid_sim PUBLIC cxx_std_20 )
set_target_properties( fluid_sim PROPERTIES CXX_EXTENSIONS OFF )

target_include_directories( fluid_sim PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    )

if(OpenMP_CXX_FOUND AND COMPILE_WITH_OPENMP)
    target_link_libraries(fluid_sim PUBLIC OpenMP::OpenMP_CXX)
endif()
