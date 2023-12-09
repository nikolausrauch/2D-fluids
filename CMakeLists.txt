#########################################
#         2D Fluid Simulations          #
#---------------------------------------#
# @author Nikolaus Rauch                #
# @date   04.11.2023                    #
#---------------------------------------#
#########################################
cmake_minimum_required( VERSION 3.9 )
project( 2D-Fluids )

message("\n * 2D Fluid Simulations ")

message(STATUS "${PROJECT_NAME} build environment --")
message(STATUS "Build system:     ${CMAKE_SYSTEM_NAME}-${CMAKE_HOST_SYSTEM_PROCESSOR}")


#########################################
#              CMake-Stuff              #
#########################################
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)


#########################################
#        Global Compiler Options        #
#########################################
set(GCC_COMPILE_OPTIONS "-Wall")
set(GCC_COMPILE_DEBUG_OPTIONS "${GCC_COMPILE_OPTIONS};-ggdb;-O0")
set(GCC_COMPILE_RELEASE_OPTIONS "${GCC_COMPILE_OPTIONS};-O3")

set(MSVC_COMPILE_OPTIONS "/Wall")
set(MSVC_COMPILE_DEBUG_OPTIONS "${MSVC_COMPILE_OPTIONS};/Zi")
set(MSVC_COMPILE_RELEASE_OPTIONS "${MSVC_COMPILE_OPTIONS};/O2")
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

add_compile_options("$<$<AND:$<CXX_COMPILER_ID:GNU>,$<CONFIG:DEBUG>>:${GCC_COMPILE_DEBUG_OPTIONS}>")
add_compile_options("$<$<AND:$<CXX_COMPILER_ID:GNU>,$<CONFIG:RELEASE>>:${GCC_COMPILE_RELEASE_OPTIONS}>")

add_compile_options("$<$<AND:$<CXX_COMPILER_ID:MSVC>,$<CONFIG:DEBUG>>:${MSVC_COMPILE_DEBUG_OPTIONS}>")
add_compile_options("$<$<AND:$<CXX_COMPILER_ID:MSVC>,$<CONFIG:RELEASE>>:${MSVC_COMPILE_RELEASE_OPTIONS}>")


#########################################
#        Build External-Libraries       #
#########################################
add_subdirectory(external/glm)
set_property(TARGET glm_static APPEND_STRING PROPERTY COMPILE_FLAGS " -w")
target_include_directories( glm_static PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external>
    )

add_subdirectory(external/glfw)
set_property(TARGET glfw APPEND_STRING PROPERTY COMPILE_FLAGS " -w")
target_include_directories( glfw PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/glfw/include>
    )

add_subdirectory(external/imgui)
set_property(TARGET imgui_static APPEND_STRING PROPERTY COMPILE_FLAGS " -w")

add_subdirectory(external/stb_image)
set_property(TARGET stb_image APPEND_STRING PROPERTY COMPILE_FLAGS " -w")



#################################
#            OpenMP             #
#################################
option(COMPILE_WITH_OPENMP "enable openmp parallelization" ON)

find_package(OpenMP)
if(OpenMP_CXX_FOUND AND COMPILE_WITH_OPENMP)
    message(STATUS "Build with OpenMP version ${OpenMP_CXX_VERSION}")
else()
    message(STATUS "Build without OpenMP")
endif()


#########################################
#              Build Viewer             #
#########################################
add_subdirectory(viewer)

#########################################
#    Smooth-Particle-Hydrodynamics      #
#########################################
# add_subdirectory(fluid_sph)

#########################################
#         Position-based Fluids        #
#########################################
add_subdirectory(fluid_pbf)

#########################################
#         Fluid-Particle in Cell        #
#########################################
# add_subdirectory(fluid_flip)

#########################################
#             Eulerian Fluid            #
#########################################
# add_subdirectory(fluid_eulerian)
