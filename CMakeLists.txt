#########################################
#         2D Fluid Simulations          #
#---------------------------------------#
# @author Nikolaus Rauch                #
# @date   04.11.2023                    #
#---------------------------------------#
#########################################
cmake_minimum_required( VERSION 3.9 )
project( 2D-Viewer )

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


#########################################
#              Build Viewer             #
#########################################
add_subdirectory(viewer)


#########################################
#      Build Fluid Dynamic Solvers      #
#########################################
# add_subdirectory(fluid)

