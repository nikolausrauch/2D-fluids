#########################################
#          2D OpenGL2 Viewer            #
#---------------------------------------#
# @author  Nikolaus Rauch               #
# @date    05.2018                      #
#########################################
project( 2D-Viewer )

#################################
#          Libraries            #
#################################
find_package( Threads REQUIRED )

set( OpenGL_GL_PREFERENCE GLVND )
find_package( OpenGL 3.0 REQUIRED )

#################################
#     Viewer Source Files       #
#################################
set( VIEWER_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/keyboard.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/mouse.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/window.cpp"
    )

set( VIEWER_HDR
    "${CMAKE_CURRENT_SOURCE_DIR}/keyboard.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/mouse.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/window.h"

    "${CMAKE_CURRENT_SOURCE_DIR}/base.h"
    )


set( VIEWER_SRC viewer.cpp          ${VIEWER_SRC} )
set( VIEWER_HDR viewer.h viewer.inl ${VIEWER_HDR} )

source_group( TREE ${CMAKE_CURRENT_SOURCE_DIR}
    FILES ${VIEWER_SRC} ${VIEWER_HDR} )


#################################
#       Build Viewer Library    #
#################################
add_library( viewer_static STATIC ${VIEWER_SRC} ${VIEWER_HDR} )

target_link_libraries( viewer_static Threads::Threads OpenGL::GL glfw glm_static imgui_static stb_image )

target_include_directories( viewer_static PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
    )

target_compile_features( viewer_static PUBLIC cxx_std_17 )
set_target_properties( viewer_static PROPERTIES CXX_EXTENSIONS OFF )
