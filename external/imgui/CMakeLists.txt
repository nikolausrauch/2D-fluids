set(IMGUI_SRC
	"${CMAKE_CURRENT_SOURCE_DIR}/imgui.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/imgui_demo.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/imgui_draw.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/imgui_stdlib.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/imgui_tables.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/imgui_widgets.cpp"
	
	"${CMAKE_CURRENT_SOURCE_DIR}/implot.cpp"
	"${CMAKE_CURRENT_SOURCE_DIR}/implot_items.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/implot_demo.cpp"
	)
	
set(IMGUI_HDR
	"${CMAKE_CURRENT_SOURCE_DIR}/fontdefs.h"
	"${CMAKE_CURRENT_SOURCE_DIR}/IconsFontAwesome.h"
	"${CMAKE_CURRENT_SOURCE_DIR}/imconfig.h"
	"${CMAKE_CURRENT_SOURCE_DIR}/imgui.h"
	"${CMAKE_CURRENT_SOURCE_DIR}/imgui_internal.h"
	"${CMAKE_CURRENT_SOURCE_DIR}/imgui_stdlib.h"
	"${CMAKE_CURRENT_SOURCE_DIR}/imstb_rectpack.h"
	"${CMAKE_CURRENT_SOURCE_DIR}/imstb_textedit.h"
	"${CMAKE_CURRENT_SOURCE_DIR}/imstb_truetype.h"
	
	"${CMAKE_CURRENT_SOURCE_DIR}/implot.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/implot_internal.h"
	)
	
source_group( TREE ${CMAKE_CURRENT_SOURCE_DIR}
    FILES ${IMGUI_SRC} ${IMGUI_HDR} )

add_library(imgui_static STATIC ${IMGUI_SRC} ${IMGUI_HDR} )
target_include_directories( imgui_static PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    )
    
target_compile_features( imgui_static PUBLIC cxx_std_17 )
set_target_properties( imgui_static PROPERTIES CXX_EXTENSIONS OFF )
