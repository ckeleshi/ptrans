cmake_minimum_required( VERSION 3.0 )

# CloudCompare example for standard plugins

# REPLACE ALL 'ExamplePlugin' OCCURENCES BY YOUR PLUGIN NAME
# AND ADAPT THE CODE BELOW TO YOUR OWN NEEDS!

# Add an option to CMake to control whether we build this plugin or not
option( PLUGIN_COLORIMETRIC_SEGMENTER "Check to install segmenter plugin" OFF )

if ( PLUGIN_COLORIMETRIC_SEGMENTER )
	# Name the plugin
    project( ColorimetricSegmenter )

    	# we need includes from the main CC source dir
	include_directories( ${CloudCompare_SOURCE_DIR} )

    list( APPEND CC_PLUGIN_CUSTOM_HEADER_LIST ${CloudCompare_SOURCE_DIR}/../common/ccPickingListener.h )
    list( APPEND CC_PLUGIN_CUSTOM_HEADER_LIST ${CloudCompare_SOURCE_DIR}/../common/ccPickingHub.h )
	list( APPEND CC_PLUGIN_CUSTOM_SOURCE_LIST ${CloudCompare_SOURCE_DIR}/../common/ccPickingHub.cpp )
    
    # load any subdirectories (see qAdditionalIO for an example)
    # add_subdirectory( LIB1 )
    
    include( ../../CMakePluginTpl.cmake )
    
    # set dependencies to necessary libraries (see qPCV for an example)
    # target_link_libraries( ${PROJECT_NAME} LIB1 )
    # include_directories( ${LIB1_INCLUDE_DIR} )
endif()
