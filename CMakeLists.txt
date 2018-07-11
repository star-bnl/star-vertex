cmake_minimum_required(VERSION 3.6 FATAL_ERROR)

project(star-vertex LANGUAGES CXX)

# Add to path in order to pick up the FindXXX.cmake files included in this project
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# This project depends on ROOT
find_package(ROOT)

if( ${ROOT_CXX_FLAGS} MATCHES "-m32" )
	message(STATUS "Found -m32 option in $ROOT_CXX_FLAGS (root-config). Will add it to $CMAKE_CXX_FLAGS")
	set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS FALSE)
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32" )
endif()

find_package(STAR)

set(STAR_SRC ${CMAKE_CURRENT_SOURCE_DIR})

# Perform some setup standard to STAR experiment environment
include(StarCommon)

include_directories(
	${ROOT_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	${STAR_INCLUDE_DIRS}
)

star_add_library(StGenericVertexMaker)
star_add_library_vertexnosti(StGenericVertexMaker)
star_add_library(StSecondaryVertexMaker)
star_add_library(StVertexRootIO)
star_add_library(StZdcVertexMaker)
