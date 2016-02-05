cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(star-vertex)


# Add to path in order to pick up the FindXXX.cmake files included in this project
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/star-cmake")

# This project depends on ROOT
find_package(ROOT)

if(NOT ROOT_FOUND)
	message(FATAL_ERROR "Fatal error: ROOT package not found")
endif()

# Perform some setup standard to STAR experiment environment
include(StarCommon)

add_definitions(-D__ROOT__)


include_directories(
	${ROOT_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	${STAR_BASE_DIR}
	${STAR_STEVENT_DIR}
	${STAR_STDB_DIR}
	${STAR_TABLES_DIR}
	${STAR_EMC_DIR}
	${STAR_DAQ_DIR}
	${STAR_TPC_DIR}
	${STAR_TOF_DIR}
	${STAR_PAMS_DIR}/global/inc
	${STAR_STI_DIR}
	${STAR_BFCHAIN_DIR}
)


# Build StGenericVertexMaker library
root_generate_linkdef_and_dictionary(StGenericVertexMaker HEADERS
	StGenericVertexMaker/StGenericVertexMaker.h StGenericVertexMaker/Minuit/St_VertexCutsC.h)

file(GLOB StGenericVertexMaker_SOURCES "StGenericVertexMaker/*.cxx"  "StGenericVertexMaker/Minuit/*.cxx"
	"StGenericVertexMaker/StiPPVertex/*.cxx" "StGenericVertexMaker/StvPPVertex/*.cxx")

add_library(StGenericVertexMaker SHARED ${StGenericVertexMaker_SOURCES} StGenericVertexMaker_dict.cxx)


# Installation section
install(TARGETS StGenericVertexMaker
	DESTINATION "${STAR_ADDITIONAL_INSTALL_PREFIX}/lib")
