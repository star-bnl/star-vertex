# cmake -D CMAKE_CXX_FLAGS=-std=c++0x \
# -D STAR_BASE_DIR=~/star-repos/star-base \
# -D STAR_STEVENT_DIR=~/star-repos/star-stevent \
# -D STAR_STDB_DIR=~/star-repos/star-stdb \
# -D STAR_TABLES_DIR=~/star-repos/star-tables \
# -D STAR_EMC_DIR=~/star-repos/star-emc \
# -D STAR_DAQ_DIR=~/star-repos/star-daq \
# -D STAR_TPC_DIR=~/star-repos/star-tpc \
# -D STAR_TOF_DIR=~/star-repos/star-tof \
# -D STAR_PAMS_DIR=~/star-repos/star-pams \
# -D STAR_STI_DIR=~/star-repos/star-sti \
# -D STAR_BFCHAIN_DIR=~/star-soft/star-bfchain ../

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(star-vertex)


# Add to path in order to pick up the FindXXX.cmake files included in this project
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# This project depends on ROOT
find_package(ROOT)

if(NOT ROOT_FOUND)
	message(FATAL_ERROR "Fatal error: ROOT package not found")
endif()


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


# StGenericVertexMaker library
#root_generate_dictionary(StGenericVertexMaker_dict
#	StGenericVertexMaker/StGenericVertexMaker.h StGenericVertexMaker/Minuit/St_VertexCutsC.h
#	LINKDEF StGenericVertexMaker/StGenericVertexMakerLinkDef.h OPTIONS -D__ROOT__)
root_generate_linkdef_and_dictionary(StGenericVertexMaker HEADERS
	StGenericVertexMaker/StGenericVertexMaker.h StGenericVertexMaker/Minuit/St_VertexCutsC.h)

file(GLOB StGenericVertexMaker_SOURCES "StGenericVertexMaker/*.cxx"  "StGenericVertexMaker/Minuit/*.cxx"
	"StGenericVertexMaker/StiPPVertex/*.cxx" "StGenericVertexMaker/StvPPVertex/*.cxx")

add_library (StGenericVertexMaker SHARED ${StGenericVertexMaker_SOURCES} StGenericVertexMaker_Cint.cxx)


# Installation section
set(ADDITIONAL_INSTALL_PREFIX ".")

if(DEFINED ENV{STAR_HOST_SYS})
	set(ADDITIONAL_INSTALL_PREFIX ".$ENV{STAR_HOST_SYS}")
endif()

install(TARGETS StGenericVertexMaker
	DESTINATION "${ADDITIONAL_INSTALL_PREFIX}/lib")
