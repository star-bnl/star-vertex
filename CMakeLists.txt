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
	${STAR_TABLES_DIR}/tables
	${STAR_EMC_DIR}
	${STAR_DAQ_DIR}
	${STAR_TPC_DIR}
	${STAR_TOF_DIR}
	${STAR_PAMS_DIR}/global/inc
	${STAR_STI_DIR}
	${STAR_BFCHAIN_DIR}
	${STAR_MUDST_DIR}
)


# Build StGenericVertexMaker library
root_generate_linkdef_and_dictionary(StGenericVertexMaker HEADERS
	StGenericVertexMaker/StGenericVertexMaker.h StGenericVertexMaker/Minuit/St_VertexCutsC.h)

add_library(StGenericVertexMaker SHARED
	StGenericVertexMaker/StCtbUtility.cxx
	StGenericVertexMaker/StFixedVertexFinder.cxx
	StGenericVertexMaker/StGenericVertexFinder.cxx
	StGenericVertexMaker/StGenericVertexMaker.cxx
	StGenericVertexMaker/VertexTrackFilter.cxx
	StGenericVertexMaker/StppLMVVertexFinder.cxx
	StGenericVertexMaker/VertexFinderOptions.cxx
	StGenericVertexMaker/Minuit/StMinuitVertexFinder.cxx
	StGenericVertexMaker/Minuit/St_VertexCutsC.cxx
	StGenericVertexMaker/StiPPVertex/BemcHitList.cxx
	StGenericVertexMaker/StiPPVertex/BtofHitList.cxx
	StGenericVertexMaker/StiPPVertex/CtbHitList.cxx
	StGenericVertexMaker/StiPPVertex/EemcHitList.cxx
	StGenericVertexMaker/StiPPVertex/ScintHitList.cxx
	StGenericVertexMaker/StiPPVertex/StPPVertexFinder.cxx
	StGenericVertexMaker/StiPPVertex/TrackData.cxx
	StGenericVertexMaker/StiPPVertex/Vertex3D.cxx
	StGenericVertexMaker/StiPPVertex/VertexData.cxx
	StGenericVertexMaker/StvPPVertex/BemcHitList.cxx
	StGenericVertexMaker/StvPPVertex/BtofHitList.cxx
	StGenericVertexMaker/StvPPVertex/CtbHitList.cxx
	StGenericVertexMaker/StvPPVertex/EemcHitList.cxx
	StGenericVertexMaker/StvPPVertex/ScintHitList.cxx
	StGenericVertexMaker/StvPPVertex/StEventToolkit.cxx
	StGenericVertexMaker/StvPPVertex/StPPVertexFinder.cxx
	StGenericVertexMaker/StvPPVertex/TrackData.cxx
	StGenericVertexMaker/StvPPVertex/Vertex3D.cxx
	StGenericVertexMaker/StvPPVertex/VertexData.cxx
	StGenericVertexMaker_dict.cxx
)

# Build StGenericVertexMakerNoSti library
root_generate_linkdef_and_dictionary(StGenericVertexMakerNoSti HEADERS
	StGenericVertexMaker/StGenericVertexMaker.h StGenericVertexMaker/Minuit/St_VertexCutsC.h)

add_library(StGenericVertexMakerNoSti SHARED
	StGenericVertexMaker/StCtbUtility.cxx
	StGenericVertexMaker/StFixedVertexFinder.cxx
	StGenericVertexMaker/StGenericVertexFinder.cxx
	StGenericVertexMaker/StGenericVertexMaker.cxx
	StGenericVertexMaker/StppLMVVertexFinder.cxx
	StGenericVertexMaker/Minuit/StMinuitVertexFinder.cxx
	StGenericVertexMaker/Minuit/St_VertexCutsC.cxx
	StGenericVertexMakerNoSti_dict.cxx
)


# Build StSecondaryVertexMaker library
root_generate_linkdef_and_dictionary(
	StSecondaryVertexMaker
	HEADERS
	StSecondaryVertexMaker/StKinkMaker.h
	StSecondaryVertexMaker/StV0FinderMaker.h
	StSecondaryVertexMaker/StXiFinderMaker.h
	StSecondaryVertexMaker/StKinkLocalTrack.hh
)

add_library(StSecondaryVertexMaker SHARED
	StSecondaryVertexMaker/DecayVertexFinder.cxx
	StSecondaryVertexMaker/StKinkLocalTrack.cc
	StSecondaryVertexMaker/StKinkMaker.cxx
	StSecondaryVertexMaker/StV0FinderMaker.cxx
	StSecondaryVertexMaker/StXiFinderMaker.cxx
	StSecondaryVertexMaker_dict.cxx
)


# Build StVertexRootIO library
root_generate_linkdef_and_dictionary(
	StVertexRootIO
	HEADERS
	StVertexRootIO/TPrimaryVertex.h
	StVertexRootIO/TDecayVertex.h
	StVertexRootIO/TDecayVertexContainers.h
)

add_library(StVertexRootIO SHARED StVertexRootIO_dict.cxx)


# Build StZdcVertexMaker library
root_generate_linkdef_and_dictionary(StZdcVertexMaker HEADERS StZdcVertexMaker/StZdcVertexMaker.h)

add_library(StZdcVertexMaker SHARED StZdcVertexMaker/StZdcVertexMaker.cxx StZdcVertexMaker_dict.cxx)


# A collective target to build all libraries in this project. Can be used to
# build all specified targets from a parent project
add_custom_target(star-vertex DEPENDS
	StGenericVertexMaker StGenericVertexMakerNoSti StSecondaryVertexMaker StZdcVertexMaker StVertexRootIO
)


# Installation section
install(TARGETS
	StGenericVertexMaker StGenericVertexMakerNoSti StSecondaryVertexMaker StZdcVertexMaker StVertexRootIO
	DESTINATION "${STAR_ADDITIONAL_INSTALL_PREFIX}/lib" OPTIONAL
)
