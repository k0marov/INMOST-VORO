cmake_minimum_required(VERSION 3.5)
project(InmostVoro)

find_package(inmost REQUIRED)

add_library(VoroLib
		planar_cutter.cpp
		voronoi_builder.cpp
		helpers.cpp
)

target_include_directories(VoroLib
		PUBLIC ${INMOST_INCLUDE_DIRS}
)

target_link_libraries(VoroLib
		PUBLIC ${INMOST_LIBRARIES}
)

target_compile_definitions(VoroLib
		PUBLIC ${INMOST_DEFINITIONS}
)

set_property(TARGET VoroLib PROPERTY PUBLIC_HEADER voronoi_builder.h)


add_executable(Voro main.cpp)
target_link_libraries(Voro PRIVATE VoroLib)


# Installation
include(GNUInstallDirs)

install(TARGETS Voro
		RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(TARGETS VoroLib EXPORT inmost-targets
		LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
		ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
		PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT inmost-targets
		FILE InmostVoroTargets.cmake
		NAMESPACE InmostVoro::
		DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/InmostVoro
)