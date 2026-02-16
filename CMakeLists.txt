cmake_minimum_required(VERSION 3.5)
project(InmostVoro)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
  add_compile_options(/O2 /Oi /Ot)
else()
  add_compile_options(-O3 -march=native -ffast-math -funroll-loops)
endif()

find_package(inmost REQUIRED)

add_library(VoroLib
		voronoi_builder.cpp
		helpers.cpp
		voroqh/voronoi.cpp
)

target_include_directories(VoroLib
		PUBLIC ${INMOST_INCLUDE_DIRS}
		PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/voroqh
)

target_link_libraries(VoroLib
		PUBLIC ${INMOST_LIBRARIES}
)

target_compile_definitions(VoroLib
		PUBLIC ${INMOST_DEFINITIONS}
)

set_property(TARGET VoroLib PROPERTY PUBLIC_HEADER voronoi_builder.h)


add_executable(Voro main.cpp)
target_include_directories(Voro PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/voroqh)
target_link_libraries(Voro PRIVATE VoroLib)

add_executable(Benchmark benchmark.cpp)
target_include_directories(Benchmark PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/voroqh)
target_link_libraries(Benchmark PRIVATE VoroLib)


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