#[[
########################################################################################
##
##  This file is part of Sandpile Simulator, a cellular automaton for sandpile dynamics.
##  Copyright (C) 2021, 2024 M. Frohne
##
##  Sandpile Simulator is free software: you can redistribute it and/or modify
##  it under the terms of the GNU Affero General Public License as published
##  by the Free Software Foundation, either version 3 of the License,
##  or (at your option) any later version.
##
##  Sandpile Simulator is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##  See the GNU Affero General Public License for more details.
##
##  You should have received a copy of the GNU Affero General Public License
##  along with Sandpile Simulator. If not, see <https://www.gnu.org/licenses/>.
##
########################################################################################
#]]

cmake_minimum_required(VERSION 3.10.1)

include(GNUInstallDirs)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose type of build (Debug or Release)")
endif(NOT CMAKE_BUILD_TYPE)

set(SANDSIM_ENABLE_PLOTTER 0 CACHE BOOL "Build plotting class. Requires SFML 2.x!")
set(SANDSIM_BUILD_DOCUMENTATION 0 CACHE BOOL "Build documentation with Doxygen.")

project(SandpileSimulator VERSION 0.2.0)

set(EXECUTABLE_NAME ${PROJECT_NAME})

set(BUILD_SHARED_LIBS ON)

if(SANDSIM_BUILD_DOCUMENTATION)
    find_package(Doxygen)
endif()

if(SANDSIM_ENABLE_PLOTTER)
    find_package(SFML 2 REQUIRED graphics system window)
endif()

set(FILENAMES
    argumentparser
    aux
    avalanche
    avalanchestatistics
    logger
    momentanalysis
    noncopyable
    randomizer
    sandbox
    simulationmanager
    simulationmodel
    simulationmodel_btw
    simulationmodel_fwm
    )
set(FILENAMES_sfReq
    plotter
    )
set(FILENAMES_sfMod
    simulator
    )
set(FILENAMES_hdOnl
    version
    )

foreach(filename ${FILENAMES})
    set(SOURCES "${SOURCES}" "${PROJECT_SOURCE_DIR}/src/sources/${filename}.cpp")
endforeach(filename)

foreach(filename ${FILENAMES_sfReq})
    set(FILENAMES ${FILENAMES} ${filename})

    if(SANDSIM_ENABLE_PLOTTER)
        set(SOURCES "${SOURCES}" "${PROJECT_SOURCE_DIR}/src/sources/${filename}.cpp")
    endif(SANDSIM_ENABLE_PLOTTER)
endforeach(filename)

foreach(filename ${FILENAMES})
    set(HEADERS "${HEADERS}" "${filename}.h")
    configure_file("${PROJECT_SOURCE_DIR}/src/headers/${filename}.h" "${PROJECT_BINARY_DIR}/include/${filename}.h" COPYONLY)
endforeach(filename)

foreach(filename ${FILENAMES_hdOnl})
    set(HEADERS "${HEADERS}" "${filename}.h")
    configure_file("${PROJECT_SOURCE_DIR}/src/headers/${filename}.h" "${PROJECT_BINARY_DIR}/include/${filename}.h" COPYONLY)
endforeach(filename)

#Always include these, but headers need to be configured
foreach(filename ${FILENAMES_sfMod})
    set(HEADERS "${HEADERS}" "${filename}.h")
    set(SOURCES "${SOURCES}" "${PROJECT_SOURCE_DIR}/src/sources/${filename}.cpp")
    configure_file("${PROJECT_SOURCE_DIR}/src/headers/${filename}.in.h" "${PROJECT_BINARY_DIR}/include/${filename}.h" @ONLY)
endforeach(filename)

set(LIB_SOURCES "${SOURCES}")
set(SOURCES "${SOURCES}" "${PROJECT_SOURCE_DIR}/src/sources/main.cpp")

add_executable("${EXECUTABLE_NAME}-bin" ${SOURCES})
add_library("${EXECUTABLE_NAME}-lib" ${LIB_SOURCES})

set_target_properties("${EXECUTABLE_NAME}-bin" "${EXECUTABLE_NAME}-lib" PROPERTIES OUTPUT_NAME ${EXECUTABLE_NAME})
set_target_properties("${EXECUTABLE_NAME}-lib" PROPERTIES VERSION ${PROJECT_VERSION})

include_directories("${PROJECT_BINARY_DIR}/include")

set(LINKER_FLAGS_REQUIRED "-pthread -lgslcblas -lgsl -lm")

if(SANDSIM_ENABLE_PLOTTER)
    include_directories(${SFML_INCLUDE_DIR})
    target_link_libraries("${EXECUTABLE_NAME}-bin" ${SFML_LIBRARIES} ${SFML_DEPENDENCIES})
    set(LINKER_FLAGS_REQUIRED "${LINKER_FLAGS_REQUIRED} -lsfml-graphics -lsfml-system -lsfml-window")
endif(SANDSIM_ENABLE_PLOTTER)

set(LINKER_FLAGS_REQUIRED "${LINKER_FLAGS_REQUIRED}" CACHE STRING "Linker flags that are required to link the dependencies." FORCE)

set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_LINKER "g++")

set(CMAKE_CXX_FLAGS_RELEASE "" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS "" CACHE STRING "" FORCE)

set(CXX_FLAGS_RELEASE "-pipe -O2 -march=corei7 -mtune=native -fstack-protector-strong -fno-plt -fPIC -std=c++17 -Wall -Wextra -fopenmp"
        CACHE STRING "Flags used by the compiler for release builds.")
set(CXX_FLAGS_DEBUG "-pipe -O0 -ggdb -march=corei7 -mtune=native -fstack-protector-strong -fno-plt -fPIC -std=c++17 -Wall -Wextra -fopenmp"
        CACHE STRING "Flags used by the compiler for debug builds.")
set(LINKER_FLAGS_RELEASE "-Wl,-O1"
        CACHE STRING "Some flags used by the linker for release builds.")
set(LINKER_FLAGS_DEBUG "-Wl,-O1"
        CACHE STRING "Some flags used by the linker for debug builds.")

set(CMAKE_CXX_FLAGS_RELEASE "${CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_DEBUG "${CXX_FLAGS_DEBUG}")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${LINKER_FLAGS_RELEASE} ${LINKER_FLAGS_REQUIRED}")
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${LINKER_FLAGS_DEBUG} ${LINKER_FLAGS_REQUIRED}")


if(SANDSIM_BUILD_DOCUMENTATION AND DOXYGEN_FOUND)
  set(doxyfile_in "${PROJECT_SOURCE_DIR}/doc/Doxyfile.in")
  set(doxyfile "${PROJECT_BINARY_DIR}/Doxyfile")
  set(TEMPLATE_DIR "${PROJECT_SOURCE_DIR}/doc/template")
  set(HEADERS_DIR "${PROJECT_BINARY_DIR}/include")
  set(SOURCES_DIR "${PROJECT_SOURCE_DIR}/src/sources")
  configure_file("${doxyfile_in}" "${doxyfile}" @ONLY)
  add_custom_target(Documentation ALL ${DOXYGEN_EXECUTABLE} "${doxyfile}" WORKING_DIRECTORY "${PROJECT_BINARY_DIR}" VERBATIM)
endif()

install(TARGETS "${EXECUTABLE_NAME}-bin" RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(TARGETS "${EXECUTABLE_NAME}-lib" LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY "${PROJECT_BINARY_DIR}/include/" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}")
if(SANDSIM_BUILD_DOCUMENTATION AND DOXYGEN_FOUND)
    install(DIRECTORY "${PROJECT_BINARY_DIR}/doc/" DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/doc")
endif()
