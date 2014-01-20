cmake_minimum_required(VERSION 2.6)


project(LevelSet_IMM C CXX)

set(CMAKE_BUILD_DIR "build")

set(MY_SOURCES vektor.cpp mymath.cpp io.cpp random.cpp applic.cpp grainhdl.cpp box.cpp ggLS.cpp domainCl.cpp utilities.cpp weightmap.cpp marchingSquares.cpp)

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
set(CMAKE_CXX_FLAGS "-g")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}") 

set(VORO++_FOUND true)
set(VORO++_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/voro++/)
set(VORO++_LIBRARIES ${CMAKE_SOURCE_DIR}/voro++/lib/libvoro++.a)
set(VORO++_DEFINITIONS "-g -O2")


# OPEN MP
FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#die nächste Zeile sollte eigentlich funktionieren, tut es aber nicht, da die lokale Installation von voro++ die FindVoro++.cmake unnötigt macht,und die Pfade lokale zum Source bleiben (jaja, kann man ändern).
#find_package(Voro++ REQUIRED)


#INCLUDE(FindBoost.cmake)
#find_package(Boost COMPONENTS serialization mpi REQUIRED)
#find_package(Boost COMPONENTS serialization mpi REQUIRED)

#include_directories(${Boost_INCLUDE_PATH})
include_directories(${VORO++_INCLUDE_DIRS})
#include_directories(${MPI_INCLUDE_PATH})

link_directories(${Boost_LIBRARY_DIR})
# dies hier hinzufügen...

#find_package(MPI REQUIRED)


ADD_DEFINITIONS("-std=c++0x")



add_executable(LevelSet_IMM ${MY_SOURCES})

target_link_libraries(LevelSet_IMM ${MPI_LIBRARIES}  ${VORO++_LIBRARIES} ${Boost_LIBRARIES} -lfftw3 -lm)


ADD_CUSTOM_TARGET (clearAll
	COMMAND ${CMAKE_COMMAND} -E remove *.png *.gnu *.ori Level* *.txt *.gif
	COMMAND ${CMAKE_COMMAND} -E echo "rm *.png *.gnu *.ori *.gif"
)

ADD_CUSTOM_TARGET (rm_cmake
	COMMAND ${CMAKE_COMMAND} -E remove CMakeFiles/
)

ADD_CUSTOM_TARGET (doc
	COMMAND doxygen
)
