cmake_minimum_required(VERSION 2.6)


project(LevelSet_IMM CXX)

set(CMAKE_BUILD_DIR "build")

set(MY_SOURCES box.cpp functions.cpp levelsetproject.cpp matrix.cpp outOfBoundsException.cpp utilities.cpp vektor.cpp)

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
set(CMAKE_CXX_FLAGS "-fno-stack-protector")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}") 

set(VORO++_FOUND true)
set(VORO++_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/voro++/)
set(VORO++_LIBRARIES ${CMAKE_SOURCE_DIR}/voro++/lib/libvoro++.a)
set(VORO++_DEFINITIONS "")
#die n�chste Zeile sollte eigentlich funktionieren, tut es aber nicht, da die lokale Installation von voro++ die FindVoro++.cmake unn�tigt macht,und die Pfade lokale zum Source bleiben (jaja, kann man �ndern).
#find_package(Voro++ REQUIRED)


# dies hier hinzuf�gen...
#find_package(Boost COMPONENTS serialization mpi REQUIRED)
#find_package(MPI REQUIRED)


find_package(Boost REQUIRED)


ADD_DEFINITIONS("-std=c++0x")

include_directories(${Boost_INCLUDE_PATH})
include_directories(${VORO++_INCLUDE_DIRS})
include_directories(${MPI_INCLUDE_PATH})


add_executable(LevelSet_IMM ${MY_SOURCES})

target_link_libraries(LevelSet_IMM ${MPI_LIBRARIES}  ${VORO++_LIBRARIES} ${Boost_LIBRARIES} -lfftw3 -lm)



ADD_CUSTOM_TARGET (clearAll
	COMMAND ${CMAKE_COMMAND} -E remove *.png *.gnu
	COMMAND ${CMAKE_COMMAND} -E echo "rm *.png *.gnu"
)
