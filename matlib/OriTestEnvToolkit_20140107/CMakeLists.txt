cmake_minimum_required(VERSION 2.6)
project(oritestenv_MPI)

#find_package(MPI REQUIRED)
#include_directories(${MPI_INCLUDE_PATH})
add_executable(oritestenv 
	../src/applic.cpp 
	../src/io.cpp 
	../src/main.cpp 
	../src/mymath.cpp 
	../src/random.cpp
)
target_link_libraries(oritestenv ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(oritestenv PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(oritestenv PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()