# Matrix makefile
#
# Author : Christian Miessen
# Date   : Juli 2012

CPPSRC =   utilities.cpp outOfBoundsException.cpp vektor.cpp matrix.cpp box.cpp functions.cpp levelsetproject.cpp  
CPPOBJ = $(CPPSRC:%.cpp=%.o)

PROG_NAME = levelset_project
#CXX=g++
INCLUDES = -I./vorolib/include/voro++
LINKER_FLAGS = -L./vorolib/lib 
C_FLAGS=-Wall -std=c++0x -ansi -pedantic -O3 -g



all: ${CPPOBJ}
	 ${CXX} -o ${PROG_NAME}  ${INCLUDES} ${LINKER_FLAGS} ${CPPOBJ} -lvoro++ -lfftw3 -lm -fno-stack-protector


%.o: %.cpp
	 ${CXX} ${C_FLAGS} ${INCLUDES} ${LINKER_FLAGS} -c $<


clean: 
	rm -f *.o


clearAll:
	rm -f *.gnu
	rm -f *.png

