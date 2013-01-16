# Matrix makefile
#
# Author : Christian Miessen
# Date   : Juli 2012

CPPSRC =   utilities.cpp outOfBoundsException.cpp vektor.cpp matrix.cpp functions.cpp levelsetproject.cpp  
CPPOBJ = $(CPPSRC:%.cpp=%.o)

PROG_NAME = levelset_project
#CXX=g++
INCLUDES = -I./vorolib/include/voro++
LINKER_FLAGS = -L./vorolib/lib 
C_FLAGS=-Wall -ansi -pedantic -O3 -g



all: ${CPPOBJ}
	 ${CXX} -o ${PROG_NAME}  ${INCLUDES} ${LINKER_FLAGS} ${CPPOBJ} -lvoro++ -lfftw3 -lm -fno-stack-protector


%.o: %.cpp
	 ${CXX} ${C_FLAGS} ${INCLUDES} ${LINKER_FLAGS} -c $<


clean: 
	rm -f *.o

<<<<<<< HEAD
remove:
	rm *.gnu
	rm *.png
=======
clearAll:
	rm -f *.gnu
	rm -f *.png
>>>>>>> b101026bcaecc9e071efe64ea32bb37e7518099c
