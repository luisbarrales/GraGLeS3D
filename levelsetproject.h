#ifndef LEVELSETPROJECT_h
#define LEVELSETPROJECT_h

#define DIM 2

#define M 100 //gridpoints in each direction
#define EPS 1e-6
#define DELTA 5 * 1/double(M)
// #define _USE_MATH_DEFINES
#define PI 3.14159265358979323846

#define TIMESTEPS 100
#define PRINTSTEP 5
#define PARTICLES 50

#define INTERIMVAL -sqrt(1. / PARTICLES)

#define SAFEFILES true
#define IMAGEOUT true
#define PLOTGNU false
#define DISCRETE_CONVOLUTION false

#pragma once
#include <vector>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip> 
#include <list>
#include <sstream>
#include <string.h>
#include <algorithm>
#include <fftw3.h>
#include <time.h>  
// #include <boost/mpi.hpp>
#include <omp.h>
// #include <mpi.h>
// #include <boost/mpi/environment.hpp>
// #include <boost/mpi/communicator.hpp>

#include "voro++/include/voro++/voro++.hh"

#include "utilities.h"
#include "outOfBoundsException.h"
#include "vektor.h"
#include "box.h"
#include "matrix.h"
#include "functions.h"


#endif

