#ifndef GGLS_h
#define GGLS_h

#define DIM 2

#define M 200 //gridpoints in each direction
#define EPS 1e-6
#define DELTA 5 * 1/double(M)
// #define _USE_MATH_DEFINES
#define PI 3.14159265358979323846
#define INTERIMVAL -sqrt(1. / PARTICLES)

#define TIMESTEPS 1000
#define PRINTSTEP 1000
#define ANALYSESTEP 100


#define PARTICLES 1000

#define ISOTROPIC false

#define SAVEIMAGE false
#define IMAGEOUT true

#define SAVEREDIST true
#define SAVECOMP false
#define SAVECONV false

#pragma once
#include <map>
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
#include <omp.h>

// g++ -I/home/user/bin/R-devel/include -L/home/user/bin/R-devel/lib -lR embed.cpp

#include "voro++/include/voro++/voro++.hh"
#include "utilities.h"
#include "outOfBoundsException.h"
#include "vektor.h"
#include "box.h"c
#include "domainCl.h"
#include "functions.h"
#include "weightmap.h"
#include "grainhdl.h"
#include "mymath.h"

#endif

