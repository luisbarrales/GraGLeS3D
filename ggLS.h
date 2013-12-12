#ifndef GGLS_h
#define GGLS_h

#define DIM 2

#define M 200 //gridpoints in each direction
#define EPS 1e-6
#define BORDER 6
#define DELTA BORDER * 1/double(M)
// #define _USE_MATH_DEFINES
#define PI 3.14159265358979323846
#define INTERIMVAL -sqrt(1. / PARTICLES)

#define TIMESTEPS 10
#define PRINTSTEP 1
#define ANALYSESTEP 1
#define MODE 2 // 2 für lesen;  1 für erzeugen der mikrostrukture
#define NDEBUG

#define PARTICLES 50

#define ISOTROPIC true

#define SAVEIMAGE false
#define IMAGEOUT true

#define SAVEREDIST true
#define SAVECOMP true
#define SAVECONV true

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

#include "vektor.h"
#include "outOfBoundsException.h"
#include "box.h"
#include "weightmap.h"
#include "grainhdl.h"
#include "mymath.h"

#endif

