#ifndef GGLS_h
#define GGLS_h

#define DIM 2

#define M 200 //gridpoints in each direction
#define EPS 1e-6
#define DELTA 5 * 1/double(M)
// #define _USE_MATH_DEFINES
#define PI 3.14159265358979323846


#define TIMESTEPS 1000
#define PRINTSTEP 50
#define PRINTNOW 1000
#define PARTICLES 10

#define ISOTROPIC false
#define TRIPLEPUNKT false
#define FIX_BOUNDARY false

#define INTERIMVAL -sqrt(1. / PARTICLES)

#define SAFEFILES true
#define IMAGEOUT true
#define PLOTGNU true
#define EULER true
#define DOMAINCOMPARISON false
#define DRAW_PARTICLES false

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
#include "box.h"
#include "domainCl.h"
#include "functions.h"
#include "weightmap.h"
#include "grainhdl.h"

#endif

