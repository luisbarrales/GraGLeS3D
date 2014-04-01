#ifndef GGLS_h
#define GGLS_h



// #define _USE_MATH_DEFINES
#define PI 3.14159265358979323846
#define INTERIMVAL -sqrt(1. / PARTICLES)

#define IMAGEOUT true

// switch if double or single precision!
#define PRECISION single//double 
#define DimensionalBufferVar DimensionalBuffer<double>
#define dataprecision double
#define fftwp_complex fftw_complex
#define fftwp_plan fftw_plan
#define PooledDimensionalBuffer PooledDimensionalBufferDouble

#if PRECISION == single
  #undef dataprecision
  #undef DimensionalBufferVar
  #undef PooledDimensionalBuffer
  #undef  fftwp_complex
  #undef  fftwp_plan 
  
  #define dataprecision float
  #define DimensionalBufferVar DimensionalBufferReal
  #define PooledDimensionalBuffer PooledDimensionalBufferReal
  #define fftwp_complex fftwf_complex
  #define fftwp_plan fftwf_plan
#endif


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
#include "box.h"
#include "weightmap.h"
#include "grainhdl.h"
#include "mymath.h"
#include "random.h"
#include "applic.h"
#include "io.h"
#include "dimensionalBuffer.h"
#include "marchingSquares.h"

#endif

