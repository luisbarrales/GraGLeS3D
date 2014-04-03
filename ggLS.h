#ifndef GGLS_h
#define GGLS_h



// #define _USE_MATH_DEFINES
#define PI 3.14159265358979323846
#define INTERIMVAL -sqrt(1. / PARTICLES)

#define IMAGEOUT true

// // switch between double or single precision!
#define PRECISION 0//0= double 1=single

#if PRECISION >0
  #define fftwp_malloc fftwf_malloc
  #define dataprecision float
  #define DimensionalBufferVar DimensionalBufferReals
  #define PooledDimensionalBuffer PooledDimensionalBufferReal
  #define fftwp_complex fftwf_complex
  #define fftwp_plan fftwf_plan
  #define fftwp_free fftwf_free
#elif PRECISION <1
  #define DimensionalBufferVar DimensionalBufferDouble
  #define fftwp_complex fftw_complex
  #define fftwp_plan fftw_plan
  #define fftwp_free fftw_free
  #define fftwp_malloc fftw_malloc
  #define PooledDimensionalBuffer PooledDimensionalBufferDouble
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

