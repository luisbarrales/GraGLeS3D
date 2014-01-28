#ifndef GGLS_h
#define GGLS_h

#define DIM 2
#define PARTICLES 50
 //gridpoints in each direction

#define BORDER 6

// #define _USE_MATH_DEFINES
#define PI 3.14159265358979323846
#define INTERIMVAL -sqrt(1. / PARTICLES)

#define TIMESTEPS 100
#define ANALYSESTEP 10
#define DISCRETESAMPLING 25

#define MODE 1 // 2 for read Microstructure;  1 for use Voro++
#define NDEBUG

#define HAGB 0.6;

#define ISOTROPIC false
#define TEXTURE true 	// generates a Texture round a bunge orientation whith a deviation -> see grainhandler
#define SAVEIMAGE false
#define IMAGEOUT true


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

