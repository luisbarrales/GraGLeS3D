#ifndef LEVELSETPROJECT_h
#define LEVELSETPROJECT_h

#define DIM 2

#define M 500
#define EPS 1e-6
#define DELTA 5 * 1/double(M)
// #define _USE_MATH_DEFINES
#define PI 3.14159265358979323846

#define TIMESTEPS 0
#define PRINTSTEP 5
#define PARTICLES 50
#define INTERIMVAL -sqrt(1. / PARTICLES)


#define SAFEFILES true
#define IMAGEOUT true
#define PLOTGNU true
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

#include "./vorolib/include/voro++/voro++.hh"
#include "utilities.h"
#include "outOfBoundsException.h"
#include "vektor.h"
#include "box.h"
#include "matrix.h"
#include "functions.h"


#endif

