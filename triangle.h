#ifndef		__TRINAGLE_H__
#define 	__TRINAGLE_H__

#include "Eigen/Dense"
using namespace Eigen;

struct Triangle
{
	Vector3d 	points[3];
	int additionaldata;
	vector<int> adjacentGrainIDs;
};

#endif		//__TRIANGLE_H__
