#ifndef		__TRINAGLE_H__
#define 	__TRINAGLE_H__

#include "Eigen/Dense"
using namespace Eigen;

struct Triangle
{
	Vector3d 	points[3];
	int additionalData;
	Vector3d computeBarycenter(){
		Vector3d barycenter;
		//TODO;
		return barycenter;
	}
	Vector3d computeOuterUnitNormal(){
		Vector3d outerUnitNormal;
		//TODO:
		return outerUnitNormal;
	}
};

#endif		//__TRIANGLE_H__
