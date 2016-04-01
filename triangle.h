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
		barycenter[0]=(points[0][0]+points[1][0]+points[2][0])/3;
		barycenter[1]=(points[0][1]+points[1][1]+points[2][1])/3;
		barycenter[2]=(points[0][2]+points[1][2]+points[2][2])/3;
		//TODO;
		return barycenter;
	}
	Vector3d computeOuterUnitNormal(){
		Vector3d outerUnitNormal;
		Vector3d R;
		Vector3d S;
		R[1]=points[0][0]-points[1][0]; R[1]=points[0][1]-points[1][1]; R[3]=points[0][2]-points[1][2];
		S[1]=points[0][0]-points[2][0]; S[1]=points[0][1]-points[2][1]; S[3]=points[0][2]-points[2][2];
		 outerUnitNormal[0]=R[2]*S[3]-R[3]*S[2];
		 outerUnitNormal[1]=R[3]*S[1]-R[1]*S[3];
		 outerUnitNormal[2]=R[1]*S[2]-R[2]*S[1];
		//TODO:
		return outerUnitNormal;
	}
};

#endif		//__TRIANGLE_H__
