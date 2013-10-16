
#ifndef _mymath_h_
#define _mymath_h_

#include <math.h>
#include "applic.h"
#include "random.h"

#define SQR(a) ((a)*(a))
#define CUBE(a) ((a)*(a)*(a))
#define MIN(X,Y) ((X) < (Y)) ? (X) : (Y)
#define MAX(X,Y) ((X) > (Y)) ? (X) : (Y)
#define SIGN(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)

typedef double Real;

double dcos(double x);
double dsin(double x);

double qroot(double x);
double kante(double x);
double WLS(double r);

double wicksell(double u, double *srad, long sCount);

float Efsqrt(float k);
void fsqrt2(float *ka, float *kb);
void fsqrt3(float *ka, float *kb, float *kc);
void fsqrt4(float *ka, float *kb, float *kc, float *kd);

double Edsqrt(double k);
void dsqrt2(double *ka, double *kb);
void dsqrt3(double *ka, double *kb, double *kc);

void _fsqrt2(float *ka, float *kb);
void _fsqrt3(float *ka, float *kb, float *kc);

/*void _dsqrt3(double *ka, double *kb, double *kc);
void _dsqrt2(double *ka, double *kb);
double _dsqrt(double a);
float _fsqrt(float a);*/

char inequalitieFulfil( double a, double b, double c, double d );
double misorientationCubic( double pa1, double Pa, double pa2, double pb1, double Pb, double pb2 );
double misorientationCubicQxQ( double q01, double q11, double q21, double q31, double q02, double q12, double q22, double q32 );
void misorientationQuaternionCubic( double* p, double* q, double* quat );
double areaPolygon( int *poly, int counter );
int pnpoly(int nvert, double *vertx, double *verty, double testx, double testy);
void bubbleSort ( Real arr [ ], int size );
void swap ( Real& x, Real& y );
void randomMisorientation( double theta, double* qr  );
void randomMisorientationAxisConsidered( double * qref, double * qr, double maxDev );
void multiplyQuaternions( double *q, double* p, double* r );
void newOrientationFromReference( double *oriOri, double deviation, double *newOri);
void newOrientationFromReferenceFixedAngularCone(double * oriOri, double maxDev,double angle, double u, double v, double w,double * newOri);
double angleBetweenQuaternions( double *, double * );
void rotateOrientation( double *oriOri, double angle, double u, double v, double w, double *newOri );
void quaternion2Euler( double * quat, double * euler );
void euler2quaternion( double * euler, double * q );


/*inline
void _dsqrt3(double *ka, double *kb, double *kc)
{
#if LENDIAN
	*ka = 1 / sqrt( *ka );
	*kb = 1 / sqrt( *kb );
	*kc = 1 / sqrt( *kc );
#else
	double a = *ka;
	double b = *kb;
	double c = *kc;
	
    double half = 0.5, one = 1.0;
    double y01;
	double t01, t11; 
    double y02; 
	double t02, t12;
    double y03; 
	double t03, t13;

    asm volatile("frsqrte %0,%1" : "=f" (y01) : "f" (a));
    asm volatile("frsqrte %0,%1" : "=f" (y02) : "f" (b));
    asm volatile("frsqrte %0,%1" : "=f" (y03) : "f" (c));

    // Do y1 = y0 + 0.5*y0*(1-x*y0*y0)
    asm volatile("fmul %0,%1,%2" : "=f" (t01) : "f" (y01), "f" (y01));
    asm volatile("fmul %0,%1,%2" : "=f" (t02) : "f" (y02), "f" (y02));
    asm volatile("fmul %0,%1,%2" : "=f" (t03) : "f" (y03), "f" (y03));
    asm volatile("fmul %0,%1,%2" : "=f" (t11) : "f" (y01), "f" (half));
    asm volatile("fmul %0,%1,%2" : "=f" (t12) : "f" (y02), "f" (half));
    asm volatile("fmul %0,%1,%2" : "=f" (t13) : "f" (y03), "f" (half));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t01) : "f" (a), "f" (t01), "f" (one));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t02) : "f" (b), "f" (t02), "f" (one));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t03) : "f" (c), "f" (t03), "f" (one));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y01) : "f" (t01), "f" (t11), "f" (y01));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y02) : "f" (t02), "f" (t12), "f" (y02));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y03) : "f" (t03), "f" (t13), "f" (y03));

    asm volatile("fmul %0,%1,%2" : "=f" (t01) : "f" (y01), "f" (y01));
    asm volatile("fmul %0,%1,%2" : "=f" (t02) : "f" (y02), "f" (y02));
    asm volatile("fmul %0,%1,%2" : "=f" (t03) : "f" (y03), "f" (y03));
    asm volatile("fmul %0,%1,%2" : "=f" (t11) : "f" (y01), "f" (half));
    asm volatile("fmul %0,%1,%2" : "=f" (t12) : "f" (y02), "f" (half));
    asm volatile("fmul %0,%1,%2" : "=f" (t13) : "f" (y03), "f" (half));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t01) : "f" (a), "f" (t01), "f" (one));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t02) : "f" (b), "f" (t02), "f" (one));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t03) : "f" (c), "f" (t03), "f" (one));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y01) : "f" (t01), "f" (t11), "f" (y01));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y02) : "f" (t02), "f" (t12), "f" (y02));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y03) : "f" (t03), "f" (t13), "f" (y03));

    asm volatile("fmul %0,%1,%2" : "=f" (t01) : "f" (y01), "f" (y01));
    asm volatile("fmul %0,%1,%2" : "=f" (t02) : "f" (y02), "f" (y02));
    asm volatile("fmul %0,%1,%2" : "=f" (t03) : "f" (y03), "f" (y03));
	asm volatile("fmul %0,%1,%2" : "=f" (t11) : "f" (y01), "f" (half));
	asm volatile("fmul %0,%1,%2" : "=f" (t12) : "f" (y02), "f" (half));
	asm volatile("fmul %0,%1,%2" : "=f" (t13) : "f" (y03), "f" (half));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t01) : "f" (a), "f" (t01), "f" (one));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t02) : "f" (b), "f" (t02), "f" (one));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t03) : "f" (c), "f" (t03), "f" (one));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y01) : "f" (t01), "f" (t11), "f" (y01));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y02) : "f" (t02), "f" (t12), "f" (y02));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y03) : "f" (t03), "f" (t13), "f" (y03));

    asm volatile("fmul %0,%1,%2" : "=f" (t01) : "f" (y01), "f" (y01));
    asm volatile("fmul %0,%1,%2" : "=f" (t02) : "f" (y02), "f" (y02));
    asm volatile("fmul %0,%1,%2" : "=f" (t03) : "f" (y03), "f" (y03));
//	asm volatile("fmul %0,%1,%2" : "=f" (t11) : "f" (y01), "f" (half));
//	asm volatile("fmul %0,%1,%2" : "=f" (t12) : "f" (y02), "f" (half));
//	asm volatile("fmul %0,%1,%2" : "=f" (t13) : "f" (y03), "f" (half));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t01) : "f" (a), "f" (t01), "f" (one));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t02) : "f" (b), "f" (t02), "f" (one));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t03) : "f" (c), "f" (t03), "f" (one));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y01) : "f" (t01), "f" (t11), "f" (y01));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y02) : "f" (t02), "f" (t12), "f" (y02));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y03) : "f" (t03), "f" (t13), "f" (y03));

	*ka = y01;
	*kb = y02;
	*kc = y03;
#endif
}


inline
void _dsqrt2(double *ka, double *kb)
{
#if LENDIAN
	*ka = 1 / sqrt( *ka );
	*kb = 1 / sqrt( *kb );
#else
	double a = *ka;
	double b = *kb;
	
        double half = 0.5, one = 1.0;
        double y01;
        double t01, t11; 
        double y02; 
	double t02, t12;

    asm volatile("frsqrte %0,%1" : "=f" (y01) : "f" (a));
    asm volatile("frsqrte %0,%1" : "=f" (y02) : "f" (b));

    // Do y1 = y0 + 0.5*y0*(1-x*y0*y0)
    asm volatile("fmul %0,%1,%2" : "=f" (t01) : "f" (y01), "f" (y01));
    asm volatile("fmul %0,%1,%2" : "=f" (t02) : "f" (y02), "f" (y02));
    asm volatile("fmul %0,%1,%2" : "=f" (t11) : "f" (y01), "f" (half));
    asm volatile("fmul %0,%1,%2" : "=f" (t12) : "f" (y02), "f" (half));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t01) : "f" (a), "f" (t01), "f" (one));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t02) : "f" (b), "f" (t02), "f" (one));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y01) : "f" (t01), "f" (t11), "f" (y01));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y02) : "f" (t02), "f" (t12), "f" (y02));

    asm volatile("fmul %0,%1,%2" : "=f" (t01) : "f" (y01), "f" (y01));
    asm volatile("fmul %0,%1,%2" : "=f" (t02) : "f" (y02), "f" (y02));
    asm volatile("fmul %0,%1,%2" : "=f" (t11) : "f" (y01), "f" (half));
    asm volatile("fmul %0,%1,%2" : "=f" (t12) : "f" (y02), "f" (half));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t01) : "f" (a), "f" (t01), "f" (one));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t02) : "f" (b), "f" (t02), "f" (one));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y01) : "f" (t01), "f" (t11), "f" (y01));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y02) : "f" (t02), "f" (t12), "f" (y02));

    asm volatile("fmul %0,%1,%2" : "=f" (t01) : "f" (y01), "f" (y01));
    asm volatile("fmul %0,%1,%2" : "=f" (t02) : "f" (y02), "f" (y02));
    asm volatile("fmul %0,%1,%2" : "=f" (t11) : "f" (y01), "f" (half));
    asm volatile("fmul %0,%1,%2" : "=f" (t12) : "f" (y02), "f" (half));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t01) : "f" (a), "f" (t01), "f" (one));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t02) : "f" (b), "f" (t02), "f" (one));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y01) : "f" (t01), "f" (t11), "f" (y01));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y02) : "f" (t02), "f" (t12), "f" (y02));

    asm volatile("fmul %0,%1,%2" : "=f" (t01) : "f" (y01), "f" (y01));
    asm volatile("fmul %0,%1,%2" : "=f" (t02) : "f" (y02), "f" (y02));
//	asm volatile("fmul %0,%1,%2" : "=f" (t11) : "f" (y01), "f" (half));
//	asm volatile("fmul %0,%1,%2" : "=f" (t12) : "f" (y02), "f" (half));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t01) : "f" (a), "f" (t01), "f" (one));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t02) : "f" (b), "f" (t02), "f" (one));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y01) : "f" (t01), "f" (t11), "f" (y01));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y02) : "f" (t02), "f" (t12), "f" (y02));

	*ka = y01;
	*kb = y02;
#endif
}


inline
void _fsqrt2(float *ka, float *kb)
{
#if LENDIAN
	*ka = 1 / sqrt( *ka );
	*kb = 1 / sqrt( *kb );
#else
	float a = *ka;
	float b = *kb;
	
        float half = 0.5f, one = 1.0f;
        float y01;
	float t01, t11; 
        float y02; 
	float t02, t12;

    asm volatile("frsqrte %0,%1" : "=f" (y01) : "f" (a));
    asm volatile("frsqrte %0,%1" : "=f" (y02) : "f" (b));

    // Do y1 = y0 + 0.5*y0*(1-x*y0*y0)
    asm volatile("fmul %0,%1,%2" : "=f" (t01) : "f" (y01), "f" (y01));
    asm volatile("fmul %0,%1,%2" : "=f" (t02) : "f" (y02), "f" (y02));
    asm volatile("fmul %0,%1,%2" : "=f" (t11) : "f" (y01), "f" (half));
    asm volatile("fmul %0,%1,%2" : "=f" (t12) : "f" (y02), "f" (half));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t01) : "f" (a), "f" (t01), "f" (one));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t02) : "f" (b), "f" (t02), "f" (one));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y01) : "f" (t01), "f" (t11), "f" (y01));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y02) : "f" (t02), "f" (t12), "f" (y02));

    asm volatile("fmul %0,%1,%2" : "=f" (t01) : "f" (y01), "f" (y01));
    asm volatile("fmul %0,%1,%2" : "=f" (t02) : "f" (y02), "f" (y02));
    asm volatile("fmul %0,%1,%2" : "=f" (t11) : "f" (y01), "f" (half));
    asm volatile("fmul %0,%1,%2" : "=f" (t12) : "f" (y02), "f" (half));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t01) : "f" (a), "f" (t01), "f" (one));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t02) : "f" (b), "f" (t02), "f" (one));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y01) : "f" (t01), "f" (t11), "f" (y01));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y02) : "f" (t02), "f" (t12), "f" (y02));

	*ka = y01;
	*kb = y02;

#endif
}


inline
float _fsqrt(float a)
{
#if LENDIAN
	return 1 / sqrt( a );
#else
    double half = 0.5, one = 1.0;
    double y01;
	double t01, t11; 

    asm volatile("frsqrte %0,%1" : "=f" (y01) : "f" (a));

    asm volatile("fmul %0,%1,%2" : "=f" (t01) : "f" (y01), "f" (y01));
    asm volatile("fmul %0,%1,%2" : "=f" (t11) : "f" (y01), "f" (half));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t01) : "f" (a), "f" (t01), "f" (one));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y01) : "f" (t01), "f" (t11), "f" (y01));

    asm volatile("fmul %0,%1,%2" : "=f" (t01) : "f" (y01), "f" (y01));
	asm volatile("fmul %0,%1,%2" : "=f" (t11) : "f" (y01), "f" (half));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t01) : "f" (a), "f" (t01), "f" (one));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y01) : "f" (t01), "f" (t11), "f" (y01));

	return y01;
#endif
}

inline
double _dsqrt(double a)
{
#if LENDIAN
	return 1 / sqrt( a );
#else
    double half = 0.5, one = 1.0;
    double y01;
    double t01, t11; 

    asm volatile("frsqrte %0,%1" : "=f" (y01) : "f" (a));

    asm volatile("fmul %0,%1,%2" : "=f" (t01) : "f" (y01), "f" (y01));
    asm volatile("fmul %0,%1,%2" : "=f" (t11) : "f" (y01), "f" (half));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t01) : "f" (a), "f" (t01), "f" (one));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y01) : "f" (t01), "f" (t11), "f" (y01));

    asm volatile("fmul %0,%1,%2" : "=f" (t01) : "f" (y01), "f" (y01));
    asm volatile("fmul %0,%1,%2" : "=f" (t11) : "f" (y01), "f" (half));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t01) : "f" (a), "f" (t01), "f" (one));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y01) : "f" (t01), "f" (t11), "f" (y01));

    asm volatile("fmul %0,%1,%2" : "=f" (t01) : "f" (y01), "f" (y01));
	  asm volatile("fmul %0,%1,%2" : "=f" (t11) : "f" (y01), "f" (half));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t01) : "f" (a), "f" (t01), "f" (one));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y01) : "f" (t01), "f" (t11), "f" (y01));

    asm volatile("fmul %0,%1,%2" : "=f" (t01) : "f" (y01), "f" (y01));
//	asm volatile("fmul %0,%1,%2" : "=f" (t11) : "f" (y01), "f" (half));
    asm volatile("fnmsub %0,%1,%2,%3" : "=f" (t01) : "f" (a), "f" (t01), "f" (one));
    asm volatile("fmadd %0,%1,%2,%3" : "=f" (y01) : "f" (t01), "f" (t11), "f" (y01));

	return y01;
#endif
}*/


#define _PI_ (3.14159265358979323841218)
#define _2PI_ (2*_PI_)
#define _4PI_ (4*_PI_)
#define __4PI_ 1/(_4PI_)
#define SQRT3 (1.732050807568877)
#define SQRT2 (1.414213562373095)
#define _SQRT2 (1/SQRT2)
#define _SQRT3 (1/SQRT3)

#endif // _mymath_h_
