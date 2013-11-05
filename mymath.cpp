
#include "mymath.h"
#include <math.h>
#include <time.h>


/*#ifdef __ALTIVEC__


void fsqrt4( float *x1, float *x2, float *x3, float *x4 )
{
	typedef union
	{
		vector float v;
		float s[4];
	} vec;
	
	vec a,b;
	
	a.s[0] = *x1;
	a.s[1] = *x2;
	a.s[2] = *x3;
	a.s[3] = *x4;

	vector float one = vec_ctf(vec_splat_u32(1), 0 );
	vector float half = vec_ctf(vec_splat_u32(1), 1 );
	vector unsigned long vr1 = vec_splat_u32( -1 );
	vector float null = (vector float) vec_sl( vr1, vr1 );

	vector float est = vec_rsqrte( a.v );
	vector float estq = vec_madd( est, est, null );
	vector float korrfac = vec_madd( est, half, null );
	vector float delta = vec_nmsub( estq, a.v, one );

	vector float y0   = vec_re ( est ); 
	est = vec_madd( delta, korrfac, est );		// jetzt genaue wurzel
	vector float temp = vec_nmsub ( y0, est, one );  
	b.v   = vec_madd ( y0, temp, y0 );     

//	b.v = vec_madd( delta, korrfac, est );		// fuer inverse wurzel

	*x1 = b.s[0];
	*x2 = b.s[1];
	*x3 = b.s[2];
	*x4 = b.s[3];
}
#endif



void fsqrt3(float *ka, float *kb, float *kc)
{
#if LENDIAN
	*ka = sqrt(*ka);
	*kb = sqrt(*kb);
	*kc = sqrt(*kc);
#else

	const float k2 = -0.07021533f;
	const float k1 = 0.6224936f;
	const float k0 = 0.449045f;
	const float x2 = 0.29f;
	const float x1 = -1.05f;
	const float sqrt2 = 1.4142135623730950488f;

	long xpa = ((long*)ka)[0];
	long xpb = ((long*)kb)[0];
	long xpc = ((long*)kc)[0];
	
	long exa = xpa & 0x807fffff;	
	long exb = xpb & 0x807fffff;	
	long exc = xpc & 0x807fffff;	

	long xwa = xpa - 0x3f800000;
	long xwb = xpb - 0x3f800000;
	long xwc = xpc - 0x3f800000;

	xwa >>= 1;
	xwb >>= 1;
	xwc >>= 1;

	((long*)ka)[0] = exa | 0x3f800000;
	((long*)kb)[0] = exb | 0x3f800000;
	((long*)kc)[0] = exc | 0x3f800000;

	xwa &= 0xff800000;
	xwb &= 0xff800000;
	xwc &= 0xff800000;

	float ak = *ka;
	float bk = *kb;
	float ck = *kc;

	float xa = ak*k2 + k1;
	float xb = bk*k2 + k1;
	float xc = ck*k2 + k1;
	float ha = 1.26f;
	float hb = 1.26f;
	float hc = 1.26f;

	xa = k0 + ak*xa;
	xb = k0 + bk*xb;
	xc = k0 + ck*xc;

	ha += xa*(x2*xa+x1);
	hb += xb*(x2*xb+x1);
	hc += xc*(x2*xc+x1);

	xa += (ak-xa*xa)*ha;
	xb += (bk-xb*xb)*hb;
	xc += (ck-xc*xc)*hc;

	xwa += 0x3f800000;
	xwb += 0x3f800000;
	xwc += 0x3f800000;
	xa += (ak-xa*xa)*ha;
	xb += (bk-xb*xb)*hb;
	xc += (ck-xc*xc)*hc;

	*ka = xa;
	*kb = xb;
	*kc = xc;

	exa = ((long*)ka)[0];
	exb = ((long*)kb)[0];
	exc = ((long*)kc)[0];

	exa &= 0x807fffff;
	exb &= 0x807fffff;
	exc &= 0x807fffff;
	exa |= xwa;
	exb |= xwb;
	exc |= xwc;
	xpa &= 0x00800000;	
	xpb &= 0x00800000;	
	xpc &= 0x00800000;	

	((long*)ka)[0] = exa;
	((long*)kb)[0] = exb;
	((long*)kc)[0] = exc;

	if( !xpa )	*ka *= sqrt2;
	if( !xpb )	*kb *= sqrt2;
	if( !xpc )	*kc *= sqrt2;

#endif
}


void fsqrt2(float *ka, float *kb)
{
#if LENDIAN
	*ka = sqrt(*ka);
	*kb = sqrt(*kb);
#else
	
	const float k2 = -0.07021533f;
	const float k1 = 0.6224936f;
	const float k0 = 0.449045f;
	const float x2 = 0.29f;
	const float x1 = -1.05f;
	const float sqrt2 = 1.4142135623730950488f;

	long xpa = ((long*)ka)[0];
	long xpb = ((long*)kb)[0];
	
	long exa = xpa & 0x807fffff;	
	long exb = xpb & 0x807fffff;	

	long xwa = xpa - 0x3f800000;
	long xwb = xpb - 0x3f800000;

	xwa >>= 1;
	xwb >>= 1;

	((long*)ka)[0] = exa | 0x3f800000;
	((long*)kb)[0] = exb | 0x3f800000;

	xwa &= 0xff800000;
	xwb &= 0xff800000;

	float ak = *ka;
	float bk = *kb;

	float xa = ak*k2 + k1;
	float xb = bk*k2 + k1;
	float ha = 1.26f;
	float hb = 1.26f;

	xa = k0 + ak*xa;
	xb = k0 + bk*xb;


	ha += xa*(x2*xa+x1);
	hb += xb*(x2*xb+x1);

	xa += (ak-xa*xa)*ha;
	xb += (bk-xb*xb)*hb;
	xwa += 0x3f800000;
	xwb += 0x3f800000;
	xa += (ak-xa*xa)*ha;
	xb += (bk-xb*xb)*hb;

	*ka = xa;
	*kb = xb;

	exa = ((long*)ka)[0];
	exb = ((long*)kb)[0];

	exa &= 0x807fffff;
	exb &= 0x807fffff;
	exa |= xwa;
	exb |= xwb;
	xpa &= 0x00800000;	
	xpb &= 0x00800000;	

	((long*)ka)[0] = exa;
	((long*)kb)[0] = exb;

	if( !xpa )	*ka *= sqrt2;
	if( !xpb )	*kb *= sqrt2;

#endif
} */


/*
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
*/

/*void _fsqrt3(float *ka, float *kb, float *kc)
{
#if LENDIAN
	*ka = 1 / sqrt( *ka );
	*kb = 1 / sqrt( *kb );
	*kc = 1 / sqrt( *kc );
#else
	float a = *ka;
	float b = *kb;
	float c = *kc;
	
    float half = 0.5, one = 1.0;
    float y01;
	float t01, t11; 
    float y02; 
	float t02, t12;
    float y03; 
	float t03, t13;

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

	*ka = y01;
	*kb = y02;
	*kc = y03;
#endif
}



void dsqrt2(double *ka, double *kb)
{
#if LENDIAN
	*ka = sqrt(*ka);
	*kb = sqrt(*kb);
#else
	const double sqrt2 = 1.4142135623730950488;
	const double k0 = 0.449045;
	const double k1 = 0.6224936;
	const double k2 = -0.07021533;
	const double x2 = 0.29;
	const double x1 = -1.05;

	long xpa = ((long*)ka)[LENDIAN];
	long xpb = ((long*)kb)[LENDIAN];

	double ha = 1.26;
	double hb = 1.26;

	long exa = xpa & 0x800fffff;	
	long exb = xpb & 0x800fffff;	

	long xwa = xpa - 0x3ff00000;
	long xwb = xpb - 0x3ff00000;

	((long*)ka)[LENDIAN] = exa | 0x3ff00000;
	((long*)kb)[LENDIAN] = exb | 0x3ff00000;

	xwa >>= 1;
	xwb >>= 1;
	double ak = *ka;
	double bk = *kb;

	double xa = k1 + ak * k2;
	double xb = k1 + bk * k2;

	xa = k0 + ak*xa;
	xb = k0 + bk*xb;

	ha += xa*(x2*xa + x1);
	hb += xb*(x2*xb + x1);

	xa += (ak-xa*xa)*ha;
	xb += (bk-xb*xb)*hb;
	xa += (ak-xa*xa)*ha;
	xb += (bk-xb*xb)*hb;

	xwa &= 0xfff00000;
	xwb &= 0xfff00000;
	
	xa += (ak-xa*xa)*ha;
	xb += (bk-xb*xb)*hb;
	xa += (ak-xa*xa)*ha;
	xb += (bk-xb*xb)*hb;

	xwa += 0x3ff00000;
	xa += (ak-xa*xa)*ha;
	xwb += 0x3ff00000;
	xb += (bk-xb*xb)*hb;

	*ka = xa;
	*kb = xb;

	exa = ((long*)ka)[LENDIAN];
	exb = ((long*)kb)[LENDIAN];

	exa &= 0x800fffff;
	exb &= 0x800fffff;

	xpa &= 0x00100000;
	exa |= xwa;
	xpb &= 0x00100000;
	exb |= xwb;

	((long*)ka)[LENDIAN] = exa;
	((long*)kb)[LENDIAN] = exb;

	if( !xpa )	*ka *= sqrt2;
	if( !xpb )	*kb *= sqrt2;

#endif
}



void dsqrt3(double *ka, double *kb, double *kc)
{
#if LENDIAN
	*ka = sqrt(*ka);
	*kb = sqrt(*kb);
	*kc = sqrt(*kc);
#else

	const double sqrt2 = 1.4142135623730950488;
	const double k0 = 0.449045;
	const double k1 = 0.6224936;
	const double k2 = -0.07021533;
	const double x2 = 0.29;
	const double x1 = -1.05;

	long xpa = ((long*)ka)[LENDIAN];
	long xpb = ((long*)kb)[LENDIAN];
	long xpc = ((long*)kc)[LENDIAN];

	double ha = 1.26;
	double hb = 1.26;
	double hc = 1.26;

	long exa = xpa & 0x800fffff;	
	long exb = xpb & 0x800fffff;	
	long exc = xpc & 0x800fffff;	

	long xwa = xpa - 0x3ff00000;
	long xwb = xpb - 0x3ff00000;
	long xwc = xpc - 0x3ff00000;

	((long*)ka)[LENDIAN] = exa | 0x3ff00000;
	((long*)kb)[LENDIAN] = exb | 0x3ff00000;
	((long*)kc)[LENDIAN] = exc | 0x3ff00000;

	xwa >>= 1;
	xwb >>= 1;
	xwc >>= 1;
	double ak = *ka;
	double bk = *kb;
	double ck = *kc;

	double xa = k1 + ak * k2;
	double xb = k1 + bk * k2;
	double xc = k1 + ck * k2;

	xa = k0 + ak*xa;
	xb = k0 + bk*xb;
	xc = k0 + ck*xc;

	ha += xa*(x2*xa + x1);
	hb += xb*(x2*xb + x1);
	hc += xc*(x2*xc + x1);

	xa += (ak-xa*xa)*ha;
	xb += (bk-xb*xb)*hb;
	xc += (ck-xc*xc)*hc;
	xa += (ak-xa*xa)*ha;
	xb += (bk-xb*xb)*hb;
	xc += (ck-xc*xc)*hc;
	xwa &= 0xfff00000;
	xwb &= 0xfff00000;
	xwc &= 0xfff00000;

	xa += (ak-xa*xa)*ha;
	xb += (bk-xb*xb)*hb;
	xc += (ck-xc*xc)*hc;
	xa += (ak-xa*xa)*ha;
	xb += (bk-xb*xb)*hb;
	xc += (ck-xc*xc)*hc;

	xwa += 0x3ff00000;
	xwb += 0x3ff00000;
	xwc += 0x3ff00000;
	xa += (ak-xa*xa)*ha;
	xb += (bk-xb*xb)*hb;
	xc += (ck-xc*xc)*hc;

	*ka = xa;
	*kb = xb;
	*kc = xc;

	exa = ((long*)ka)[LENDIAN];
	exb = ((long*)kb)[LENDIAN];
	exc = ((long*)kc)[LENDIAN];

	exa &= 0x800fffff;
	exb &= 0x800fffff;
	exc &= 0x800fffff;

	exa |= xwa;
	exb |= xwb;
	exc |= xwc;

	xpa &= 0x00100000;
	xpb &= 0x00100000;
	xpc &= 0x00100000;

	((long*)ka)[LENDIAN] = exa;
	((long*)kb)[LENDIAN] = exb;
	((long*)kc)[LENDIAN] = exc;

	if( !xpa )	*ka *= sqrt2;
	if( !xpb )	*kb *= sqrt2;
	if( !xpc )	*kc *= sqrt2;

#endif
}


float Efsqrt(float k)
{
	const float k2 = -0.07021533f;
	const float k1 = 0.6224936f;
	const float k0 = 0.449045f;
	const float x2 = 0.29f;
	const float x1 = -1.05f;
	const float sqrt2 = 1.4142135623730950488f;
	long xp;

	xp = ((long*)&k)[0];
	
	long xp2 = xp & 0x807fffff;	
	long xw = xp - 0x3f800000;

	xw >>= 1;
	((long*)&k)[0] = xp2 | 0x3f800000;

	xw &= 0xff800000;

	float x = k*k2 + k1;
	
	float h = 1.26f;

	x = k0 + k*x;
	xw += 0x3f800000;
	h += x*(x2*x+x1);
	x += (k-x*x)*h;
	x += (k-x*x)*h;

	float hx = x;

	long ex = ((long*)&hx)[0];
	ex &= 0x807fffff;
	xp &= 0x00800000;	
	ex |= xw;

	((long*)&hx)[0] = ex;

	return xp ? hx : hx*sqrt2;
}



double Edsqrt(double k)
{
	const double sqrt2 = 1.4142135623730950488;
	const double k0 = 0.449045;
	const double k1 = 0.6224936;
	const double k2 = -0.07021533;
	const double x2 = 0.29;
	const double x1 = -1.05;
	double h = 1.26;
	long xp;

	xp = ((long*)&k)[LENDIAN];

	long xp2 = xp & 0x800fffff;	
	long xw = xp - 0x3ff00000;
	((long*)&k)[LENDIAN] = xp2 | 0x3ff00000;
	xw >>= 1;

	double x = k1 + k*k2;
	x = k0 + k*x;
	double tt = x2*x + x1;
	h += x*tt;
	x += (k-x*x)*h;
	x += (k-x*x)*h;
	
	xw &= 0xfff00000;
	
	x += (k-x*x)*h;
	x += (k-x*x)*h;
	xw += 0x3ff00000;
	x += (k-x*x)*h;

	double hx = x;

	long ex = ((long*)&hx)[LENDIAN];
	ex &= 0x800fffff;
	xp &= 0x00100000;
	ex |= xw;
	((long*)&hx)[LENDIAN] = ex;
	
	return xp ? hx : sqrt2*hx;
}*/


double dcos(double x)	{ return cos(x*(_PI_/180)); }
double dsin(double x)	{ return sin(x*(_PI_/180)); }

// double qroot(double x)
// {
// 	QUICKASSERT(x>=0);
// 
// 	return pow( x, 1/3.0 );
// }


double kante(double x)
{
//	double x3 = (x*x)*(x*8);
//	return 0.5 + x * pow( 1 + x3*x3 , -1/6.0 ) ;

	double x4 = (x*x)*(x*x)*16;
	return 0.5 + x * pow( 1 + x4*x4 , -1/8.0 ) ;
}


double WLS(double r)
{
	if( r>= 1.5 )	return 0;
	if( r < 0 )		return 0;
	double h;
	h  = (7.0/3.0) * log( 3.0 / (3.0+r) );
	h += (11.0/3.0) * log( 1.5 / (1.5-r) );
	h += r / (r-1.5);
	return (4.0/9.0) *r*r * exp(h);
}

/*
double wicksell(double u, double *srad, long sCount)
{
	double sum = 0;

	if( u < 0 )		u = 0;
	double u2 = SQR(u);
	long i;

	for(i=0;i<sCount;i++)
	{
		double rs = srad[i];
		
		if( rs > u )
		{
			double nen2 = SQR(rs) - u2;
			sum += 1 / sqrt( nen2 );
		}
	}
	eventLoop();
	sum *= (2 / _PI_) / sCount;
	sum = 1-sum;

	return sum;
}*/

// double misorientationCubicQxQ( double q01, double q11, double q21, double q31, double q02, double q12, double q22, double q32 )
// {
//         int i;
// 
//         Real p[4] = {q01,q11,q21,q31};
// 	Real q[4] = {q02,q12,q22,q32};
// 
// 	Real qm1[4];    //Inverse of quaternion q
// 
// 	for(i=0;i<4;i++)       //Inverting unit quaternion
//         {
// 	        qm1[i]=q[i];
//                 if( i>0 ) qm1[i]*=-1;
//         }
// 
// 	Real r[4];     //Resulting quaternion, rotation of the two previous quaternions pq-1
// 
//         multiplyQuaternions( p, qm1, r );
// 
//         //Now, we have to determine the smallest angle.
// 
// 	Real r0[6][4];    //There are 12 possible angles
// 
//         //Note: The notation r0 is due to the definition of the quaternion which lie
//         //in the fundamental zone, this vector possesses the smallest angle, in such a way
//         //that r0 is actually the scalar part of this quaternion
// 
// 	double a,b,c,d;
// 	Real rt3=sqrt(3.0);
// 
// 	a=r[0]; b=r[1]; c=r[2]; d=r[3];
// 
// 	Real fac=0.70710678;
// 
// 	r0[0][0]=(r[0]+r[1])*fac; r0[0][1]=(r[0]-r[1])*fac; r0[0][2]=(r[2]+r[3])*fac; r0[0][3]=(r[2]-r[3])*fac;
// 	r0[1][0]=(r[0]+r[2])*fac; r0[1][1]=(r[0]-r[2])*fac; r0[1][2]=(r[1]+r[3])*fac; r0[1][3]=(r[1]-r[3])*fac;
// 	r0[2][0]=(r[0]+r[3])*fac; r0[2][1]=(r[0]-r[3])*fac; r0[2][2]=(r[1]+r[2])*fac; r0[2][3]=(r[1]-r[2])*fac;
// 	r0[3][0]=(r[0]+r[1]+r[2]+r[3])*0.5; r0[3][1]=(r[0]+r[1]-r[2]-r[3])*0.5; r0[3][2]=(r[0]-r[1]+r[2]-r[3])*0.5; r0[3][3]=(r[0]-r[1]-r[2]+r[3])*0.5;
// 	r0[4][0]=(r[0]+r[1]+r[2]-r[3])*0.5; r0[4][1]=(r[0]+r[1]-r[2]+r[3])*0.5; r0[4][2]=(r[0]-r[1]+r[2]+r[3])*0.5; r0[4][3]=(r[0]-r[1]-r[2]-r[3])*0.5;
// 	r0[5][0]=r[0];r0[5][1]=r[1];r0[5][2]=r[2];r0[5][3]=r[3];
// 
// 
// 	Real omega=0.0;
// 
// 	for(i=0;i<6;i++)
// 		for( int j=0;j<4;j++ )
// 			if( fabs(r0[i][j]) > omega )
// 				omega=fabs(r0[i][j]);
// 
// 	QUICKASSERT( omega < 1.01 );
// 
// 	if( omega > 1.0 )
// 		omega = (Real) (int) omega;
// 
// 	omega=2*acos(omega);
// 	QUICKASSERT( omega <= 1.099 );
// 	return omega;
// }

// void randomMisorientationAxisConsidered(  double * qref, double * qr, double maxDev  )
// {
//         Real theta = cos( 0.5 * _PI_ );
// 
//         double q[4] = {0};
//         maxDev *= _PI_/180;
// 
//         double dev = cos(0.5 * maxDev);
// 
//         double refNorm = sqrt( SQR(qref[0]) + SQR(qref[1]) + SQR(qref[2]) + SQR(qref[3]) );
// 
//         qref[0] /= refNorm;
//         qref[1] /= refNorm;
//         qref[2] /= refNorm;
//         qref[3] /= refNorm;
// 
//         while( theta < dev  )
//         {
//                 double s = applic.myRandom.parkMiller();
//                 double sigma1 = sqrt(1-s);
//                 double sigma2 = sqrt(s);
//                 double theta1 = 2*_PI_*applic.myRandom.parkMiller();
//                 double theta2 = 2*_PI_*applic.myRandom.parkMiller();
// 
//                 q[0]=fabs(sigma2*cos(theta2));
//                 q[1]=fabs(sigma1*sin(theta1));
//                 q[2]=fabs(sigma1*cos(theta1));
//                 q[3]=fabs(sigma2*sin(theta2));
// 
//                 double norm = sqrt(SQR(q[0])+SQR(q[1])+SQR(q[2])+SQR(q[3]));
// 
//                 q[0] /= norm;
//                 q[1] /= norm;
//                 q[2] /= norm;
//                 q[3] /= norm;
// 
//                 bubbleSort( q, 4 );
// 
//                 theta = q[3]*qref[0] + q[2]*qref[1] + q[1]*qref[2] + q[0]*qref[3];
//         }
//         qr[0]=q[3];
//         qr[1]=q[2];
//         qr[2]=q[1];
//         qr[3]=q[0];
// 
// }

double angleBetweenQuaternions( double * q, double * p )
{
        double _qNorm = 1 / sqrt(SQR(q[1])+SQR(q[2])+SQR(q[3])+SQR(q[0]));
        double _pNorm = 1 / sqrt(SQR(p[1])+SQR(p[2])+SQR(p[3])+SQR(p[0]));
        return acos( _qNorm * _pNorm * ( q[0]*p[0] + q[1]*p[1] + q[2]*p[2] + q[3]*p[3] ) );
}

// void randomMisorientation( double theta, double* qr  )
// {
//         double q[4]={0,0,0,0};
//         theta*=(_PI_/180);
//         double qcrit = cos(0.5*theta);
// 
//         while( q[3]<qcrit )
//         {
//                 double s = applic.myRandom.parkMiller();
//                 double sigma1 = sqrt(1-s);
//                 double sigma2 = sqrt(s);
//                 double theta1 = 2*_PI_*applic.myRandom.parkMiller();
//                 double theta2 = 2*_PI_*applic.myRandom.parkMiller();
//                 q[0]=fabs(sigma2*cos(theta2));
//                 q[1]=fabs(sigma1*sin(theta1));
//                 q[2]=fabs(sigma1*cos(theta1));
//                 q[3]=fabs(sigma2*sin(theta2));
//                 bubbleSort( q,4 );
//         }
//         qr[0]=q[3];
//         qr[1]=q[2];
//         qr[2]=q[1];
//         qr[3]=q[0];
// }

void multiplyQuaternions( double *q, double* p, double* r )
{
        r[0]=q[0]*p[0]-q[1]*p[1]-q[2]*p[2]-q[3]*p[3];
        r[1]=q[1]*p[0]+q[0]*p[1]-q[3]*p[2]+q[2]*p[3];
        r[2]=q[2]*p[0]+q[3]*p[1]+q[0]*p[2]-q[1]*p[3];
        r[3]=q[3]*p[0]-q[2]*p[1]+q[1]*p[2]+q[0]*p[3];
}

// void misorientationQuaternionCubic( double* p, double* q, double* quat  )
// {
// 
// 	Real qm1[4];
// 	int i;
// 
// 	for(i=0;i<4;i++)                 //inverting unit quaternion q
//         {				//Copy quaternion; not really necessary
// 		qm1[i]=q[i];
//                 if( i>0 ) qm1[i]*=-1;
//         }
// 
// 	Real r1[4];
// 
//         multiplyQuaternions( p, qm1, r1 );
// 
// /***********From this point on, the calculation differs for different lattices due to symmetry***********/
// 
// 
// 	Real sqrt2=1/sqrt(2.0);
// 	Real rot[6][4];
// 	Real a, b, c ,d;
// 
// 	//The six fundamental quaternions describing the misorientation
// 
// 	a=r1[0]; b=r1[1]; c=r1[2]; d=r1[3];
// 
// 	rot[0][0] = a; rot[0][1] = b; rot[0][2] = c; rot[0][3] = d;
// 
// 	rot[1][0] = sqrt2 * ( a + b ); rot[1][1] = sqrt2 * ( a - b ); rot[1][2] = sqrt2 * ( c + d ); rot[1][3] = sqrt2 * ( c - d );
// 
// 	rot[2][0] = sqrt2 * ( a + c ); rot[2][1] = sqrt2 * ( a - c ); rot[2][2] = sqrt2 * ( b + d ); rot[2][3] = sqrt2 * ( b - d );
// 
// 	rot[3][0] = sqrt2 * ( a + d ); rot[3][1] = sqrt2 * ( a - d ); rot[3][2] = sqrt2 * ( b + c ); rot[3][3] = sqrt2 * ( b - c );
// 
// 	rot[4][0] = 0.5 * ( a + b + c + d ); rot[4][1] = 0.5 * ( a + b - c - d ); rot[4][2] = 0.5 * ( a - b + c - d ); rot[4][3] = 0.5 * ( a - b - c + d );
// 
// 	rot[5][0] = 0.5 * ( a + b + c - d ); rot[5][1] = 0.5 * ( a + b - c + d ); rot[5][2] = 0.5 * ( a - b + c + d ); rot[5][3]= 0.5 * ( a - b - c - d );
// 
// 	Real rq[4];
// 
// 	int mi;
// 	Real max=0.0;
// 	int j=0;
// 
// 	for( i=0;i<6;i++ )						//Determing the quaternion with the maximal component and the component itself
// 		for( j=0;j<4;j++ )
// 		{
// 			if( fabs(rot[i][j]) > max )
// 			{
// 				max=fabs(rot[i][j]);
// 				mi=i;
// 				//mj=j;
// 			}
// 		}
// 
// 	rq[0] = fabs( rot[mi][0] );					//Desorientation requires all components positive
// 	rq[1] = fabs( rot[mi][1] );
// 	rq[2] = fabs( rot[mi][2] );
// 	rq[3] = fabs( rot[mi][3] );
// 
// 	bubbleSort( rq,4 );						//Sorting into ascending order, because a desorientation in the SST
// 									//requires a quaternion with q0>q1>q2>q3 which represents a minimal 
// 	quat[0] = rq[3];						//rotation angle and an axis fulfilling h>k>l
// 	quat[1] = rq[2];
// 	quat[2] = rq[1];
// 	quat[3] = rq[0];
// 
// }

double misorientationCubic( double pa1, double Pa, double pa2, double pb1, double Pb, double pb2 )
{
	int i;

        /*pa1*=_PI_/180; Pa*=_PI_/180; pa2*=_PI_/180;
        pb1*=_PI_/180; Pb*=_PI_/180; pb2*=_PI_/180;

        Real p1 = pa1 - 0.5 * _PI_;             Real p12 = pb1 - 0.5 * _PI_;
	Real t = Pa;                            Real t2 = Pb;
	Real p2 = pa2 - 1.5 * _PI_;             Real p22 = pb2 - 1.5 * _PI_;

	Real co1 = cos(t/2); Real co2 = cos(t2/2);
	Real s1 = sin(t/2); Real s2 = sin(t2/2);*/

        Real oria[3] = { pa1, Pa, pa2 };
        Real orib[3] = { pb1, Pb, pb2 };

        Real p[4] /*= {co1*cos((p1+p2)/2),-s1*sin((p1-p2)/2),s1*cos((p1-p2)/2),co1*sin((p1+p2)/2)}*/;
	Real q[4] /*= {co2*cos((p12+p22)/2),-s2*sin((p12-p22)/2),s2*cos((p12-p22)/2),co2*sin((p12+p22)/2)}*/;

        euler2quaternion( oria, p );
        euler2quaternion( orib, q );

	Real qm1[4];    //Inverse of quaternion q

	for(i=0;i<4;i++)               //Inverting unit quaternion
        {
	        qm1[i]=q[i];
                if( i>0 ) qm1[i]*=-1;
        }

	Real r[4];     //Resulting quaternion, rotation of the two previous quaternions pq-1

        multiplyQuaternions( p, qm1, r );

        //Now, we have to determine the smallest angle.

	Real r0[6][4];    //There are 12 possible angles

        //Note: The notation r0 is due to the definition of the quaternion which lie
        //in the fundamental zone, this vector possesses the smallest angle, in such a way
        //that r0 is actually the scalar part of this quaternion

	double a,b,c,d;
	Real rt3=sqrt(3.0);

	a=r[0]; b=r[1]; c=r[2]; d=r[3];

	Real fac=0.70710678;

	r0[0][0]=(r[0]+r[1])*fac; r0[0][1]=(r[0]-r[1])*fac; r0[0][2]=(r[2]+r[3])*fac; r0[0][3]=(r[2]-r[3])*fac;
	r0[1][0]=(r[0]+r[2])*fac; r0[1][1]=(r[0]-r[2])*fac; r0[1][2]=(r[1]+r[3])*fac; r0[1][3]=(r[1]-r[3])*fac;
	r0[2][0]=(r[0]+r[3])*fac; r0[2][1]=(r[0]-r[3])*fac; r0[2][2]=(r[1]+r[2])*fac; r0[2][3]=(r[1]-r[2])*fac;
	r0[3][0]=(r[0]+r[1]+r[2]+r[3])*0.5; r0[3][1]=(r[0]+r[1]-r[2]-r[3])*0.5; r0[3][2]=(r[0]-r[1]+r[2]-r[3])*0.5; r0[3][3]=(r[0]-r[1]-r[2]+r[3])*0.5;
	r0[4][0]=(r[0]+r[1]+r[2]-r[3])*0.5; r0[4][1]=(r[0]+r[1]-r[2]+r[3])*0.5; r0[4][2]=(r[0]-r[1]+r[2]+r[3])*0.5; r0[4][3]=(r[0]-r[1]-r[2]-r[3])*0.5;
	r0[5][0]=r[0];r0[5][1]=r[1];r0[5][2]=r[2];r0[5][3]=r[3];

	Real omega=0.0;

	for(i=0;i<6;i++)
		for( int j=0;j<4;j++ )
			if( fabs(r0[i][j]) > omega )
				omega=fabs(r0[i][j]);

// 	QUICKASSERT( omega < 1.01 );

	if( omega > 1.0 )
		omega = (Real) (int) omega;

	omega=2*acos(omega);
// 	QUICKASSERT( omega <= 1.099 );
	return omega;
}

// void newOrientationFromReferenceFixedAngularCone(double * oriOri, double maxDev,double angle, double u, double v, double w,double * newOri)
// {
//         Real qr[4];
//         Real ori[4], qideal[4];
// 
//        /* Real p1=_PI_*(oriOri[0]/180.0);
//         Real teta=_PI_*(oriOri[1]/180.0);
//         Real p2=_PI_*(oriOri[2]/180.0);
// 
//         Real co1=cos(teta/2);
// 	Real s1=sin(teta/2);*/
// 
//         angle *= _PI_/180;
//         Real _norm = 1.0 / sqrt(SQR(u)+SQR(v)+SQR(w));
// 
//         //Real qideal[4]={co1*cos((p1+p2)/2),-s1*sin((p1-p2)/2),s1*cos((p1-p2)/2),co1*sin((p1+p2)/2)};
// 
//         euler2quaternion( oriOri, qideal );
//         Real qref[4] = { cos( 0.5* angle ), u * _norm * sin( 0.5 * angle ), v * _norm * sin( 0.5 * angle ), w * _norm * sin( 0.5 * angle ) };
// 
//         randomMisorientationAxisConsidered(  qref, qr, maxDev  );
//         multiplyQuaternions( qr,qideal,ori );
// 
//         Real euler[3];
//         quaternion2Euler( ori, euler );
// 
//         newOri[0] = euler[0] * 180 / _PI_;
//         newOri[1] = euler[1] * 180 / _PI_;
//         newOri[2] = euler[2] * 180 / _PI_;
//               //Real mis = misorientationCubicQxQ(qideal[0],qideal[1],qideal[2],qideal[3],newOri[0],newOri[1],newOri[2],newOri[3]);
// 
// }

void rotateOrientation( double *oriOri, double angle, double u, double v, double w, double *newOri )
{
        Real qr[4];
        Real ori[4], qideal[4];

        angle *= _PI_/180;
        Real _norm = 1.0 / sqrt(SQR(u)+SQR(v)+SQR(w));

        /*Real p1=_PI_*(oriOri[0]/180.0);
        Real teta=_PI_*(oriOri[1]/180.0);
        Real p2=_PI_*(oriOri[2]/180.0);

        Real co1=cos(teta/2);
	Real s1=sin(teta/2);*/

        //Real qideal[4]={co1*cos((p1+p2)/2),-s1*sin((p1-p2)/2),s1*cos((p1-p2)/2),co1*sin((p1+p2)/2)};

        euler2quaternion( oriOri, qideal );

        Real qref[4] = { cos( 0.5* angle ), u * _norm * sin( 0.5 * angle ), v * _norm * sin( 0.5 * angle ), w * _norm * sin( 0.5 * angle ) };

        multiplyQuaternions( qref,qideal,ori );

        Real euler[3] = {0};

        quaternion2Euler( ori, euler );

        newOri[0] = euler[0] * 180 / _PI_;
        newOri[1] = euler[1] * 180 / _PI_;
        newOri[2] = euler[2] * 180 / _PI_;
}

void euler2quaternion( double * euler, double * q )
{
        double p1 = euler[0]*_PI_/180;
        double t  = euler[1]*_PI_/180;
        double p2 = euler[2]*_PI_/180;

        Real co1 = cos(t/2);
	Real s1 = sin(t/2);

        Real p[4] = {co1*cos((p1+p2)/2),s1*cos((p1-p2)/2),s1*sin((p1-p2)/2),co1*sin((p1+p2)/2)};
        q[0] = p[0];
        q[1] = p[1];
        q[2] = p[2];
        q[3] = p[3];
}

void quaternion2Euler( double * quat, double * euler )
{

        Real q0 = quat[0];
        Real q1 = quat[1];
        Real q2 = quat[2];
        Real q3 = quat[3];
        Real PHI, sP, phi1, phi2;

        Real cosPHI = SQR(q3)-SQR(q2)-SQR(q1)+SQR(q0);

        Real y0 = 2*q1*q3-2*q0*q2;
        Real x0 = 2*q2*q3+2*q0*q1;
        Real y1 = 2*q1*q3+2*q0*q2;
        Real x1 = -2*q2*q3+2*q0*q1;

        if( cosPHI > 1. ) cosPHI = 1.;

        if( SQR( 1. - cosPHI ) <= 1e-20 )
                PHI = 0.;
        else
                PHI = acos( cosPHI );

        sP = sin(PHI);

        if( sP != 0 )
        {
                phi2 = atan2( y0 / sP, x0 / sP );
                phi1 = atan2( y1 / sP, x1 / sP );
        }else
        {
                phi1 = atan2( (2*q1*q2+2*q0*q3),SQR(q0)+SQR(q1)-SQR(q2)-SQR(q3) );
                phi2 = 0.;
        }

        euler[0]=phi1;
        euler[1]=PHI;
        euler[2]=phi2;
}

// void newOrientationFromReference( double *oriOri, double deviation, double *newOri)
// {
//         Real qr[4];
//         Real ori[4], qideal[4];
// 
//         /*Real p1=_PI_*(oriOri[0]/180.0);
//         Real teta=_PI_*(oriOri[1]/180.0);
//         Real p2=_PI_*(oriOri[2]/180.0);*/
// 
//         /*Real co1=cos(teta/2);
// 	Real s1=sin(teta/2);
// 
//         Real qideal[4]={co1*cos((p1+p2)/2),-s1*sin((p1-p2)/2),s1*cos((p1-p2)/2),co1*sin((p1+p2)/2)};*/
// 
//         euler2quaternion( oriOri, qideal );
// 
//         randomMisorientation( deviation, qr  );
//         multiplyQuaternions( qr,qideal,ori );
// 
//         Real euler[3];
// 
//         quaternion2Euler( ori, euler );
// 
//         newOri[0] = euler[0] * 180 / _PI_;
//         newOri[1] = euler[1] * 180 / _PI_;
//         newOri[2] = euler[2] * 180 / _PI_;
// }

char inequalitieFulfil( double a, double b, double c, double d )
{
	if( b >= c && c >= d && d >=0 )
	{
		if( b<=(sqrt(2)-1)*a )
		{
			if( b+c+d <= a )
				return 3;
			return 2;
		}
		return 1;
	}
	return 0;
}

double areaPolygon( int *poly, int counter )
{
        double sum=0;

        for( int i=0;i<(2*counter-2);i+=2 )
        {
                sum += (poly[i]+poly[i+2])*(poly[i+3]-poly[i+1]);
        }
        if( sum<0 ) sum*=-1;
        return 0.5*sum;
}

int pnpoly(int nvert, double *vertx, double *verty, double testx, double testy)
{
        int i, j, c = 0;
        for (i = 0, j = nvert-1; i < nvert; j = i++)
        {
                if ( ((verty[i]>testy) != (verty[j]>testy)) && (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
                        c = !c;
        }
        return c;
}

void bubbleSort ( Real arr [ ], int size ) // Sort components of a quaternion
 { 
    int last = size - 2; 
    int isChanged = 1; 

    while ( last >= 0 && isChanged ) 
    { 
            isChanged = 0; 
            for ( int k = 0; k <= last; k++ ) 
                if ( arr[k] > arr[k+1] ) 
                { 
                    swap ( arr[k], arr[k+1] ); 
                    isChanged = 1; 
                } 
            last--; 
    } 
 } 
 
void swap ( Real& x, Real& y ) //Required for the sorting
{
    Real temp;
    temp = x;
    x = y;
    y = temp;
}



