#ifndef FUNCTIONS_h
#define FUNCTIONS_h

#include "ggLS.h"
class vektor;

// LU Decomposition of A=L*U with partial pivoting
//void LUDecomp(int n, domainCl& A);

// For A*x=b we get L*U*x=b and set y = U*x
// Hence to get x, we first solve L*y = b
void FSubst(int n, const domainCl& LU, const vektor& b, vektor& y);

// Having y we solve U*x=y to get x
void BSubst(int n, const domainCl& LU, const vektor& y, vektor& x);

#endif
