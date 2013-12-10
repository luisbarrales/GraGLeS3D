#include "functions.h"

using namespace std;


// LU Decomposition of A=L*U with partial pivoting
//void LUDecomp(int n, domainCl& A) {
 //   for (int k = 0; k < n; k++){
 //       for (int i = k+1; i < n; i++)
 //           if (A[i][k] != 0) {
 //               A[i][k] /= A[k][k];
 //           }
  ///      for (int j=k+1;j<n;j++)
  //          for (int i=k+1;i<n;i++)
  //              if (A[i][k]!=0) 
  //                  A[i][j]-=A[i][k]*A[k][j];
  //  }
//}

// For A*x=b we get L*U*x=b and set y = U*x
// Hence to get x, we first solve L*y = b
void FSubst(int n, const domainCl& LU, const vektor& b, vektor& y) {
    for (int i=0;i<n;i++) {
        y[i]=b[i];
        for (int j=0;j<i;j++)
            y[i]-=LU[i][j]*y[j];
    }
}

// Having y we solve U*x=y to get x
void BSubst(int n, const domainCl& LU, const vektor& y, vektor& x) {
    for (int i=n-1;i>=0;i--) {
        x[i]=y[i];
        for (int j=n-1;j>i;j--)
            x[i]-=LU[i][j]*x[j];
        x[i]/=LU[i][i];
    }
}

