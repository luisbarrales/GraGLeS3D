#include<iostream>
#include<cmath>
#include"matrix.h"
#include"vector.h"
#include"outOfBoundsException.h"
#include"functions.h"

using namespace std;


int main(int c, char*v[]) {
    
    const int n=3;
    matrix A(n,n); 
    A[0][0]=3;A[0][1]=5; A[0][2]=-8;
    A[1][0]=1;A[1][1]=2; A[1][2]=-1;
    A[2][0]=0;A[2][1]=-2; A[2][2]=2;
    
    cout << "Initial Matrix A:" << endl;
    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++)
            try {
                cout << A[i][j] << " ";
            }
        catch (outOfBoundsException e) {
            cout << "invalid index " << e.what() << endl;
        }
        cout << endl;
    }
    
    LUDecomp(n,A);  
    
    cout << "LU Version of Matrix A:" << endl;
    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++)
            cout << A[i][j] << " ";
        cout << endl;
    }
    
    vector b(n);
    b[0]=3; b[1]=1; b[2]=2;	
    cout << "Vector b:"	<< endl;
    for (int i=0;i<n;i++)
        cout << b[i] << " ";
    cout << endl;
    
    vector x(n);
    vector y(n);
    FSubst(n, A, b, y);
    cout << "Vector y gained by Solving L*y=b : "	<< endl;
    for (int i=0;i<n;i++)
        cout << y[i] << " ";
    cout << endl;
    
    BSubst(n, A, y, x);
    cout << "Vector x gained by Solving U*x=y : "	<< endl;
    for (int i=0;i<n;i++)
        cout << x[i] << " ";
    cout << endl;
    
    // Bsp Dereferenzierungsoperator
    cout << "Bsp Dereferenzierungsoperator" << endl << "A[0][2] = " << A[0][2] << endl;
    
    // Bsp Addition/ Subtraktion
    cout << "Vektoraddition/-subtraktion" << endl;
    vector vec1(n);
    cout << "v = ";
    for (int i = 0; i < n; i++) {
        vec1[i] = i;
        cout << vec1[i] << " ";
    }
    cout << endl;
    vector vec2(n);
    cout << "w = ";
    for (int i = 0; i < n; i++) {
        vec2[i] = i+1;
        cout << vec2[i] << " ";
    }
    cout << endl;
    vector erg1(n);
    erg1 = vec1 + vec2;
    vector erg2(n);
    erg2 = vec1 - vec2;
    cout << "v + w = ";
    for (int i = 0; i<n; i++) cout << erg1[i] << " ";
    cout << endl;
    cout << "v - w = ";
    for (int i = 0; i<n; i++) cout << erg2[i] << " ";
    cout << endl;
    
    // Bsp Multiplikation
    cout << "Skalarprodukt" << endl << "v * w = " << vec1*vec2 << endl;
    
    cout << "Matrixmultiplikation" << endl << "A = " << endl;
    matrix A1(n,n);
    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++){
            A1[i][j] = (i+1)*j;
            cout << A1[i][j] << " ";
        }
        cout << endl;
    }
    
    matrix B(n,n);
    cout << "B = " << endl;
    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++){
            B[i][j] = i*j + 1;
            cout << B[i][j] << " ";
        }
        cout << endl;
    }
    cout << "C = A * B" << endl << "C = " << endl;
    matrix ERG(n,n);
    ERG = A1*B;
    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++) cout << ERG[i][j] << " ";
        cout << endl;
    }
    cout << endl;
    
    
    // Bsp Matrix * Vektor
    
    cout << "Matrix * Vektor" << endl << "A * v =" << endl;
    vector erg3(n);
    erg3 = A1*vec1;
    for (int i=0;i<n;i++) {
        cout << erg3[i] << endl;
    }
    cout << endl;
    
    // Bsp Vektor * Matrix
    
    cout << "Vektor * Matrix" << endl << "v * A =" << endl;
    vector erg4(n);
    erg4 = vec1*A1;
    for (int i=0;i<n;i++) {
        cout << erg4[i] << endl;
    }
    cout << endl;
    
    
    // Bsp Division Vektor/ Matrix
    
    cout << "Vektor / Matrix" << endl << "v = ";
    for (int i = 0; i < n; i++) cout << vec1[i] << " ";
    cout << "Matrix A: " << endl;
    matrix A2(n,n);
    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++){
            A2[i][j] = i+j+1;
            cout << A2[i][j] << " ";
        }
        cout << endl;
    }
    cout << "v / A =" << endl;
    vector erg5(n);
    erg5 = vec1/A;
    for (int i=0;i<n;i++) {
        cout << erg5[i] << endl;
    }
    cout << endl;
    
    // Matrix / Matrix
    cout << "Matrix / Matrix" << endl << "B / A =" << endl;
    matrix ERG2(n,n);
    ERG2 = B/A;
    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++) cout << ERG2[i][j] << " ";
        cout << endl;
    }
    cout << endl;
    
    
    return 0;
}