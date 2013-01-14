#ifndef MATRIX_h
#define MATRIX_h

#include "levelsetproject.h"

using namespace std;

class vektor;
class particle_information;

class matrix{
    int m,n,id;
    vektor **x;
	
public:
    friend class vektor;
    friend ostream& operator<<(ostream &os, const matrix& A);
	
    matrix(int m, int n);
	matrix(int m, int n, int id);
	~matrix();
    matrix(const matrix& v);
	
    vektor& operator[](int i);
    const vektor& operator[](int i) const;
	
    double& operator=(const matrix& A);
    matrix operator+(const matrix& A);
    matrix operator-(const matrix& A);
    matrix operator*(const matrix& A);
    vektor operator*(const vektor& v);
    matrix operator/(matrix& A);
	double entry(const int i, const int j);

	matrix distancefunction(voro::voronoicell_neighbor& c, int *ID_mat, double *part_pos, int grid_blowup, double h);

	
	void maximum(const matrix &A, const matrix &B);
    int minimumInPoint(std::list<matrix> distances, int m, int n, int neglect);
	void matrix_to_array(double *u);
	void array_to_matrix(double *u);
	void makeFFTPlans(double *u, fftw_complex *fftTemp, fftw_plan *fftplan1, fftw_plan *fftplan2);

	void conv_generator(double *u, fftw_complex *fftTemp, fftw_plan fftplan1, fftw_plan fftplan2, double dt);
	void convolution(double dt);
	
    void redistancing(double h, int grid_blowup, std::list<matrix> distances, double** borderSlopes, double** slopeField);
    void fastSweepingStepX(int mPos, int nPos, double h, char sign );
    void fastSweepingStepY(int mPos, int nPos, double h, char sign );
    
	void five_point_formula(const double dt, const double dx);
	bool discrete_convolution(const double dt, const double h, const int grid_blowup, double (*kernel)(double,int,int,int));
	bool comparison(std::list<matrix> distances, int grid_blowup);
	void save_matrix(const char* filename ); //schreibt datei für gnuplotskript
	

	int get_m() const;
    int get_n() const;
	int get_id() const;
	
	template <class T>
	void mult_with_scalar(const T d){
	  for (int i = 0; i < m; i++) 
		for (int j = 0; j < n; j++) 
		  (*x[i])[j] *= d;
	}
    
};

#endif
