#ifndef MATRIX_h
#define MATRIX_h

#include "levelsetproject.h"

using namespace std;

class vektor;
class LSbox;
class particle_information;

class matrix{
    int m,n,id;
    vektor **x;
    vector<LSbox*> grains;
	
public:
    friend class vektor;
    friend class LSbox;
    friend ostream& operator<<(ostream &os, const matrix& A);
	matrix();
    matrix(int m, int n);
	matrix(int m, int n, int id);
	matrix(int m, int n, int id, double startval);
	~matrix();
    matrix(const matrix& v);
	
    vektor& operator[](int i);
    const vektor& operator[](int i) const;
	
	double& operator=(const matrix& A);
    matrix	operator+(const matrix& A);
    matrix	operator-(matrix& A);
    matrix	operator*(const matrix& A);
    vektor	operator*(const vektor& v);
    matrix	operator/(matrix& A);
	double entry(const int i, const int j);

	matrix distancefunction(voro::voronoicell_neighbor& c, int *gridIDs, double *part_pos, int grid_blowup, double h);
	matrix energy_correction(const LSbox ***&ST);
	// grain check: Erklärung anfügen!!
    bool grainCheck(double h,  int grid_blowup, vector<LSbox*>& buffer);
	bool addBox(LSbox* aBox);
    void redistancing_for_all_boxes(double h, int grid_blowup); // ruft boxweise auf
	vector<LSbox*> getBoxList();
    
	void maximum(const matrix &A, const matrix &B);
    int  minimumInPoint(std::list<matrix> distances, int m, int n, int neglect);
	void matrix_to_array(double *u);
	void array_to_matrix(double *u);
	void makeFFTPlans(double *u, fftw_complex *fftTemp, fftw_plan *fftplan1, fftw_plan *fftplan2);

	void conv_generator(double *u, fftw_complex *fftTemp, fftw_plan fftplan1, fftw_plan fftplan2, double dt);
	void convolution(double dt);
	
    void redistancing_advanced(double h, int grid_blowup, std::list<matrix> distances, double** borderSlopes, double** slopeField);
    void redistancing(double h, int grid_blowup);
	void redistancing_2(double h, int grid_blowup);

	bool discrete_convolution(const double dt, const double h, const int grid_blowup, double (*kernel)(double,int,int,int));
	void comparison(std::list<matrix> distances, int grid_blowup);
	void set_border_to_INTERIMVAL(int grid_blowup);
	void save_matrix(const char* filename ); //schreibt datei für gnuplotskript
	void clear_domain(double value);
	void euler(double dt, double h);

	int get_m() const;
	int get_n() const;
	int get_id() const;
	int get_nr_of_grains();
	void mult_with_scalar(const double d);
	

};

#endif
