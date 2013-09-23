#ifndef domainCl_h
#define domainCl_h

#include "levelsetproject.h"

using namespace std;

// class vektor;
class LSbox;
class weightmap;

class domainCl{
    int m,n,id;
    double **x;
    double *val;
    fftw_complex *fftTemp;
    vector<LSbox*> grains;
	grainhdl* owner;
	
	
public:
    friend class vektor;
    friend class LSbox;
    friend ostream& operator<<(ostream &os, const domainCl& A);
	domainCl();
    domainCl(int m, int n);
	domainCl(int m, int n, int id);
	domainCl(int m, int n, int id, double startval);
	domainCl(int m, int n, int id, double startval, grainhdl* owner);
	~domainCl();
    domainCl(const domainCl& v);
	
    double* operator[](int i);
    const double* operator[](int i) const;
//     double& operator()(int x,int y) ;
//     const double& operator()(int x,int y) const ;
  domainCl operator-(domainCl& A);

	double& operator=(const domainCl& A);
//        domainCl	operator*(const domainCl& A);
//     vektor	operator*(const vektor& v);
   	double entry(const int i, const int j);

	domainCl distancefunction(voro::voronoicell_neighbor& c, int *gridIDs, double *part_pos, int grid_blowup, double h);
	domainCl energy_correction(const LSbox ***&ST);
	// grain check: Erklärung anfügen!!
    bool grainCheck(double h,  int grid_blowup, vector<LSbox*>& buffer, int loop);
	bool addBox(LSbox* aBox);
    void redistancing_for_all_boxes(double h, int grid_blowup); // ruft boxweise auf
	vector<LSbox*> getBoxList();
    
	void maximum(const domainCl &A, const domainCl &B);
    int  minimumInPoint(std::list<domainCl> distances, int m, int n, int neglect);
	void makeFFTPlans(double *u, fftw_complex *fftTemp, fftw_plan *fftplan1, fftw_plan *fftplan2);

	void conv_generator(double *u, fftw_complex *fftTemp, fftw_plan fftplan1, fftw_plan fftplan2, double dt);
	void convolution(double dt, double *ST, LSbox ***ID, domainCl &ref, LSbox* zeroBox, int grid_blowup, weightmap& my_weights);
	
    void redistancing_advanced(double h, int grid_blowup, std::list<domainCl> distances, double** borderSlopes, double** slopeField);
    void redistancing(double h, int grid_blowup);
	void redistancing_2(double h, int grid_blowup);

	void comparison(std::list<domainCl> distances, int grid_blowup);
	void set_border_to_INTERIMVAL(int grid_blowup);
	void save_domainCl(const char* filename ); //schreibt datei für gnuplotskript
	void clear_domain(double value);
	void euler(double dt, double h);

	int get_m() const;
	int get_n() const;
	int get_id() const;
	int get_nr_of_grains();
	void mult_with_scalar(const double d);
};

#endif
