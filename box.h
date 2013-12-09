#ifndef BOX_H
#define BOX_H

#include <vector>
#include "ggLS.h"
// #include "contour.h"
using namespace std;

class LSbox;
class grainhdl;
class weightmap;


struct SPoint
{
   SPoint(double x,double y){this->x=x;this->y=y;}
   SPoint(){}
   double x,y;
};


struct pointVal {
    int x,y;
    double val;
	int direction;
    pointVal(int yy, int xx, double aVal, int dir):x(xx), y(yy), val(aVal), direction(dir){}
};



//box is ambiguous, so L(evel)S(et)box...
class LSbox {
    unsigned int id;
    int xmin, xmax, ymin, ymax;
    int old_xmin, old_xmax, old_ymin, old_ymax;
    vector<double> distance_new, distance_current;
    domainCl* domain;

	
	double* distance_2neighbor;
	bool exist;
	vector<vector<LSbox*>> IDLocal;
	weightmap* local_weights;
//     LSbox **IDLocal[2]; 	// local array to asign a cell id to each grid point
    int nvertices;
	double phi1;
	double PHI;
	double phi2;
	double volume;
	double energy;
	grainhdl* handler;
	
public:
	friend class grainhdl;
    LSbox();
    ~LSbox();
	
    vector<LSbox*> neighbors;
	vector<LSbox*> neighbors_old;
	vector<LSbox*> neighbors_2order;
	LSbox(int id, int xmin, int xmax, int ymin, int ymax, double phi1, double PHI, double phi2, grainhdl* owner);
    LSbox(int aID, voro::voronoicell_neighbor& c, double *part_pos, grainhdl* owner);
	LSbox(int id, int nvertex, double* vertices, double phi1, double PHI, double phi2, grainhdl* owner);
	void distancefunction(int nvertex, double* vertices);
    void distancefunction(voro::voronoicell_neighbor& c, double *part_pos);
    void redist_box();
	void find_contour();
	int  getID();
    void setDomain(domainCl* aDomain);
    void comparison();
    void set_comparison(LSbox ***ID, int grid_blowup);
    void add_n2o();
	void maximum(const domainCl &A, const domainCl &B);
    bool checkIntersect(LSbox* box2);   	
	void free_memory_distance();

	void convolution();

	void plot_box(bool distanceplot);
	double mis_ori(LSbox* grain_2);
	void checkIntersect_zero_grain();
	void shape_distance();
	
	void makeFFTPlans(double *in, double* out,fftw_complex *fftTemp, fftw_plan *fftplan1, fftw_plan *fftplan2);
	void conv_generator(fftw_complex *fftTemp, fftw_plan fftplan1, fftw_plan fftplan2);

		
	inline bool get_status() { return exist;}
	inline int get_id() { return id; }
	inline double get_vol() {return volume ;}
	inline double get_phi1() {return phi1;}
	inline double get_PHI() {return PHI;}
	inline double get_phi2() {return phi2;}
};






#endif