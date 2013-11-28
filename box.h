#ifndef BOX_H
#define BOX_H

#include <vector>
#include "ggLS.h"
// #include "contour.h"
using namespace std;

class domainCl;
class LSbox;
class grainhdl;

struct neighbor{
	LSbox* who;
	double sigma;
};
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
    vector<pointVal> zeros;
	double* distance;
    domainCl* domain;
	bool exist;
    LSbox **IDLocal[2]; 	// local array to asign a cell id to each grid point
    int nvertices;
	double phi1;
	double PHI;
	double phi2;
	double volume;
	double energy;
	grainhdl* handler;
		    
public:
	friend class domainCl;
	friend class grainhdl;
    LSbox();
    ~LSbox();
	vector<neighbor> weights;
    vector<LSbox*> neighbors;
	vector<LSbox*> neighbors_2order;
	LSbox(int id, int xmin, int xmax, int ymin, int ymax, double phi1, double PHI, double phi2);
    LSbox(int aID, voro::voronoicell_neighbor& c, double *part_pos, int grid_blowup, double h, grainhdl* owner);
	LSbox(int id, int nvertex, double* vertices, double phi1, double PHI, double phi2, int grid_blowup, double h, grainhdl* owner);
	LSbox distancefunction(int nvertex, double* vertices, int grid_blowup, double h);
    LSbox distancefunction(voro::voronoicell_neighbor& c, int *gridIDs, double *part_pos, int grid_blowup, double h);
    void copy_distances();
	void copy_distances_to_domain();
	void sweeping (double h, int start_i, int start_j, int direction);
	void redist_box(double h, int grid_blowup );
	void setZeros(double h, int grid_blowup, int loop);
	void sweep(pointVal zero, double h);
    int  getID();
    void setDomain(domainCl* aDomain);
    void comparison(const domainCl &domain_copy, int loop );
    void comparison_set_to_domain(LSbox ***ID, int grid_blowup);
    void add_n2o();
	void maximum(const domainCl &A, const domainCl &B);
    bool checkIntersect(LSbox* box2);   	
	void free_memory_distance();
    double curvature (int x, int y, double h);
    void euler_forward(double dt, double h);
	void plot_box(bool distanceplot);
	double mis_ori(LSbox* grain_2);
	void checkIntersect_zero_grain(domainCl* temp);
		
	inline bool get_status() { return exist;}
	inline int get_id() { return id; }
	inline double get_vol() {return volume ;}
	inline double get_phi1() {return phi1;}
	inline double get_PHI() {return PHI;}
	inline double get_phi2() {return phi2;}
};






#endif