#ifndef BOX_H
#define BOX_H

#include <vector>
#include "levelsetproject.h"

using namespace std;

class matrix;
class LSbox;

struct neighbor{
	LSbox* who;
	double sigma;
};

struct pointVal {
    int x,y;
    double val;
	int direction;
    pointVal(int yy, int xx, double aVal, int dir):x(xx), y(yy), val(aVal), direction(dir){}
};

//box is ambiguous, so L(evel)S(et)box...
class LSbox {
    int id;
    int xmin, xmax, ymin, ymax;
    vector<pointVal> zeros;
	double* distance;
    matrix* domain;
	bool exist;
    LSbox **IDLocal[2]; 	// local array to asign a cell id to each grid point
    
public:
	friend class matrix;
    LSbox();
    ~LSbox();
	vector<neighbor> weights;
    vector<LSbox*> neighbors;
	vector<LSbox*> neighbors_2order;
    LSbox(int aID, voro::voronoicell_neighbor& c, double *part_pos, int grid_blowup, double h);
    LSbox distancefunction(voro::voronoicell_neighbor& c, int *gridIDs, double *part_pos, int grid_blowup, double h);
    void copy_distances();
	void copy_distances_to_domain();
	void sweeping (double h, int start_i, int start_j, int direction);
	void redistancing(double h, int grid_blowup );
	void setZeros(double h, int grid_blowup, int loop);
	void sweep(pointVal zero, double h);
    int  getID();
    void setDomain(matrix* aDomain);
    void comparison(const matrix &domain_copy, int loop );
    void comparison_set_to_domain(LSbox ***ID, int resized_m, int grid_blowup);
    void add_n2o();
	void maximum(const matrix &A, const matrix &B);
    bool checkIntersect(LSbox* box2);   	
	void free_memory_distance();
    double curvature (int x, int y, double h);
    void euler_forward(double dt, double h);
	void plot_box(bool distanceplot);
	inline bool get_status() { return exist;}
	inline int get_id() { return getID(); }
};






#endif