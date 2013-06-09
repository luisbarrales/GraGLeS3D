#ifndef BOX_H
#define BOX_H

#include <vector>
#include "levelsetproject.h"

using namespace std;

class matrix;

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
    

public:
    LSbox();
    ~LSbox();
    vector<LSbox*> neighbors;
	vector<LSbox*> neighbors_2order;
    LSbox(int aID, voro::voronoicell_neighbor& c, double *part_pos, int grid_blowup, double h);
//     void determineNeighbors()
    LSbox distancefunction(voro::voronoicell_neighbor& c, LSbox ***&ID_mat, double *part_pos, int grid_blowup, double h);
    void copy_distances();
	void copy_distances_to_domain();
	void sweeping (double h, int start_i, int start_j, int direction);
	void redistancing(double h, int grid_blowup /*, std::list<matrix> distances, double** borderSlopes, double** slopeField*/);
	bool setZeros(double h, int grid_blowup);
	void sweep(pointVal zero, double h);
    int  getID();
    void setDomain(matrix* aDomain);
    void comparison(const matrix &domain_copy);
      void comparison_set_to_domain();
    void maximum(const matrix &A, const matrix &B);
    bool checkIntersect(LSbox* box2);   	
	void free_memory_distance();
    double curvature (int x, int y, double h);
    void euler_forward(double dt, double h);
};






#endif