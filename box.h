#ifndef BOX_H
#define BOX_H

#include <vector>
#include "levelsetproject.h"

using namespace std;

class matrix;

struct pointVal {
    int x,y;
    double val;
    
    pointVal(int xx, int yy, double aVal):x(xx), y(yy), val(aVal) {}
};

//box is ambiguous, so L(evel)S(et)box...
class LSbox {
    
    int id;
    double xmin, xmax, ymin, ymax;
    vector<pointVal> zeros;
    matrix* domain;

public:
    LSbox();
    ~LSbox();
    LSbox(int aID, voro::voronoicell_neighbor& c, double *part_pos, int grid_blowup, double h);
    LSbox distancefunction(voro::voronoicell_neighbor& c, int *ID_mat, double *part_pos, int grid_blowup, double h);
    void setZeros(double h);
    bool checkIntersect(LSbox* box2);
    
};






#endif