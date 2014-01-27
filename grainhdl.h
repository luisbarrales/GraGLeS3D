#ifndef GRAINHDL_h
#define GRAINHDL_h

#include "ggLS.h"
#include <omp.h>

#define xsect(p1,p2) (h[p2]*xh[p1]-h[p1]*xh[p2])/(h[p2]-h[p1])
#define ysect(p1,p2) (h[p2]*yh[p1]-h[p1]*yh[p2])/(h[p2]-h[p1])
// #define min(x,y) (x<y?x:y)
// #define max(x,y) (x>y?x:y)
	

using namespace voro;
using namespace std;
	
class Weightmap;
class LSbox;
class mathMethods;

class grainhdl{

	int ngrains;
	double dt;
	double h;
	
	int realDomainSize;
	int ngridpoints;
	int grid_blowup;
	
	int Mode;   	
		
	
	public:
	
	mathMethods* mymath;
	int loop;
	vector<double> totalenergy;
	vector<int> nr_grains;
	
	double tubeRadius;
	double *ST;
	double *part_pos;	
	double delta;
	double *bunge;
	double deviation;
	
	vector<LSbox*> grains;		

	LSbox* boundary;
	grainhdl();
	~grainhdl();
		
	void setSimulationParameter(); 
	void generateRandomEnergy();
	
	void VOROMicrostructure();
	void readMicrostructurefromVertex();
	
	void find_neighbors();
	 
	void convolution();
	void save_conv_step();
	void comparison_box();
	void updateSecondOrderNeighbors();
	void determineIDs();

	void level_set();
	void redistancing();
	
	void run_sim();
	void save_sim();
	void clear_mem(); 
	void save_texture();
	void plot_contour();
	
	void compute_Boundary_Energy();
	void construct_boundary();
	void saveAllContourLines();
	void saveAllContourEnergies();
	void switchDistancebuffer();
	
	void saveSpecialContourEnergies(int id);

	

	
// 	wrapper functions:
	
	inline long get_ngrains(){ return ngrains; }
	inline int get_realDomainSize() { return realDomainSize; }
	inline int get_ngridpoints() { return ngridpoints; }
	inline double get_h() { return h; }
	inline int get_grid_blowup() { return grid_blowup; }
	inline int get_loop() { return loop; }
	inline double get_dt() { return dt; }
};
#endif
