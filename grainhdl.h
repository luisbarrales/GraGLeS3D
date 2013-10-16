#ifndef GRAINHDL_h
#define GRAINHDL_h

#include "ggLS.h"


#define xsect(p1,p2) (h[p2]*xh[p1]-h[p1]*xh[p2])/(h[p2]-h[p1])
#define ysect(p1,p2) (h[p2]*yh[p1]-h[p1]*yh[p2])/(h[p2]-h[p1])
// #define min(x,y) (x<y?x:y)
// #define max(x,y) (x>y?x:y)
	

using namespace voro;
using namespace std;

class domainCl;
class weightmap;
class LSbox;

class grainhdl{

	int ngrains;
	double dt;
	double h;
	int realDomainSize;
	int ngridpoints;
	int grid_blowup;
	int Mode;
	int compare_mod;
	int loop;
	
	int *gridIDs;
		
	weightmap* my_weights;
	std::list<double*> vol_list;
	std::list<double*> modf;
	std::list<double*> energy;
	std::list<domainCl> domains;
	std::list<domainCl> domains_copy;
    	
	vector<int> nr_grains;
	
	LSbox ***ID;
	double *ST;
	double *part_pos;
	
	public:
	
			
	LSbox* zeroBox;
	grainhdl();
	~grainhdl();
	void readInit(); 
	void setSimulationParameter(); 
	void generateRandomEnergy();
	
	void VOROMicrostructure();
	void readMicrostructurefromVertex();
	
	void find_neighbors();
	 
	void convolution();
	void save_conv_step();
	void comparison_domain();
	void comparison_box();
	void swap_grains();
	void redistancing();
	
	void run_sim();
	void save_sim();
	void clear_mem(); 
	
	void compute_grain_vol();
	int  conrec();


	
// 	wrapper functions:
	
	inline long get_ngrains(){ return ngrains; }
	inline int get_realDomainSize() { return realDomainSize; }
	inline int get_ngridpoints() { return ngridpoints; }
	inline double get_h() { return h; }
	inline int get_grid_blowup() { return grid_blowup; }
};
#endif