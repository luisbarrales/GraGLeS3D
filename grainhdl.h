#ifndef GRAINHDL_h
#define GRAINHDL_h

#include "ggLS.h"

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
	
	LSbox* zeroBox;
	
	std::list<domainCl> domains;
	std::list<domainCl> domains_copy;
    	
	vector<int> nr_grains;
	vector<LSbox*> buffer;
	
	LSbox ***ID;
	double *ST;
	double *part_pos;
	
	public:
	
	grainhdl();
	~grainhdl();
	void readInit(); //to do!!!
	void setSimulationParameter(); 
	void generateRandomEnergy();
	
	void VOROMicrostructure();
	void readMicrostructurefromVertex();
	
	 
	void convolution();
	void save_conv_step();
	void comparison_domain();
	void comparison_box();
	void swap_grains();
	void redistancing();
	
	void run_sim();
	void save_sim();
	void clear_mem(); 
	
	
// 	wrapper functions:
	
	inline long get_ngrains(){ return ngrains; }
	inline int get_realDomainSize() { return realDomainSize; }
	inline int get_ngridpoints() { return ngridpoints; }

};
#endif