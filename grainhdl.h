#ifndef GRAINHDL_h
#define GRAINHDL_h

#include levelsetproject.h
using namespace voro;
using namespace std;

class matrix;
class weightmap;
class LSbox;

class grainhdl{

	long ngrains;
	double dt;
	double h;
	int realDomainSize;
	int ngridpoints;
	int grid_blowup;
	int Mode;
	
	int *gridIDs;
	ID = new LSbox**[3];
	
	weightmap my_weights;
	
	
	std::list<matrix> domains, domains_copy;
    std::list<matrix>::iterator it, itc, it_domain;
	
	vector<int> nr_grains(TIMESTEPS+1);
	vector<LSbox*> buffer;
	
	LSbox ***ID;
	double *ST;
	double *part_pos;
	
	vector<LSbox*> grains;gridpoint
    vector<LSbox*>::iterator it2,it2c;
	
	public:
	
	void readInit(); //to do!!!
	void setSimulationParameter(); 
	void generateRandomEnergy();
	
	void VOROMicrostructure();
	void readMicrostructurefromVertex();
	
	void save_conv_step();
	
	
	void my_sim.run_sim();
	void my_sim.save_sim();
	void my_sim.clear_mem(); 
	
	
// 	wrapper functions:
	
	inline long ngrains(){ return ngrains; }
	inline int realDomainSize() { return realDomainSize; }
	inline int ngridpoints() { return ngridpoints; }

};