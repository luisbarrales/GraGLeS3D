#ifndef WEIGHTMAP_H
#define WEIGHTMAP_H
#include "ggLS.h"

class grainhdl;

class weightmap{
	
	grainhdl* handler;
	std::map<key, double> weights_table;
		
	public:
	friend class LSbox;
	weightmap(grainhdl* handler);
	~weightmap();
	void find_representer(key& rep, vector<LSbox*> IDs);
	double load_weights(vector<LSbox*> IDs,LSbox* me, double* ST,key rep);
	double* compute_weights( key rep, LSbox* me, double* ST);
	void plot_weightmap(int length, LSbox*** ID, double* ST, LSbox * zeroBox);
};

#endif