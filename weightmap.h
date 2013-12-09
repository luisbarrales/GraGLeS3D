#ifndef WEIGHTMAP_H
#define WEIGHTMAP_H
#include "ggLS.h"

class grainhdl;

struct mapkey {
  LSbox* first;
  LSbox* second;
  LSbox* third; 
  
   bool operator<( const mapkey & c ) const {
	if (this->first < c.first)
	  return true;
	
	else if(this->first == c.first){
	  if (this->second < c.second){
	    return true;
	  }
	  else if ( this->second == c.second){
	      
	    if (this->third < c.third)
		return true;
	    else return false;
	  
	  }
	  return false;
	}
	return false;
  }

};

class weightmap{
	
	grainhdl* handler;
	std::map<mapkey, double> weights_table;
		
	public:
	friend class LSbox;
	weightmap(grainhdl* handler);
	~weightmap();
	void find_representer(mapkey& rep, vector<LSbox*> IDs);
	double load_weights(vector<LSbox*> IDs,LSbox* me, double* ST,mapkey rep);
	double* compute_weights( mapkey rep, LSbox* me, double* ST);
	void plot_weightmap(int length, LSbox*** ID, double* ST, LSbox * zeroBox);
};

#endif