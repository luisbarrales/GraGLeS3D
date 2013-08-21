#ifndef WEIGHTMAP_H
#define WEIGHTMAP_H
#include "levelsetproject.h"

typedef std::map<int, double* > storage_map;
typedef std::map<int, storage_map* > inner_map;
typedef std::map<int, inner_map* > outer_map;

class weightmap{
	
// 	map<int , map<int, map<int,double*>* >* > maptable;
	outer_map weights_table;
		
	public:
	friend class LSbox;
	weightmap();
	~weightmap();
	LSbox** find_representer(int length, LSbox ***ID,int i, int j);
	double load_weights(int length, double *ST, LSbox*** ID, int i, int j);
	double* compute_weights( double* ST, int* ids);
	void add_weights(int* ids, double* sigma);
	void add_weights(int* ids, outer_map::iterator it, double* sigma);
	void add_weights(int* ids, inner_map::iterator it2, double* sigma);
};

#endif