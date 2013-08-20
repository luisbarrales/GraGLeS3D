#ifndef WEIGHTMAP_H
#define WEIGHTMAP_H
#include "levelsetproject.h"
class weightmap{
	
	map<int , map<int, map<int,double*>* >* > maptable;
		
	public:
	friend class LSbox;
	weightmap();
	~weightmap();
	LSbox** find_representer(int ***ID,int i, int j);
	double* load_weights(LSbox** rep);
	double* compute_weights(LSbox** rep);
};

#endif