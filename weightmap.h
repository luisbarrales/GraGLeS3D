#ifndef WEIGHTMAP_H
#define WEIGHTMAP_H
#include "levelsetproject.h"
class weightmap{
	
	map<int , map<int, map<int,double*>* >* > maptable;
		
	public:
	friend class LSbox;
	weightmap();
	~weightmap();
	LSbox** find_representer(LSbox ***ID,int i, int j);
	double load_weights(double *ST, LSbox** rep);
	double* compute_weights(double* ST, int& ids);
	void add_weights(int& ids, double* sigma);
	void add_weights(int& ids, map<int, map<int, map<int, double*>* >* > :: iterator it, double* sigma);
	void add_weights(int& ids, map<int, map<int, double*>* >* :: iterator it2, double* sigma);
};

#endif