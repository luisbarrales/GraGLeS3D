#include "weightmap.h"

weightmap::weightmap(){}
weightmap::~weightmap(){}

LSbox** weightmap::find_representer(LSbox ***ID,int i, int j){
	LSbox** rep = new LSbox* [3];
	int length = sqrt(sizeof(ID[0])/sizeof(ID[0][0]));
	if( ID[0][i*length +j]->get_id() < ID[1][i*length +j]->get_id() ){
		rep[0]=ID[0][i*length +j];
		rep[2]=ID[1][i*length +j];
	}
	else {
		rep[2]=ID[0][i*length +j];
		rep[0]=ID[1][i*length +j];
	}
	if( rep[2]->get_id() < ID[2][i*length +j]->get_id() ){
		rep[1]=rep[2];
		rep[2]=ID[2][i*length +j];
	}
	else {
		if( rep[0]->get_id() > ID[2][i*length +j]->get_id() ){
			rep[1]=rep[0];
			rep[0]=ID[2][i*length +j];
		}
		else {
			rep[1]=ID[2][i*length +j];
		}
	}
	return rep;
}

void weightmap::add_weights(int& ids, double* sigma){
	map<int, map<int, double*>* >* temp = new map<int, map<int, double*>* > ;
	maptable[ids[0]] = temp;
	(*temp)[ids[1]] = new map<int, double*>;
	(*(*temp)[ids[1]])[ids[2]] = sigma;
}

void weightmap::add_weights(int& ids, map<int, map<int, map<int, double*>* >* > :: iterator it, double* sigma){
	(**it)[ids[1]] = new map<int, double*>;
	(*(**it)[ids[1]])[ids[2]] = sigma;
}

void weightmap::add_weights(int& ids, map<int, map<int, double*>* >* :: iterator it2, double* sigma){
	(**it2)[ids[2]] = sigma;
}



double weightmap::load_weights(double* ST, int*** ID, int i, int j){
	LSbox** rep = find_representer(ID,i,j);	
	double* sigma;
	map<int, map<int, map<int, double*>* >* > :: iterator it;
	map<int, map<int, double*>* >* :: iterator it2;
	map<int, double*>* :: iterator it3;
	double ids[3]={(*rep[0]).get_id(),(*rep[1]).get_id(), (*rep[2]).get_id()};
	
	it = maptable.find(ids[0]);
	if (it == maptable.end() ){
		sigma = compute_weights(ST, ids);
		add_weights(ids, sigma);
	}
	else {
		it2 = **it.find(ids[1]);
		if (it2 == **it.end() ){
			sigma = compute_weights(ST, ids);
			add_weights (ids, it, sigma);
		}
		else{
			it3= **it2.find(ids[2]);
			sigma = compute_weights(ST, ids);
			add_weights(ids, it2, sigma);
		}
	}
// 	if exist
	int length = sqrt(sizeof(ID[0])/sizeof(ID[0][0]));
	int i=0;
	while(ids[i] != ID[0][i*length+j].get_id()){
		i++;
	}
	return sigma[i];
	
}

double* weightmap::compute_weights(double *ST, LSbox** rep){
	double* sigma = new double[3];
	double gamma[3];
	gamma[0] = ST[ (ids[0]-1) + (PARTICLES* ( ids[1]-1) ) ];
	gamma[1] = ST[ (ids[0]-1) + (PARTICLES* ( ids[2]-1) ) ];
	gamma[2] = ST[ (ids[1]-1) + (PARTICLES* ( ids[2]-1) ) ];
	// wähle gamma oder lade aus ST-Feld
	sigma[0]= 	gamma[0] + gamma[1] - gamma[2];
	sigma[1]= 	gamma[0] - gamma[1] + gamma[2];
	sigma[2]= -	gamma[0] + gamma[1] + gamma[2];
// 	if exist

// speichere sigma
	return sigma;
}