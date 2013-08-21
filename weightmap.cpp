#include "weightmap.h"

weightmap::weightmap(){}
weightmap::~weightmap(){}



LSbox** weightmap::find_representer(int length, LSbox ***ID,int i, int j){
	LSbox** rep = new LSbox* [3];
	
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



void weightmap::add_weights(int* ids, double* sigma){
	inner_map* temp = new inner_map;	
	weights_table[ids[0]] = temp;
	(*temp)[ids[1]] = new storage_map;
	(*(*temp)[ids[1]])[ids[2]] = sigma;
}


void weightmap::add_weights(int* ids, outer_map::iterator it, double* sigma){
	(*(*it).second)[ids[1]] = new storage_map;
	(*(*(*it).second)[ids[1]])[ids[2]] = sigma;
}

void weightmap::add_weights(int* ids, inner_map::iterator it2, double* sigma){
	(*(*it2).second)[ids[2]] = sigma;
}



double weightmap::load_weights(int length, double* ST, LSbox*** ID, int i, int j){
	LSbox** rep = find_representer(length,ID,i,j);	
	double* sigma;
	outer_map::iterator it;
	inner_map::iterator it2;
	storage_map::iterator it3;
	int ids[3]={(*rep[0]).get_id(),(*rep[1]).get_id(), (*rep[2]).get_id()};
	
	it = weights_table.find(ids[0]);
	if (it == weights_table.end() ){		
		sigma = compute_weights(ST, ids);
		add_weights(ids, sigma);
	}
	else {
		it2 = (*(*it).second).find(ids[1]);
		if (it2 == (*(*it).second).end() ){
			sigma = compute_weights(ST, ids);
			add_weights (ids, it, sigma);
		}
		else{
			it3= (*(*it2).second).find(ids[2]);
			if(it3 == (*(*it2).second).end()){
				sigma = compute_weights(ST, ids);
				add_weights(ids, it2, sigma);
			}
			else { sigma = (*it3).second;
			cout << "read map :";
			utils::print_2dim_array(ids,1,3);
			utils::print_2dim_array(sigma,1,3);
			}
		}
	}
// 	if exist
	int ii=0;
	while(ids[ii] != ID[0][i*length+j]->get_id()){
		ii++;
	}
	return sigma[ii];
	
}

double* weightmap::compute_weights(double *ST,  int* ids){
	double* sigma = new double[3];
	double gamma[3];
	gamma[0] = ST[ (ids[0]-1) + (PARTICLES* ( ids[1]-1) ) ];
	gamma[1] = ST[ (ids[0]-1) + (PARTICLES* ( ids[2]-1) ) ];
	gamma[2] = ST[ (ids[1]-1) + (PARTICLES* ( ids[2]-1) ) ];
	
	
	// wähle gamma oder lade aus ST-Feld
	sigma[0]= 	gamma[0] + gamma[1] - gamma[2];
	sigma[1]= 	gamma[0] - gamma[1] + gamma[2];
	sigma[2]= -	gamma[0] + gamma[1] + gamma[2];
	
	if( sigma[0]!= 1.0 ) {
		cout << "copute weights" << endl;
		utils::print_2dim_array(gamma,1,3);
		utils::print_2dim_array(ids,1,3);
		utils::print_2dim_array(sigma,1,3);
		char c;
		cin >> c;
	}
// 	if exist

// speichere sigma
	return sigma;
}