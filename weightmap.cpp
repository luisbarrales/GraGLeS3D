#include "weightmap.h"

weightmap::weightmap(grainhdl* handler): handler(handler){}
weightmap::~weightmap(){}



void weightmap::find_representer(LSbox** rep, int length, LSbox ***ID,int i, int j){
		
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



double weightmap::load_weights(int length, double* ST, LSbox*** ID, int i, int j, int id){
	LSbox** rep= new LSbox* [3];
	find_representer(rep, length,ID,i,j);	
	double* sigma;
	outer_map::iterator it;
	inner_map::iterator it2;
	storage_map::iterator it3;
	
	int ids[3]={(*rep[0]).get_id(),(*rep[1]).get_id(), (*rep[2]).get_id()};		
		if(ids[0]==ids[1] || ids[1]==ids[2]) return 1.0;
		else{
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
// 						cout << "entry  not exist "<< endl; 
// 						utils::print_2dim_array(ids,1,3);
// 						utils::print_2dim_array(sigma,1,3);
						add_weights(ids, it2, sigma);
					}
					else { 
						sigma = (*it3).second;
					}
				}
			}
		}
	int ii=0;
	while(ids[ii] != id ){
		ii++;
	}
	delete [] rep;
	return (sigma[ii]*sigma[3]);
	
}

double* weightmap::compute_weights(double *ST,  int* ids){
	double* sigma = new double[4];
	double gamma[3];
	double gamma_hagb = 0.6;
	double theta_ref = 15.0;
	double theta_mis;
	double drag = 0.5;
	
	if (ids[0] == ids[1]  ||  ids[1] == ids[2]){
		sigma[0]= 1.0;
		sigma[1]= 1.0;
		sigma[2]= 1.0;	
		sigma[3]= 1.0;
	}
	else{
		gamma[0] = ST[ (ids[0]-1) + (handler->get_ngrains() * ( ids[1]-1) ) ];
		gamma[1] = ST[ (ids[0]-1) + (handler->get_ngrains() * ( ids[2]-1) ) ];
		gamma[2] = ST[ (ids[1]-1) + (handler->get_ngrains() * ( ids[2]-1) ) ];
		
// 		theta_mis = (*handler->grains)[ids[0]]->mis_ori((*handler->grains)[ids[1]]);
// 		if (theta_mis <= theta_ref)	gamma[0] = gamma_hagb * ( theta_mis / theta_ref) * (1.0 - log( theta_mis / theta_ref));
// 		else gamma[0] = gamma_hagb;
// 		
// 		theta_mis = (*handler->grains)[ids[0]]->mis_ori((*handler->grains)[ids[2]]);
// 		if (theta_mis <= theta_ref) gamma[1] = gamma_hagb * ( theta_mis / theta_ref) * (1.0 - log( theta_mis / theta_ref));
// 		else gamma[1] = gamma_hagb;
// 		
// 		theta_mis = (*handler->grains)[ids[1]]->mis_ori((*handler->grains)[ids[2]]);
// 		if (theta_mis <= theta_ref) gamma[2] = gamma_hagb * ( theta_mis / theta_ref) * (1.0 - log( theta_mis / theta_ref));
// 		else gamma[2] = gamma_hagb;
		
		// wähle gamma oder lade aus ST-Feld
		sigma[0]= 	gamma[0] + gamma[1] - gamma[2];
		sigma[1]= 	gamma[0] - gamma[1] + gamma[2];
		sigma[2]= -	gamma[0] + gamma[1] + gamma[2];
// 		if(ids[0]==1 && ids[1]==3) sigma[3]=drag;
// 		else sigma[3]=1.0;
		sigma[3]=0.1;
	}

	return sigma;
}

void weightmap::plot_weightmap(int length, LSbox*** ID, double* ST, LSbox* zeroBox){
	domainCl *temp = new domainCl(length,length);
	domainCl *id_0 = new domainCl(length,length);
	domainCl *id_1 = new domainCl(length,length);
	domainCl *id_2 = new domainCl(length,length);
	double weight;
	for(int i=0; i<length; i++)
		for(int j=0; j<length; j++){
			if( ID[0][i*length +j] != zeroBox ){
				weight = load_weights(length,ST,ID,i,j, ID[0][i*length +j]->get_id() );
// 				cout << weight << endl;
				if( weight != 1.0 )(*temp)[i][j] = weight;
				else (*temp)[i][j] = 0.0;
				(*id_0)[i][j]= ID[0][i*length+j]->get_id();
				(*id_1)[i][j]= ID[1][i*length+j]->get_id();
				(*id_2)[i][j]= ID[2][i*length+j]->get_id();
			}
			else (*temp)[i][j]=0.0;
		}
		(*id_0).save_domainCl("id_0.gnu");
		(*id_1).save_domainCl("id_1.gnu");
		(*id_2).save_domainCl("id_2.gnu");
		(*temp).save_domainCl("weights.gnu");
		
		delete temp;
		delete id_0;
		delete id_1;
		delete id_2;
		
}
