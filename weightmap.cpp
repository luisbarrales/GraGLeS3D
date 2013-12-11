#include "weightmap.h"

weightmap::weightmap(grainhdl* handler): handler(handler){}
weightmap::~weightmap(){}
 
 


void weightmap::find_representer(mapkey& rep, vector<LSbox*> IDs){
//
  if ( IDs.size() <= 3){
	if(IDs[0]->get_id()< IDs[1]->get_id()){
		rep.first=IDs[0];
		rep.third=IDs[1];
	}
	else {
		rep.third=IDs[0];
		rep.first=IDs[1];
	}
	if( rep.third->get_id()< IDs[2]->get_id()){
		rep.second=rep.third;
		rep.third=IDs[2];
	}
	else {
		if( rep.first->get_id()> IDs[2]->get_id()){
			rep.second=rep.first;
			rep.first=IDs[2];
		}
		else {
			rep.second=IDs[2];
		}
	}
  }
  else cout << "More than 3 IDs" ; 
}


double weightmap::load_weights(vector<LSbox*> IDs, LSbox* me, double* ST){
    mapkey rep;
    find_representer(rep,IDs);	
    double sigma;
  
    
    if(rep.first==rep.second || rep.second==rep.third) 
      return 1.0;
    else{
		std::map<mapkey,double>::iterator it = weights_table.find(rep);
		if (it == weights_table.end() ){		
		sigma = compute_weights(rep, me, ST);
		weights_table[rep]= sigma;
		}
    }
    return sigma;
}

double weightmap::compute_weights(mapkey rep,  LSbox* me, double* ST){
	double sigma;
	double gamma[3];
	double gamma_hagb = 0.6;
// 	double theta_ref = 15.0;
// 	double theta_mis;
	double drag = 0.5;
	
	if (rep.first == rep.second  ||  rep.second == rep.third){
		sigma= 1.0;
	}
	else{
		gamma[0] = ST[ (rep.first->get_id()-1) + (handler->get_ngrains() * ( rep.second->get_id()-1) ) ];
		gamma[1] = ST[ (rep.first->get_id()-1) + (handler->get_ngrains() * ( rep.third->get_id()-1 ) ) ];
		gamma[2] = ST[ (rep.second->get_id()-1)+ (handler->get_ngrains() * ( rep.third->get_id()-1 ) ) ];
		
// 		theta_mis = (*handler->grains)[rep.first]->mis_ori((*handler->grains)[rep.second]);
// 		if (theta_mis <= theta_ref)	gamma[0] = gamma_hagb * ( theta_mis / theta_ref) * (1.0 - log( theta_mis / theta_ref));
// 		else gamma[0] = gamma_hagb;
// 		
// 		theta_mis = (*handler->grains)[rep.first]->mis_ori((*handler->grains)[rep.third]);
// 		if (theta_mis <= theta_ref) gamma[1] = gamma_hagb * ( theta_mis / theta_ref) * (1.0 - log( theta_mis / theta_ref));
// 		else gamma[1] = gamma_hagb;
// 		
// 		theta_mis = (*handler->grains)[rep.second]->mis_ori((*handler->grains)[rep.third]);
// 		if (theta_mis <= theta_ref) gamma[2] = gamma_hagb * ( theta_mis / theta_ref) * (1.0 - log( theta_mis / theta_ref));
// 		else gamma[2] = gamma_hagb;
		
		// wähle gamma oder lade aus ST-Feld
		
		if ( me == rep.first)
		  sigma = 	gamma[0] + gamma[1] - gamma[2];
		else if ( me == rep.second)
		  sigma= 	gamma[0] - gamma[1] + gamma[2];
		else if ( me == rep.third)
		  sigma= -	gamma[0] + gamma[1] + gamma[2];
/*		
		if (  std::isnan(sigma[0])||std::isnan(sigma[1])||std::isnan(sigma[2]))	{
		  cout << "IS NAN " << endl;
		  cout << sigma[0] << "\t" << sigma[1] << "\t"<< sigma[2] << "\t";
		  char buffin;
		  cin >> buffin;
		}

		*/
		// 		if(rep.first==1 && rep.second==3) sigma[3]=drag;
// 		else sigma[3]=1.0;
// 		sigma[3]=1.0; // could be used as a dragfactor
	}

	return sigma;
}

void weightmap::plot_weightmap(int length, LSbox*** ID, double* ST, LSbox* zeroBox){
	
	ofstream myfile;
	myfile.open ("weightmap.txt");/*
	for (it=weights_table.begin(); it!=weights_table.end(); it++)
		for (it2=(*it).second->begin(); it2!=(*it).second->end() ;it2++)
			for (it3=(*it2).second->begin(); it3!=(*it2).second->end() ;it3++){
// 				int ids[3]={(*it).first, (*it2).first, (*it3).first};
// // 				double* weight = (*it3).second;
// 				cout << rep.first << "\t" << rep.second << "\t" << rep.third << "\t" << weight[0] << "\t" << weight[1] << "\t" << weight[2] << "\n";
// 				myfile << rep.first << "\t" << rep.second << "\t" << rep.third << "\t" << weight[0] << "\t" << weight[1] << "\t" << weight[2] << "\n";
			}*/
	
	myfile.close();
		
}
