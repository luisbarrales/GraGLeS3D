#include "weightmap.h"
Weightmap::mapkey::mapkey(vector<LSbox*> IDs) :
	first(NULL),
	second(NULL)
{
	//Sort & store first three elements in the member fields.
	//Sort of 3 elements based on a binary decision tree.
	if(IDs[0] < IDs[1])
	{
		first = IDs[0];
		second = IDs[1];
	}
	else{
		second = IDs[0];
		first = IDs[1];	
	}	
}

bool Weightmap::mapkey::operator<(const Weightmap::mapkey & other) const
{
	if(this->first < other.first)
		return true;
	else if (this->first > other.first)
		return false;
	else if (this->second < other.second)
		return true;
	else return false;	
}

Weightmap::Weightmap(grainhdl* owner) :
	m_pHandler(owner)
{
}
Weightmap::~Weightmap()
{
}

double Weightmap::loadWeights(vector<LSbox*> IDs, LSbox* me, double* ST)
{
	if(IDs.size() == 2)
	{
			Weightmap::mapkey key_tuple(IDs);
		double weight;
		std::map<mapkey, double>::iterator it = m_Weights.find(key_tuple);
		
		if (it == m_Weights.end())	//If value is not present, calculate and store it.
		{
			weight = computeWeights(key_tuple, me, ST);
			m_Weights[key_tuple] = weight;
		}
		else						//If present just fetch.
		{
			weight = (*it).second;
		}

		return weight;
	}
// 	else if (IDs.size() >= 2)
// TODO
	else return m_pHandler->hagb;

}

double Weightmap::computeWeights(Weightmap::mapkey rep, LSbox* me, double* ST)
{
	
	double sigma;
	double gamma[3];
	double gamma_hagb = m_pHandler->hagb;//0.6;
	double theta_ref = 15.0 * PI /180;
	double theta_mis;
// 	double drag = 0.5;

	if(MODE == 2){
		gamma[0] = ST[(me->get_id() - 1)
				+ (m_pHandler->get_ngrains() * (rep.first->get_id() - 1))];
		gamma[1] = ST[(rep.first->get_id() - 1)
				+ (m_pHandler->get_ngrains() * (rep.second->get_id() - 1))];
		gamma[2] = ST[(me->get_id() - 1)
				+ (m_pHandler->get_ngrains() * (rep.second->get_id() - 1))];
	}	
	
	if(MODE == 1){
		if(!MOBILITY){
			theta_mis = me->mis_ori(rep.first);		
			if (theta_mis <= theta_ref)	gamma[0] = gamma_hagb * ( theta_mis / theta_ref) * (1.0 - log( theta_mis / theta_ref));
			else gamma[0] = gamma_hagb;
			
			
			theta_mis = rep.first->mis_ori(rep.second);			
			if (theta_mis <= theta_ref) gamma[1] = gamma_hagb * ( theta_mis / theta_ref) * (1.0 - log( theta_mis / theta_ref));
			else gamma[1] = gamma_hagb;
			
			
			theta_mis = me->mis_ori(rep.second);
			if (theta_mis <= theta_ref) gamma[2] = gamma_hagb * ( theta_mis / theta_ref) * (1.0 - log( theta_mis / theta_ref));
			else gamma[2] = gamma_hagb;
			
// 			gamma[0] = me->getGBEnergy(rep.first);
// 			gamma[1] = rep.first->getGBEnergy(rep.second);
// 			gamma[2] = me->getGBEnergy(rep.second);
		}
		else {
			gamma[0] = me->getGBEnergyTimesGBMobility(rep.first);
			gamma[1] = rep.first->getGBEnergyTimesGBMobility(rep.second);
			gamma[2] = me->getGBEnergyTimesGBMobility(rep.second);
		}
	}
	sigma = gamma[0] - gamma[1] + gamma[2];
	
	if(sigma < 0.05) {
		cout << sigma << endl;
		cout <<"gamma in weighmap :  ";
		cout << gamma[0] << "    " << gamma[1] << "    "  << gamma[2] << endl;
		sigma=0.05;
		char buf;
		cin >> buf;
	}
	return sigma;
}

void Weightmap::plotWeightmap(int length, LSbox*** ID, double* ST,
		LSbox* zeroBox)
{

	ofstream myfile;
	myfile.open("weightmap.txt");
//	for (it = weights_table.begin(); it != weights_table.end(); it++)
//		for (it2 = (*it).second->begin(); it2 != (*it).second->end(); it2++)
//			for (it3 = (*it2).second->begin(); it3 != (*it2).second->end();
//					it3++)
//			{
//				int ids[3] = { (*it).first, (*it2).first, (*it3).first };
//				double* weight = (*it3).second;
//				cout << rep.first << "\t" << rep.second << "\t" << rep.third
//						<< "\t" << weight[0] << "\t" << weight[1] << "\t"
//						<< weight[2] << "\n";
//				myfile << rep.first << "\t" << rep.second << "\t" << rep.third
//						<< "\t" << weight[0] << "\t" << weight[1] << "\t"
//						<< weight[2] << "\n";
//			}
	myfile.close();
}
