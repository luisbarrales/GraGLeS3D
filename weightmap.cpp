#include "weightmap.h"
#include "Settings.h"
Weightmap::mapkey::mapkey(LSbox** IDs) :
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
	return this->loadWeights(&IDs[0], IDs.size(), me, ST);
}

double Weightmap::loadWeights(LSbox** IDs, int length, LSbox* me, double* ST)
{
	if(length == 2)
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
	else return m_pHandler->hagb;
}

double Weightmap::isTriplePoint(vector<LSbox*> IDs){
	Weightmap::mapkey key_tuple(&IDs[0]);
	std::map<mapkey, double>::iterator it = m_Weights.find(key_tuple);
	if (it == m_Weights.end())	//If value is not present, calculate and store it.
		{
			return 0.0;
		}
	else return (*it).second;
}


double Weightmap::computeWeights(Weightmap::mapkey rep, LSbox* me, double* ST)
{
	double sigma;
	double gamma[3];
	double gamma_hagb = m_pHandler->hagb;//0.6;
	double theta_ref = 15.0 * PI /180;
	double theta_mis;


	if(Settings::MicrostructureGenMode == 2){
		gamma[0] = ST[(me->get_id() - 1)
				+ (m_pHandler->get_ngrains() * (rep.first->get_id() - 1))];
		gamma[1] = ST[(rep.first->get_id() - 1)
				+ (m_pHandler->get_ngrains() * (rep.second->get_id() - 1))];
		gamma[2] = ST[(me->get_id() - 1)
				+ (m_pHandler->get_ngrains() * (rep.second->get_id() - 1))];
	}	
	
	if(Settings::MicrostructureGenMode == 1){
 		if(Settings::ResearchMode == 0){
			theta_mis = me->mis_ori(rep.first);		
			if (theta_mis <= theta_ref)	gamma[0] = gamma_hagb * ( theta_mis / theta_ref) * (1.0 - log( theta_mis / theta_ref));
			else gamma[0] = gamma_hagb;
			
			
			theta_mis = rep.first->mis_ori(rep.second);			
			if (theta_mis <= theta_ref) gamma[1] = gamma_hagb * ( theta_mis / theta_ref) * (1.0 - log( theta_mis / theta_ref));
			else gamma[1] = gamma_hagb;
			
			
			theta_mis = me->mis_ori(rep.second);
			if (theta_mis <= theta_ref) gamma[2] = gamma_hagb * ( theta_mis / theta_ref) * (1.0 - log( theta_mis / theta_ref));
			else gamma[2] = gamma_hagb;
 		}
 		else {
 			theta_ref = 42 * PI /180;
 			theta_mis = me->mis_ori(rep.first);
			if (theta_mis <= theta_ref)	gamma[0] = 0.3;
			else gamma[0] = gamma_hagb;

			theta_mis = rep.first->mis_ori(rep.second);
			if (theta_mis <= theta_ref) gamma[1] = 0.3;
			else gamma[1] = gamma_hagb;


			theta_mis = me->mis_ori(rep.second);
			if (theta_mis <= theta_ref) gamma[2] = 0.3;
			else gamma[2] = gamma_hagb;

 		}


			
	}
	sigma = gamma[0] - gamma[1] + gamma[2];
	
	if(sigma < 0.0) {
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
