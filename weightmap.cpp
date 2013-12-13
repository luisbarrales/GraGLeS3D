#include "weightmap.h"
Weightmap::mapkey::mapkey(vector<LSbox*> IDs) :
	first(NULL),
	second(NULL),
	third(NULL)
{
	//Debug checks
	if(IDs.size() != 3)
	{
		//TODO: Replace cout with logger. Possibly investigate case of 2 IDs incoming.
		cout<<"IDs list contains more or less than 3 ids...\n";
	}
	//Sort & store first three elements in the member fields.
	//Sort of 3 elements based on a binary decision tree.
	if(IDs[0] < IDs[1])
	{
		if(IDs[1] < IDs[2])		//{0,1,2}
		{first = IDs[0]; second = IDs[1]; third = IDs[2];}
		else
		{
			if(IDs[0] < IDs[2])	//{0,2,1}
			{first = IDs[0]; second = IDs[2]; third = IDs[1];}
			else 				//{2,0,1}
			{first = IDs[2]; second = IDs[0]; third = IDs[1];}
		}
	}
	else
	{
		if(IDs[0] < IDs[2])		//{1,0,2}
		{first = IDs[1]; second = IDs[0]; third = IDs[2];}
		else
		{
			if(IDs[1] < IDs[2])	//{1,2,0}
			{first = IDs[1]; second = IDs[2]; third = IDs[0];}
			else				//{2,1,0}
			{first = IDs[2]; second = IDs[1]; third = IDs[0];}
		}
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
	else if (this->second > other.second)
		return false;
	else if (this->third < other.third)
		return true;
	else
		return false;
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
	Weightmap::mapkey key_tuple(IDs);
	double weight;

	//TODO: Why ?
	if (key_tuple.first == key_tuple.second || key_tuple.second == key_tuple.third)
		return 1.0;
	else
	{
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
	}
	return weight;
}

double Weightmap::computeWeights(Weightmap::mapkey rep, LSbox* me, double* ST)
{
	double sigma;
	double gamma[3];
	double gamma_hagb = 0.6;
// 	double theta_ref = 15.0;
// 	double theta_mis;
	double drag = 0.5;

	if (rep.first == rep.second || rep.second == rep.third)
	{
		sigma = 1.0;
	}
	else
	{
		gamma[0] = ST[(rep.first->get_id() - 1)
				+ (m_pHandler->get_ngrains() * (rep.second->get_id() - 1))];
		gamma[1] = ST[(rep.first->get_id() - 1)
				+ (m_pHandler->get_ngrains() * (rep.third->get_id() - 1))];
		gamma[2] = ST[(rep.second->get_id() - 1)
				+ (m_pHandler->get_ngrains() * (rep.third->get_id() - 1))];

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

// w�hle gamma oder lade aus ST-Feld

		if (me == rep.first)
			sigma = gamma[0] + gamma[1] - gamma[2];
		else if (me == rep.second)
			sigma = gamma[0] - gamma[1] + gamma[2];
		else if (me == rep.third)
			sigma = -gamma[0] + gamma[1] + gamma[2];
//
//		if (std::isnan(sigma[0]) || std::isnan(sigma[1])
//				|| std::isnan(sigma[2]))
//		{
//			cout << "IS NAN " << endl;
//			cout << sigma[0] << "\t" << sigma[1] << "\t" << sigma[2] << "\t";
//			char buffin;
//			cin >> buffin;
//		}
//
//		if (rep.first == 1 && rep.second == 3)
//			sigma[3] = drag;
//		else
//			sigma[3] = 1.0;
//		sigma[3] = 1.0; // could be used as a dragfactor
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
