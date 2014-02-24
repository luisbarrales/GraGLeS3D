#include "Settings.h"
#include "rapidxml.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>

using namespace std;
using namespace rapidxml;

//Initializing the static setting variables
unsigned int Settings::NumberOfParticles = 0;
unsigned int Settings::NumberOfTimesteps = 0;
unsigned int Settings::AnalysysTimestep = 0;
unsigned int Settings::DiscreteSamplingRate = 0;
unsigned int Settings::DomainBorderSize = 0;
E_MICROSTRUCTURE_GEN_MODE Settings::MicrostructureGenMode = E_INVALID_VALUE;
string Settings::ReadFromFilename;
double Settings::HAGB = 0.0;
bool Settings::UseMobilityFactor = false;
bool Settings::IsIsotropicNetwork = false;
bool Settings::UseTexture = false;
bool Settings::ExecuteInParallel = false;

void Settings::initializeParameters(string filename)
{
	if(0 == filename.compare(""))
		filename = string("parameters.xml");
	ifstream file(filename);
	if(file.fail())
	{
		cout<<"Unable to locate simulations parameters. Will now halt !" << endl;
		exit(2);
	}
	stringstream contents;
	contents << file.rdbuf();
	string xmlDocument(contents.str());

	xml_document<> tree;
	try
	{
		tree.parse<0>(&xmlDocument[0]);
	}
	catch (parse_error Error)
	{
		cout << filename.c_str() << "is not a valid XML!" << endl;
		cout << "Exception is " << Error.what() << endl;
	}
	xml_node<>* rootNode = tree.first_node();
	if( 0 != strcmp(rootNode->name(), "Parameters") )
	{
		cout<<"Root node is not 'Parameters'! Will now halt!"<<endl;
		exit(2);
	}

	//Now read all parameters if present
	if( 0 != rootNode->first_node("NumberOfParticles") )
	{
		NumberOfParticles = std::stoi(rootNode->first_node("NumberOfParticles")->value());
	}
	if( 0 != rootNode->first_node("AnalysysTimestep") )
	{
		AnalysysTimestep = std::stoi(rootNode->first_node("AnalysysTimestep")->value());
	}
	if( 0 != rootNode->first_node("NumberOfTimesteps") )
	{
		NumberOfTimesteps = std::stoi(rootNode->first_node("NumberOfTimesteps")->value());
	}
	if( 0 != rootNode->first_node("DiscreteSamplingRate") )
	{
		DiscreteSamplingRate = std::stoi(rootNode->first_node("DiscreteSamplingRate")->value());
	}
	if( 0 != rootNode->first_node("DomainBorderSize") )
	{
		DomainBorderSize = std::stoi(rootNode->first_node("DomainBorderSize")->value());
	}
	if( 0 != rootNode->first_node("MicrostructureGenMode") )
	{
		MicrostructureGenMode = (E_MICROSTRUCTURE_GEN_MODE)std::stoi(rootNode->first_node("MicrostructureGenMode")->value());
		if(MicrostructureGenMode >= E_INVALID_VALUE)
			MicrostructureGenMode = E_INVALID_VALUE;
	}
	if( 0 != rootNode->first_node("ReadFromFilename") )
	{
		ReadFromFilename = rootNode->first_node("ReadFromFilename")->value();
	}
	if( 0 != rootNode->first_node("HAGB") )
	{
		HAGB = std::stod(rootNode->first_node("HAGB")->value());
	}
	if( 0 != rootNode->first_node("UseMobilityFactor") )
	{
		UseMobilityFactor = (bool)std::stoul(rootNode->first_node("UseMobilityFactor")->value());
	}
	if( 0 != rootNode->first_node("IsIsotropicNetwork") )
	{
		IsIsotropicNetwork = (bool)std::stoul(rootNode->first_node("IsIsotropicNetwork")->value());
	}
	if( 0 != rootNode->first_node("UseTexture") )
	{
		UseTexture = (bool)std::stoul(rootNode->first_node("UseTexture")->value());
	}
	if( 0 != rootNode->first_node("ExecuteInParallel") )
	{
		ExecuteInParallel = (bool)std::stoul(rootNode->first_node("ExecuteInParallel")->value());
	}
	file.close();
}
