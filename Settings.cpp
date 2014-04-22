#include "Settings.h"
#include "rapidxml.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>

using namespace std;
using namespace rapidxml;

//Initializing the static setting variables
unsigned long Settings::StartTime=0;
unsigned long Settings::NumberOfParticles = 0;
unsigned long Settings::NumberOfPointsPerGrain = 0;
unsigned long Settings::NumberOfTimesteps = 0;
unsigned long Settings::AnalysysTimestep = 0;
unsigned long Settings::DiscreteSamplingRate = 0;
unsigned long Settings::DomainBorderSize = 0;
E_MICROSTRUCTURE_GEN_MODE Settings::MicrostructureGenMode = E_INVALID_VALUE;
string Settings::ReadFromFilename;
double Settings::HAGB = 0.0;
double Settings::TriplePointDrag=0.0;
bool Settings::UseMobilityFactor = false;
bool Settings::IsIsotropicNetwork = false;
bool Settings::UseTexture = false;
bool Settings::ExecuteInParallel = false;
bool Settings::GridCorasment = false;
unsigned long Settings::MaximumNumberOfThreads = 0;

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
	catch (parse_error& Error)
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
	if( 0 != rootNode->first_node("StartTime") )
	{
		StartTime = std::stoul(rootNode->first_node("StartTime")->value());
	}
	if( 0 != rootNode->first_node("NumberOfParticles") )
	{
		NumberOfParticles = std::stoul(rootNode->first_node("NumberOfParticles")->value());
	}
	if( 0 != rootNode->first_node("NumberOfPointsPerGrain") )
	{
		NumberOfPointsPerGrain = std::stoul(rootNode->first_node("NumberOfPointsPerGrain")->value());
	}
	if( 0 != rootNode->first_node("AnalysysTimestep") )
	{
		AnalysysTimestep = std::stoul(rootNode->first_node("AnalysysTimestep")->value());
	}
	if( 0 != rootNode->first_node("NumberOfTimesteps") )
	{
		NumberOfTimesteps = std::stoul(rootNode->first_node("NumberOfTimesteps")->value());
	}
	if( 0 != rootNode->first_node("DiscreteSamplingRate") )
	{
		DiscreteSamplingRate = std::stoul(rootNode->first_node("DiscreteSamplingRate")->value());
	}
	if( 0 != rootNode->first_node("DomainBorderSize") )
	{
		DomainBorderSize = std::stoul(rootNode->first_node("DomainBorderSize")->value());
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
	if( 0 != rootNode->first_node("TriplePointDrag") )
	{
		TriplePointDrag = std::stod(rootNode->first_node("TriplePointDrag")->value());
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
	if( 0 != rootNode->first_node("MaximumNumberOfThreads") )
	{
	MaximumNumberOfThreads = std::stoul(rootNode->first_node("MaximumNumberOfThreads")->value());

	}
	if( 0 != rootNode->first_node("GridCorasment") )
	{
		GridCorasment = std::stoul(rootNode->first_node("GridCorasment")->value());
	}
	file.close();
}

#define PUSH_PARAM(param_name) 	\
		temp_string.str("");	\
		temp_string << param_name ;	\
		params->append_node(root->allocate_node(node_element,	\
				root->allocate_string(#param_name),	\
				root->allocate_string(temp_string.str().c_str()) ));

#define PUSH_VALUE(param_name, value) 	\
		temp_string.str("");	\
		temp_string << value ;	\
		params->append_node(root->allocate_node(node_element,	\
				root->allocate_string(#param_name),	\
				root->allocate_string(temp_string.str().c_str()) ));


xml_node<>* Settings::generateXMLParametersNode(xml_document<>* root, const char* filename, int loop, int grains)
{
	xml_node<>* params = root->allocate_node(node_element, "Parameters", "");
	stringstream temp_string;

	PUSH_VALUE(StartTime, loop);
	PUSH_VALUE(NumberOfParticles, grains);
	PUSH_PARAM(NumberOfPointsPerGrain);
	PUSH_PARAM(AnalysysTimestep);
	PUSH_PARAM(NumberOfTimesteps);
	PUSH_PARAM(DiscreteSamplingRate);
	PUSH_PARAM(DomainBorderSize);
	PUSH_VALUE(MicrostructureGenMode, 3);
	PUSH_VALUE(ReadFromFilename, filename);
	PUSH_PARAM(HAGB);
	PUSH_PARAM(TriplePointDrag);
	PUSH_PARAM(UseMobilityFactor);
	PUSH_PARAM(IsIsotropicNetwork);
	PUSH_PARAM(UseTexture);
	PUSH_PARAM(ExecuteInParallel);
	PUSH_PARAM(MaximumNumberOfThreads);
	PUSH_PARAM(GridCorasment);
	return params;
}
#undef PUSH_PARAM
#undef PUSH_VALUE
