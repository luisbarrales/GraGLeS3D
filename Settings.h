#ifndef __SETTINGS_H__
#define __SETTINGS_H__

#include <string>
#include "rapidxml.hpp"


enum E_MICROSTRUCTURE_GEN_MODE
{
	E_READ_FROM_FILE,
	E_GENERATE_WITH_VORONOY,
	E_READ_VERTEX,
	E_GENERATE_TESTCASE,
	E_INVALID_VAL
};

enum E_CONVOLUTION_MODE
{
	E_LAPLACE,
	E_LAPLACE_RITCHARDSON,
	E_GAUSSIAN,
	E_INVALID_VALUE
};

class Settings
{
	public:
	static unsigned long StartTime;
	static unsigned long NumberOfParticles;
	static unsigned long NumberOfPointsPerGrain;
	static unsigned long NumberOfTimesteps;
	static unsigned long AnalysisTimestep;
	static unsigned long DiscreteSamplingRate;
	static unsigned long DomainBorderSize;
	static E_MICROSTRUCTURE_GEN_MODE MicrostructureGenMode;
	static std::string ReadFromFilename;
	static double HAGB;
	static double TriplePointDrag;
	static bool UseMobilityFactor;
	static bool IsIsotropicNetwork;
	static bool UseTexture;
	static bool ExecuteInParallel;
	static bool GridCoarsement;
	static bool ResearchMode;
	static double GridCoarsementGradient;
	static unsigned long MaximumNumberOfThreads;
	static E_CONVOLUTION_MODE ConvolutionMode;

	static void initializeParameters(std::string filename = "");
	static rapidxml::xml_node<>* generateXMLParametersNode(rapidxml::xml_document<>* root, const char* filename, int loop, int grains);
};

#endif	//__SETTINGS_H__
