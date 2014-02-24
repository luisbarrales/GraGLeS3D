#ifndef __SETTINGS_H__
#define __SETTINGS_H__

#include <string>

enum E_MICROSTRUCTURE_GEN_MODE
{
	E_READ_FROM_FILE,
	E_GENERATE_WITH_VORONOY,
	E_READMICROSTRUCTURE,
	E_INVALID_VALUE
};

class Settings
{
	public:
	static unsigned int NumberOfParticles;
	static unsigned int NumberOfTimesteps;
	static unsigned int AnalysysTimestep;
	static unsigned int DiscreteSamplingRate;
	static unsigned int DomainBorderSize;
	static E_MICROSTRUCTURE_GEN_MODE MicrostructureGenMode;
	static std::string ReadFromFilename;
	static double HAGB;
	static bool UseMobilityFactor;
	static bool IsIsotropicNetwork;
	static bool UseTexture;
	static bool ExecuteInParallel;


	void static initializeParameters(std::string filename = "");

};

#endif	//__SETTINGS_H__
