#include "ggLS.h"
#include "Settings.h"
#include "parallelHandler.h"
using namespace voro;

int main() {
	
	Settings::initializeParameters();

	grainhdl* my_sim = NULL;

	if(Settings::ExecuteInParallel)
	{
		my_sim = new parallelHandler();
	}
	else
	{
		my_sim = new grainhdl();
	}
	
	my_sim->setSimulationParameter();
	
	my_sim->run_sim();
	
	my_sim->save_sim();
	
	my_sim->clear_mem();
	
	delete my_sim;
}
