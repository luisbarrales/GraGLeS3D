#include "ggLS.h"
#include "Settings.h"
#include "parallelHandler.h"
using namespace voro;
using namespace std;

int main(int argc,char *argv[]) {

	if(argc > 1)
		Settings::initializeParameters(argv[1]);
	else
		Settings::initializeParameters();

	grainhdl* my_sim = NULL;

	if(Settings::ExecuteInParallel)
	{
		my_sim = new parallelHandler() ;
	}
	else
	{
		my_sim = new grainhdl();
	}
	
	my_sim->setSimulationParameter();
	
	my_sim->saveMicrostructure();

	clock_t begin = clock();

<<<<<<< HEAD
        my_sim->run_sim();

        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout << "elapsed secs for main loop:" << elapsed_secs << endl;

	
=======
	my_sim->run_sim();

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "elapsed secs for main loop:" << elapsed_secs << endl;

>>>>>>> e791f3603f48051d2aa8871d93a5d4ebdda5cc37
	my_sim->save_sim();
	
	my_sim->clear_mem();
	
	delete my_sim;

}
