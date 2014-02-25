#include "parallelHandler.h"
#include "Settings.h"
#include "omp.h"
void parallelHandler::run_sim()
{
	initParallelEnvironment();

	find_neighbors();

	for(loop=0; loop <= Settings::NumberOfTimesteps; loop++){
		//Switch Distance Buffers
		for (int i = 1; i < grains.size(); i++){
			if(grains[i]==NULL)
				continue;
			grains[i]->switchInNOut();
		}
		convolution();

#pragma omp parallel
{
		//Switch Distance Buffers
#pragma omp for schedule(dynamic)
		for (int i = 1; i < grains.size(); i++){
			if(grains[i]==NULL)
				continue;
			grains[i]->switchInNOut();
		}
		//Update Second Order Neighbour
#pragma omp for
		for (int i = 1; i < grains.size(); i++){
			if(grains[i]==NULL)
				continue;
			grains[i]->add_n2o_2();
		}
		//Comparison box
#pragma omp for
		for (int i = 1; i < grains.size(); i++){
			if(grains[i]==NULL)
				continue;
			grains[i]->comparison();
		}
		//Switch Distance Buffers
#pragma omp for
		for (int i = 1; i < grains.size(); i++){
			if(grains[i]==NULL)
				continue;
			grains[i]->switchInNOut();
		}
		//Level Set
#pragma omp for
		for (int i = 1; i < grains.size(); i++){
			if(grains[i]==NULL)
				continue;
			grains[i]->find_contour();
		}
		//redistancing
#pragma omp single
		currentNrGrains=0;
#pragma omp for
		for (int i = 1; i < grains.size(); i++){
			if(grains[i]==NULL)
				continue;
			grains[i]->redist_box();
#pragma omp atomic
			currentNrGrains++;
		}
}
		if ( (loop % int(Settings::AnalysysTimestep)) == 0 || loop == Settings::NumberOfTimesteps ) {
			saveAllContourEnergies();
			save_texture();
			saveMicrostructure();
		}

	}
	cout << "Simulation complete." << endl;
}

void parallelHandler::initParallelEnvironment()
{
	if (Settings::MaximumNumberOfThreads == 0)
	{
		Settings::MaximumNumberOfThreads = omp_get_max_threads();
	}
	omp_set_num_threads(Settings::MaximumNumberOfThreads);


}
