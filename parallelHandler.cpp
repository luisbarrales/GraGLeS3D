#include "parallelHandler.h"
#include "Settings.h"
#include "omp.h"
void parallelHandler::run_sim()
{
	find_neighbors();

	for(loop=Settings::StartTime; loop <= Settings::StartTime + Settings::NumberOfTimesteps; loop++){
		//Switch Distance Buffers
#pragma omp parallel
{
#pragma omp for
		for (int i = 1; i < grains.size(); i++){
			if(grains[i]==NULL)
				continue;
			grains[i]->switchInNOut();
		}
#pragma omp for
		for (int i = 1; i < grains.size(); i++){
			if(grains[i]==NULL)
				continue;
			grains[i]->convolution(m_ThreadMemPool[omp_get_thread_num()]);
		}

		//Switch Distance Buffers
#pragma omp for
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
			grains[i]->comparison(m_ThreadMemPool[omp_get_thread_num()]);
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
#pragma omp single
{
		if ( (loop % int(Settings::AnalysysTimestep)) == 0 || loop == Settings::NumberOfTimesteps ) {
			saveAllContourEnergies();
			save_texture();
			saveMicrostructure();
		}
}

	}
	cout << "Simulation complete." << endl;
}

