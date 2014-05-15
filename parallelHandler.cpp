#include "parallelHandler.h"
#include "Settings.h"
#include "omp.h"
void parallelHandler::run_sim()
{


#pragma omp parallel for
		for (int i = 1; i < grains.size(); i++){
			grains[i]->distancefunction();
		}


simulationTime =0;
find_neighbors();


for(loop=Settings::StartTime; loop <= Settings::StartTime + Settings::NumberOfTimesteps; loop++){
	//Switch Distance Buffers

#pragma omp parallel
{
if (sqrt(currentNrGrains)*Settings::NumberOfPointsPerGrain/realDomainSize < Settings::GridCoarsementGradient && loop!=0&& Settings::GridCoarsement){
	  double shrink = 1-sqrt(currentNrGrains)*Settings::NumberOfPointsPerGrain/realDomainSize;
	  #pragma omp for  
	    for (int i = 1; i < grains.size(); i++){
		if(grains[i]==NULL)
			continue;
		  grains[i]->resizeGrid(shrink);
	    }	
	 #pragma omp single
	  {
	    realDomainSize = realDomainSize * (1-shrink)+1; 
	    ngridpoints = realDomainSize+2*grid_blowup; 
	    h = 1.0/realDomainSize;
	    dt = 1.0/double(realDomainSize*realDomainSize);
	  }
    }

    else {
    #pragma omp for
		  for (int i = 1; i < grains.size(); i++){
			  if(grains[i]==NULL)
				  continue;
			  grains[i]->switchInNOut();  
		  }
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
			if(grains[i]->get_status() == false ) {
				delete grains[i];
				removeGrain(i);
			}
			else grains[i]->find_contour();
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


#pragma omp single
{
		if ( ((loop-Settings::StartTime) % int(Settings::AnalysisTimestep)) == 0 || loop == Settings::NumberOfTimesteps  ) {
			saveAllContourEnergies();
			save_texture();
			if(loop == Settings::NumberOfTimesteps) saveMicrostructure();
		}
		simulationTime += dt;
}

}
}

	cout << "Simulation complete." << endl;
	cout << "Simulation Time: " << simulationTime<< endl;
}

