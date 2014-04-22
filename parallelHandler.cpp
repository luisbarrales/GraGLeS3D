#include "parallelHandler.h"
#include "Settings.h"
#include "omp.h"
void parallelHandler::run_sim()
{
	simulationTime =0;
	find_neighbors();

	for(loop=0; loop <= Settings::NumberOfTimesteps; loop++){
		//Switch Distance Buffers
	
#pragma omp parallel
{
  

if (sqrt(currentNrGrains)*Settings::NumberOfPointsPerGrain/realDomainSize < 0.95 && loop!=0){
	  double shrink = 1-sqrt(currentNrGrains)*Settings::NumberOfPointsPerGrain/realDomainSize;
	  #pragma omp for  
	    for (int i = 1; i < grains.size(); i++){
		if(grains[i]==NULL)
			continue;
		  grains[i]->resizeGrid(shrink);
	    }	
	 #pragma omp barrier
	 #pragma omp single
	  {
	    cout << "RESIZED IN TIMESTEP: " << loop << "   with: " <<  sqrt(currentNrGrains)*Settings::NumberOfPointsPerGrain/realDomainSize<< endl;
	    cout << "OLD : " << ngridpoints << " new: " <<int(ngridpoints*(1-shrink)+1)<< endl;  
	    cout << "old dt " << dt << endl;
	    cout << "Number of remaining grains: " << currentNrGrains;
	    cout << "shrink " << shrink<< endl;		
	    realDomainSize = realDomainSize * (1-shrink)+1; 
	    ngridpoints = realDomainSize+2*grid_blowup; 
	    h = 1.0/realDomainSize;
	    dt = 1.0/double(realDomainSize*realDomainSize);
// 	    cout << "new dt " << dt << endl;
// 	    char bufferwait;
// 	    cin >> bufferwait;
	    
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
		simulationTime += dt;
}


	}

	cout << "Simulation complete." << endl;
	cout << "Simulation Time: " << simulationTime<< endl;
}

