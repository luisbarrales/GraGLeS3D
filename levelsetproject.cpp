/*----------------------------------------------------------------------------------------------------------------------*/

/*		Programm zur Simulation der Kornvergroesserung----------------------------------------------------------------------*/

/*		Ziel: Parallelisierung von bekannten Methoden---------------------------------------------------------------------*/

/*		Aufbau: Der Algorithmus wird in vier Teile gegliedert, sodass einzelne
 Module ausgetauscht werden koennen:
 1. Initialisierungsschritt (siehe "initialisation.h")
 2. Wachstumsschritt (surfacemotion.h)
 3. Korrekturschritt (correction.h)
 4. Redistanzierungsschritt (redistancing.h)
 
 Die Schritte 2 und 3 stellen den Kern des Algorithmus dar, sie werden je Zeitschritt einmal ausgefuehrt.
 Das "Redistancing" ist aus Stabilitaetsgruenden notwendig - jedoch nicht in jedem Zeitschritt.
 
 Beschreibung der Kornoberflaechen: Hierzu wird ein spezieller Level Set Ansatz gewaehlt. Wir identifizieren
 die Oberflaeche eines Korns ueber die Nullstellenmenge einer Hilfsfunktion \Phi[i,j] = dist(x_i,j , \Gamma).
 Diese Funktion wird auf einem Aequidistanten Gitter ausgewertet und traegt die Werte des kuerzesten
 vorzeichenbehafteten Abstandes zur Kornoberflaeche. Die Steigung ist in Normalenrichtung zur Oberflaeche in
 jedem Punkt 1. Wir werten dies Funktion nur auf eine eps-Schlauch um die zu beschreibende Oberflaeche aus und
 kâˆšâˆ‚nnen so mehrere/viele Koerner in einem Gitter speichern.----------------------------------------------------------*/

/*		Stand: März 2013 --------------------------------------------------------------------------------------------------*/

/*		Autor: C. Miessen--------------------------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------------------------------------------*/



#include "levelsetproject.h"



double kernel(double dt,int m, int i ,int j){
    double erg;
    double nominator = 4.0 * PI *dt;
    double dx = 1/ (double) m;
	double abs_x_sq= ((i*dx)*(i*dx)) + ((j*dx)*(j*dx));
    erg = 1/nominator * exp(-abs(abs_x_sq) /(4.0*dt));
    return(erg);
}



using namespace voro;

int main() {
	/***************/
    // Init
    /***************/
    char buffer;
    const int particles = int(PARTICLES);
    double dt = 1.0/double(M*M);
	double dt_e = dt/2 ;
    
    double x,y,z,rx,ry,rz;
    int m = M, current_cell, cell_id;
    
    int cell_order[particles]; // stores the order of computed cells (vl.pid())
    
    double *part_pos; // stores the centroids of the cells ; access by (3*Id, 3*ID +1, 3*ID +2)
    part_pos = (double*) calloc (3*particles,sizeof(double));
    stringstream filename, plotfiles;
    
    std::list<matrix> domains, domains_copy;
    std::list<matrix>::iterator it, itc;
    std::vector<LSbox*>::iterator it2;
    std::vector<LSbox*>::iterator it2c;
    
    voronoicell_neighbor c;
    
    bool randbedingung = false; // fÂ¸r false ist der container halb offen?! d.h. gitterwert mit 1 werden keinem partikel zugeordnet
    if (randbedingung == false) m = int(M)-1;
    
    const double h = 1.0/double(m); // there are m-2*gridbloup grid points in each direction in the domain
    const int grid_blowup = int(((double)DELTA / h)+1); // number of grid points to extract the domain at each boundary
    
    int resized_m = m + (2*grid_blowup); //resize m
    LSbox **ID; // array to asign a cell id to each grid point
    ID = new LSbox*[2];
	ID[0] = new LSbox[resized_m*resized_m];
	ID[1] = new LSbox[resized_m*resized_m];
    
    double  (*fp)(double,int, int, int); // function pointer
    fp = &kernel;
    
    container con(0,1,0,1,0,1,5,5,5,randbedingung,randbedingung,randbedingung,2);
    c_loop_all vl(con);
    
    
    //program options:
    cout << endl << "******* PROGRAM OPTIONS: *******" << endl << endl;
    cout << "Number of Grains: " << PARTICLES << endl;
    if (!EULER ) cout << "FFT "<< endl; else cout << "EULER FORWARD "<< endl;
    cout << "simulated Timesteps: " << TIMESTEPS << endl;
    cout << "Timestepwidth " << dt << endl;
    cout << "Number of Gridpoints: " << M << endl << endl;
    
    cout << endl << "******* start simulation: *******" << endl << endl;
    
    /******************************/
    // Randomly add particles into the container
    /******************************/
    
    for(int i=0;i<particles;i++) {
        x=utils::rnd();
        y=utils::rnd();
        z=0;
        con.put(i,x,y,z);
    }
    
    // generate random or deterministic surface tension coefficients
	double *ST;
	ST = (double*) calloc (PARTICLES*PARTICLES,sizeof(double));
	const double MIN = 0.6;
	const double MAX = 1.0;

	for(int i=0; i < PARTICLES; i++)
		for(int j=0; j < PARTICLES; j++){			
			double zahl=(double)(rand() / (((double)RAND_MAX+1)/ (double)(MAX-MIN)))+MIN;
			ST[i+(PARTICLES*j)] = zahl;
			ST[j+(PARTICLES*i)] = zahl;
		}
		
	/*utils::print_2dim_array(ST, PARTICLES, PARTICLES);
	
// 	*/
// 	cin >> buffer;

	
    
    // 	Output the Voronoi cells to a file, in the gnuplot format
    con.draw_cells_gnuplot("random_points_v.gnu");
    
    // find cell information fpr each grid point
    for(int i=0; i < m; i++) for(int j=0; j < m; j++){
        x=double(i*h);y=double(j*h); // only point within the domain
        if(con.find_voronoi_cell(x,y,z,rx,ry,rz,cell_id)){
            ID[(i+grid_blowup)*resized_m + (j+grid_blowup)]=cell_id;
            part_pos[3*cell_id]=rx;
            part_pos[3*cell_id+1]=ry;
            part_pos[3*cell_id+2]=rz;
        }
        else fprintf(stderr,"# find_voronoi_cell error for %g %g 0\n",x,y);
    }
    utils::save_2dim_array( ID, resized_m, resized_m , "IDmatrix.gnu");
    
	con.draw_cells_gnuplot("particles.gnu");
	
	
/*********************************************************************************/
// Initialisation: Voronoizellen -> Box
/*********************************************************************************/
    
	int i=0;
	// iteration over all cells in the container con:
	if(vl.start()) 
	  do {
	    // compute the current cell, taken out of the container
	    con.compute_cell(c,vl);
	    cell_order[particles-1-i]=vl.pid();
	    
	    // create a new Box for the current cell
	    LSbox* newBox = new LSbox(vl.pid(), c, part_pos, grid_blowup, h);
	    // find domain for new box
	    bool foundDomain = false;
	    cout << "trying to add a box -- ";
	    if (!domains.empty())
	      for (it = domains.begin(); it !=domains.end(); it++){
		  foundDomain = (*it).addBox(newBox);
		  if (foundDomain) break;
	      }
	  
	  
	    if (!foundDomain) {
		// create domain
		cout << "failed: creating new domain" << endl;
			    domains.emplace_back(resized_m,resized_m, i,INTERIMVAL);  // = push_back
			    i++;
		domains.back().addBox(newBox);            
	    } else cout << "success" << endl;
		    
	    // calculate distances
	    newBox->distancefunction(c, ID, part_pos, grid_blowup, h);        
	  
	  } while(vl.inc());

/*********************************************************************************/

    int j;
    for (it = domains.begin(), j = 0; it !=domains.end(); it++, j++){
        filename.str(std::string());
        filename << "Distanzmatrix";
        vector<LSbox*> grains = (*it).getBoxList();
        vector<LSbox*>::iterator it2;
        for (it2 = grains.begin(); it2 != grains.end(); it2++) {
            filename << "_" <<(*it2)->getID();
        }
        filename << ".gnu";    
		
        if (SAFEFILES){
        	cout << filename.str() << endl << endl;        
			(*it).save_matrix(filename.str().c_str());
		}
    }
/*********************************************************************************/


// Determine the box' neighbors 

for (it = domains.begin(); it !=domains.end(); it++){
	vector<LSbox*> grains = (*it).getBoxList();
	vector<LSbox*> grainstoComp;

	for (itc=it;itc!=domains.end();itc++){
		grainstoComp = itc->getBoxList();
		for (it2=grains.begin();it2!=grains.end();it2++){
			for(it2c=grainstoComp.begin();it2c!=grainstoComp.end();it2c++){
				if ((*it2)->getID()!=(*it2c)->getID()){
					if((*it2)->checkIntersect(*it2c)){
						(*it2) ->	neighbors.push_back(*it2c);
						(*it2c)->	neighbors.push_back(*it2);
						cout <<"Grain: "<< (*(*it2)).getID() << " with Grain: " << (*(*it2c)).getID()<<endl;
					}
				}
			}
		}
	}
}


	
/*********************************************************************************/
// MAIN LOOP
/*********************************************************************************/

vector<int> nr_grains(TIMESTEPS);

for(int loop=0; loop <= TIMESTEPS; loop++){

	stringstream plotfiles;
	plotfiles.str(std::string());  
	
	/*********************************************************************************/
	// Convolution simulates grain growth
	/*********************************************************************************/
	/*if(!EULER){
		for (it = domains.begin(); it !=domains.end(); it++){	
			(*it).convolution(dt);
		}
	}
	else {
		for (it = domains.begin(); it !=domains.end(); it++){	
		 (*it).euler(dt_e,h);
		}
	}*/
	// ACHTUNG hier kopieren wir die ganze LISTE!
	domains_copy=domains;
	
	for (it = domains.begin(); it !=domains.end(); it++){	
		(*it).convolution(dt);
	}
	for (it = domains.begin(); it !=domains.end(); it++){
		// Output			
		if ((loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS){
			filename.str(std::string());
			filename << "Convolutedmatrix_";
			vector<LSbox*> grains = (*it).getBoxList();
			vector<LSbox*>::iterator it2;
			filename << "T"<<loop;
			for (it2 = grains.begin(); it2 != grains.end(); it2++) {
				filename << "_"<<(*it2)->getID();
			}
			filename << ".gnu";
			if (SAFEFILES) {
				(*it).save_matrix(filename.str().c_str());
				cout << filename.str() << endl << endl;
			}
		}
	}
	

	
	/*********************************************************************************/
	// Comparison Step: step 0.5 *(A_k(x) - max A_i(x) | i!=k)
	/*********************************************************************************/
	// Create a list for storing the new distances after comparison

	    /*********************************************************************************/
	    // Comparison Step with domains
	    /*********************************************************************************/
	    if(DOMAINCOMPARISON==true){
		  domains_copy = domains;
		for (it = domains.begin(); it != domains.end(); it++){		  
		  (*it).comparison(domains_copy, grid_blowup);
			if ((loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS){
			  vector<LSbox*>::iterator it2;
			  vector<LSbox*> grains;    
			  
			  filename.str(std::string());
			  filename << "Comparedmatrix_"<< "T"<<loop;;
			  grains = (*it).getBoxList();
			 
			  for (it2 = grains.begin(); it2 != grains.end(); it2++) {
			    filename << "_"<<(*it2)->getID();
			  }
			  
			  filename << ".gnu";
			  if (SAFEFILES) {
			    (*it).save_matrix(filename.str().c_str());
			    cout << filename.str() << endl << endl;
			  }
			}
		}
	    }
	    
	    else {
	      /*********************************************************************************/
	      // Comparison Step on boxes
	      /*********************************************************************************/	
			for (it = domains.begin(); it != domains.end(); it++){		     
				grains = (*it).getBoxList();
				for (it_domain = domains_copy.begin(); it_domain != domains_copy.end(); it_domain++){
					for (it2 = grains.begin(); it2 != grains.end(); it2++){
						(*it2).comparison((*it_domain));
					}
				}
				for (it2 = grains.begin(); it2 != grains.end(); it2++){
					(*it2).copy_distances_to_domain();
				}
			}
	    }
	    
	    
/*************** 			Auskommentierter Code von Comparison ************/
/*			bool exist = false;
			exist = (*it).comparison(domains_copy, grid_blowup);
		
	// 		if (exist == false) {
	// 			cout << "now we delete domain: "<< (*itc).get_id() << endl << endl;
	// 			itc = compared_domains.erase(itc);
	// 			itc--;
	// 			it = domains.erase(it);
	// 			it--;
	// 		}
	// 		else {
*/

			

//*********************************************************************************/
// Redistancing Step
/*********************************************************************************/


	/****************************************************/
	// Nullstellenverfolgung:
	// Speichert Nullstellen als NNZ-OBjekt in jeder Box
	// Testet Boxen auf Konlikte
	// Verschiebt "Konfliktboxen" in andere Domain
	
	
	// Ist dieser Schritt in jedem Zeitschritt notwenig????
	
	vector<LSbox*> buffer;
	for (it = domains.begin(); it != domains.end(); it++) {
		bool exist=true;
		char buffer1;
		//check domain it for intersecting grains
		exist = (*it).grainCheck(h, grid_blowup, buffer); // h Gitterabstand
		if (!exist){
			cout << (*it).get_id() <<"domain leer" << endl;
// 			domains.erase(it); it--;
// 			cin >> buffer1;
		}
		
	}
	// check if buffer is empty
	while(!buffer.empty()){
		cout << "created a new domain" << endl;
		domains.emplace_back(resized_m,resized_m, i,INTERIMVAL);
		domains.back().grainCheck(h, grid_blowup, buffer);
	}
	/****************************************************/
	
	
	

	
	
	
	for (i=0, it = domains.begin(); it != domains.end(); it++, i++) {
		//Nullstellenverfolgung:
		//cout << "Rechne Redistancing auf Boxen der Domain: " << (*it).get_id() << endl << endl;
		// zugriff auf Boxen über die Domain "it"
		// Intern können verschiedenRedistancing Routinen verwendet werden

		(*it).redistancing_2(h, grid_blowup);		
		(*it).redistancing_2(h, grid_blowup);
		//doppelt?
		// (*it).clear_domain(INTERIMVAL);
		// (*it).redistancing_for_all_boxes(h, grid_blowup);
		nr_grains[loop]+=(*it).get_nr_of_grains();
		
		if ( (loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS){
				filename.str(std::string());
				filename << "Redistanced_matrix_";
				filename << it->get_id() << "_";
				vector<LSbox*> grains = (*it).getBoxList();
				vector<LSbox*>::iterator it2;
				filename << "T"<<loop;
				for (it2 = grains.begin(); it2 != grains.end(); it2++) {
						filename <<"_"<< (*it2)->getID();
				}
				filename << ".gnu";
				if (SAFEFILES) {
					(*it).save_matrix(filename.str().c_str());
					cout << filename.str() << endl << endl;
				}
				plotfiles << " \""<<filename.str();
				plotfiles << "\" matrix w l";
				if(i!=(length-1)) plotfiles << ",";
			}
		}
		if ( (loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS){
			if (PLOTGNU) {
				filename.str(std::string());
				filename << "GrainNetwork" << "_"<< loop << ".gnu";
				utils::plotGnu(filename.str().c_str(), plotfiles.str().c_str());
			}			
			if (IMAGEOUT) {
				int imgnum = (loop/PRINTSTEP);
				filename.str(std::string());
				filename << "GrainNetwork";
				if (imgnum < 100) filename << "0";
				if (imgnum < 10) filename << "0";
				filename << imgnum << ".png";
				utils::plotGnuPNG(filename.str().c_str(), plotfiles.str().c_str());
			}
		}
			
	cout << "Timestep: "<< loop << " complete" << endl;
	cout << "Number of remaining grains: "<< nr_grains[loop] << endl << endl;
   
	
	/************************Auskommentiert Redistancing ***********************/
   /*
	// 	int length = domains.size();
	// 	omp_set_dynamic(0);
	// 	omp_set_num_threads(length);
	// 	#pragma omp parallel
	// 	#pragma omp single
	// 	{
	// 	for (auto it = domains.begin(); it != domains.end(); it++) {
	// 		#pragma omp task firstprivate(it)
	// 		{
	// 		(*it).redistancing_2(h, grid_blowup);
	// 		cout << "I compute domain "<< (*it).get_id() << " in " << omp_get_thread_num() << " --- "<< omp_get_num_threads() << endl;;
	// 		}
	// 		}
	// 	#pragma omp taskwait
	// 	}
   */

}

/*******************************************************************************************/
// end of simulation
/*******************************************************************************************/
	ofstream myfile;
	myfile.open ("kinetics.txt");
	for(int i=0; i< TIMESTEPS; i++)
		myfile << nr_grains[i] << "\t";
	myfile.close();

	utils::PNGtoGIF("test.mp4");
	cout << "number of distanzmatrices: "<< domains.size() << endl;
//	//utils::print_2dim_array( ID, m, m );
	
 	free (ID);    
	return 0;
}
