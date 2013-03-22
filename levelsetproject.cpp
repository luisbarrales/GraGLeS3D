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

/*		Stand: Juni 2012 --------------------------------------------------------------------------------------------------*/

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
    
    const int particles = int(PARTICLES);
    const double dt = 1.0/double(M*M);
    
    double x,y,z,rx,ry,rz;
    int m = M, current_cell, cell_id;
    
    int cell_order[particles]; // stores the order of computed cells (vl.pid())
    
    double *part_pos; // stores the centroids of the cells ; access by (3*Id, 3*ID +1, 3*ID +2)
    part_pos = (double*) calloc (3*particles,sizeof(double));
    stringstream filename;
    
    std::list<matrix> domains, domains_copy;
    std::list<matrix>::iterator it, itc;
    
    voronoicell_neighbor c;
    
    bool randbedingung = false; // fÂ¸r false ist der container halb offen?! d.h. gitterwert mit 1 werden keinem partikel zugeordnet
    if (randbedingung == false) m = int(M)-1;
    
    const double h = 1.0/double(m); // there are m-2*gridbloup grid points in each direction in the domain
    const int grid_blowup = int(((double)DELTA / h)+1); // number of grid points to extract the domain at each boundary
    
    int resized_m = m + (2*grid_blowup); //resize m
    int *ID; // array to asign a cell id to each grid point
    ID = (int*) calloc ((resized_m)*(resized_m),sizeof(int));
    
    double  (*fp)(double,int, int, int); // function pointer
    fp = &kernel;
    
    container con(0,1,0,1,0,1,5,5,5,randbedingung,randbedingung,randbedingung,2);
    c_loop_all vl(con);
    
    
    //program options:
    cout << endl << "******* PROGRAM OPTIONS: *******" << endl << endl;
    cout << "Number of Grains: " << PARTICLES << endl;
    cout << "FFT: " << DISCRETE_CONVOLUTION << endl;
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
    
    
    /*********************************************************************************/
	// Initialisation: Voronoizellen -> Box
    /*********************************************************************************/
    
	int i=0;
	// iteration over all cells in the container con:
	if(vl.start()) do {
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
			domains.emplace_back(resized_m,resized_m, i,INTERIMVAL);
//             domains.push_back(matrix(resized_m,resized_m, i,INTERIMVAL));
			i++;
            domains.back().addBox(newBox);
            
        } else {
            cout << "success" << endl;
        }
        
        // calculate distances
        newBox->distancefunction(c, ID, part_pos, grid_blowup, h);
        
        
    } while(vl.inc());
    
    int j;
    for (it = domains.begin(), j = 0; it !=domains.end(); it++, j++){
        filename.str(std::string());
        filename << "Distanzmatrix_";
        vector<LSbox*> grains = (*it).getBoxList();
        vector<LSbox*>::iterator it2;
        for (it2 = grains.begin(); it2 != grains.end(); it2++) {
            filename << (*it2)->getID() << "_";
        }
        filename << "\b"<< ".gnu";
        
        cout << filename.str() << endl << endl;
        
        (*it).save_matrix(filename.str().c_str());
    }

/*********************************************************************************/
// MAIN LOOP
/*********************************************************************************/


for(int loop=0; loop <= TIMESTEPS; loop++){

	stringstream plotfiles;
	plotfiles.str(std::string());  


	/*********************************************************************************/
	// Convolution simulates grain growth
	/*********************************************************************************/

	for (it = domains.begin(); it !=domains.end(); it++){	
			
		(*it).convolution(dt);
			
		// Output			
		if ((loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS){
			filename.str(std::string());
			filename << "Convolutedmatrix_";
			vector<LSbox*> grains = (*it).getBoxList();
			vector<LSbox*>::iterator it2;
			filename << "T"<<loop<<"_";
			for (it2 = grains.begin(); it2 != grains.end(); it2++) {
				filename << (*it2)->getID() << "_";
			}
			filename << "\b"<< ".gnu";
			
			cout << filename.str() << endl << endl;
			
			(*it).save_matrix(filename.str().c_str());
		}
	}
	
	
	/*********************************************************************************/
	// Comparison Step: step 0.5 *(A_k(x) - max A_i(x) | i!=k)
	/*********************************************************************************/
	// Create a list for storing the new distances after comparison

		domains_copy=domains;
		for (it = domains.begin(), itc= domains_copy.begin(); it != domains.end(); it++, itc++){

			bool exist = false;
			exist = (*it).comparison(domains_copy, grid_blowup);
	// 		if (exist == false) {
	// 			cout << "now we delete domain: "<< (*itc).get_id() << endl << endl;
	// 			itc = compared_domains.erase(itc);
	// 			itc--;
	// 			it = domains.erase(it);
	// 			it--;
	// 		}
	// 		else {
				if ( (loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS){
					filename.str(std::string());
					filename << "Comparedmatrix_";
									
					vector<LSbox*> grains = (*it).getBoxList();
					vector<LSbox*>::iterator it2;
					filename << "T"<<loop<<"_";
					for (it2 = grains.begin(); it2 != grains.end(); it2++) {
						filename << (*it2)->getID() << "_";
					}
					filename << "\b"<< ".gnu";
					
					cout << filename.str() << endl << endl;
					
					(*it).save_matrix(filename.str().c_str());
				}
	}  
			
// 	cout << "Comparison complete" << endl << endl;;
	/*********************************************************************************/
	// Redistancing Step
	/*********************************************************************************/

	/*	// Slope-Field solution attempt
		// create slope-field
		for (int k = 0; k < resized_m; k++) {
			for (int l = 0; l < resized_m; l++) {
				double min1=99999, min2=99999; // just some random large numbers for first comparison
				int min1ID=-1, min2ID=-1;
				// find Minima in [k][l]
				std::list<matrix>::iterator it;
				for(it = distances.begin(); it != distances.end(); it++) {
					double val = abs((*it)[k][l]);

					if (val < min1) {
						min2 = min1; min2ID = min1ID;
						min1 = val; min1ID = (*it).get_id();
					}
					else if(val < min2) {
						min2 = val; min2ID = (*it).get_id();
					}
				}
				// assign slope
				if (min1ID != -1 && min2ID != -1) {
					slopeField[k][l] = borderSlopes[min1ID][min2ID];
				} else {
					slopeField[k][l] = 1;
				}
			}
		*/

	/****************************************************/
	// Nullstellenverfolgung:
	// Speichert Nullstellen als NNZ-OBjekt in jeder Box
	// Testet Boxen auf Konlikte
	// Verschiebt "Konfliktboxen" in andere Domain
	
	vector<LSbox*> buffer;
	for (it = domains.begin(); it != domains.end(); it++) {
		(*it).grainCheck(h, grid_blowup, buffer); // h Gitterabstand
	}
	
	/****************************************************/
	for (it = domains.begin(); it != domains.end(); it++) {
		//Nullstellenverfolgung:
// 		cout << "Rechne Redistancing auf Boxen der Domain: " << (*it).get_id() << endl << endl;
		
		// zugriff auf Boxen über die Domain "it"
		// Intern können verschiedenRedistancing Routinen verwendet werden
		(*it).redistancing_2(h, grid_blowup);
// 		(*it).clear_domain(INTERIMVAL);
// 		(*it).redistancing_for_all_boxes(h, grid_blowup);
		
		if ((loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS){
			filename.str(std::string());
			filename << "Redistanced_matrix_";
							
			vector<LSbox*> grains = (*it).getBoxList();
			vector<LSbox*>::iterator it2;
			filename << "T"<<loop<<"_";
			for (it2 = grains.begin(); it2 != grains.end(); it2++) {
				filename << (*it2)->getID() << "_";
			}
			
			filename << "\b"<< ".gnu";				
			cout << filename.str() << endl << endl;				
			(*it).save_matrix(filename.str().c_str());
		}
	}  
	cout << "loop: "<< loop << " complete" << endl;
}

/*******************************************************************************************/
// end of simulation
/*******************************************************************************************/

//    if (IMAGEOUT) {
//        // make gif
//        utils::PNGtoGIF("test.mp4");
//    }   
  
//    /*********************************************************************************/
//    /******************************************************************************/
	con.draw_cells_gnuplot("particles.gnu");
	cout << "number of distanzmatrices: "<< domains.size() << endl;
//	//utils::print_2dim_array( ID, m, m );
	
 	free (ID);    
	return 0;
}
