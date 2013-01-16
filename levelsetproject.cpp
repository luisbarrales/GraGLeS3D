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
 k√∂nnen so mehrere/viele Koerner in einem Gitter speichern.----------------------------------------------------------*/

/*		Stand: Juni 2012 --------------------------------------------------------------------------------------------------*/

/*		Autor: C. Miessen--------------------------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------------------------------------------*/



#include "levelsetproject.h"



double kernel(double dt,int m, int i ,int j){
    double erg;
    double k = 2.0 * PI / m;
    double nsq = (double) m * (double) m;
    erg = exp((-2.0 * dt) * nsq * (2.0-cos(k*i)-cos(k*j)));
    return(erg);
}



using namespace voro;

int main() {

    /*********************************************************************************/
    // Initialisation
    /*********************************************************************************/
    
	const int particles = int(PARTICLES);
	const double dt = 1./(double(M*M));
    

	double x,y,z,rx,ry,rz;
	int m = M, current_cell, cell_id;
	
    int cell_order[particles]; // stores the order of computed cells (vl.pid())
    
    double** borderSlopes; // stores slopes between two grains in format [grain1][grain2]
    borderSlopes = (double**)calloc(particles,sizeof(double*));
    for(int i = 0; i < particles; i++) borderSlopes[i] = (double*)calloc(particles,sizeof(double));

	double *part_pos; // stores the centroids of the cells ; access by (3*Id, 3*ID +1, 3*ID +2)
	part_pos = (double*) calloc (3*particles,sizeof(double));
	stringstream filename;
	
	// Create a list for storing the distance to the nearest Voronoi volumes
	std::list<matrix> distances;
	std::list<matrix>::iterator it;
    
	std::list<matrix> compared_dist;
	std::list<matrix>::iterator itc;
	std::list<matrix>::iterator it_rm;
	voronoicell_neighbor c;

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block
    
	bool randbedingung = false; // f¸r false ist der container halb offen?! d.h. gitterwert mit 1 werden keinem partikel zugeordnet
	if (randbedingung == false) m = int(M)-1;
    
	const double h = 1.0/double(m); // there are m-2*gridbloup grid points in each direction in the domain
	const int grid_blowup = int(((double)DELTA / h)+1); // number of grid points to extract the domain at each boundary
    
	int resized_m = m + (2*grid_blowup); //resize m
	int *ID; // array to asign a cell id to each grid point
	ID = (int*) calloc ((resized_m)*(resized_m),sizeof(int));
    
    double** slopeField; // stores slopes for each point in the whole domain
    slopeField = (double**)calloc(resized_m, sizeof(double*));
    for (int i = 0; i < resized_m; i++) slopeField[i] = (double*)calloc(resized_m, sizeof(double));
    
	double  (*fp)(double,int, int, int); // function pointer
	fp = &kernel;
    
	container con(0,1,0,1,0,1,5,5,5,randbedingung,randbedingung,randbedingung,2);
	c_loop_all vl(con);



	//program options:
	cout << endl << "******* PROGRAM OPTIONS: *******" << endl << endl;
	cout << "Number of Grains: " << PARTICLES << endl;
	cout << "FFT: " << MODUS << endl;
	cout << "simulated Timesteps: " << TIMESTEPS << endl;
	cout << "Timestepwidth " << dt << endl;
	cout << "Number of Gridpoints: " << M << endl << endl;
		
	cout << endl << "******* start simulation: *******" << endl << endl;



	/*********************************************************************************/
	// Randomly add particles into the container
	/*********************************************************************************/
	/*

    for(int i=0;i<particles;i++) {
		x=utils::rnd();
		y=utils::rnd();
        // 		cout << x <<y <<endl;
		z=0;
		con.put(i,x,y,z);
	}*/
    
    
    
//     ****
//     Add particles for one triple-point
//     ***
    con.put(0,0.33,0.33,0);
    con.put(1,0.33,0.66,0);
    con.put(2,0.66,0.66,0);
    
    borderSlopes[0][1] = 1.2;
    borderSlopes[0][2] = 0.8;
    borderSlopes[1][0] = 1.2;
    borderSlopes[1][2] = 1.5;
    borderSlopes[2][0] = 0.8;
    borderSlopes[2][1] = 1.5;
    
// 	borderSlopes[0][1] = 1.;
//     borderSlopes[0][2] = 1.;
//     borderSlopes[1][0] = 1.;
//     borderSlopes[1][2] = 1.;
//     borderSlopes[2][0] = 1.;
//     borderSlopes[2][1] = 1.;
	
	
    /*********************************************************************************/

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
	// Initialisation: Voronoizellen -> Distancematrices
    /*********************************************************************************/
    
	int i=0;
	// iteration over all cells in the container con:
	if(vl.start()) do {
		// compute the current cell, taken out of the container
		con.compute_cell(c,vl);
		cell_order[particles-1-i]=vl.pid();

		vector<int> v; c.neighbors(v);
		utils::print_vector(v);
		// vl.pid(): holds current cell_id
		// c: current voronoicell
		// ID: matrix with Information grid point -> cell_id
		// part_pos: array with length 3 times particles,
		// 		holds information about the shift to the local coordinatesystem
		distances.push_front(matrix(resized_m,resized_m, vl.pid()));

		(*distances.begin()).distancefunction(c, ID, part_pos, grid_blowup, h);
		
        vector<int> cellNeighbors;
        c.neighbors(cellNeighbors);
        

        
//         for (vector<int>::iterator it = cellNeighbors.begin();
//              it != cellNeighbors.end(); ++it) {
//             if (*it >= 0) {
//                 borderSlopes[vl.pid()][*it] = 1; // ((rand() % 100) / 100.0) + 1;
//             }
//         }
        

        
		compared_dist.push_front(matrix(resized_m,resized_m, vl.pid()));
		// write the distancematrices to outputfile:
		filename.str(std::string());
		filename << "Distanzmatrix" << vl.pid() << ".gnu";
		cout << filename.str() << endl << endl;
		if (SAFEFILES)(*distances.begin()).save_matrix(filename.str().c_str());
	}
	while(vl.inc());
    
    
    /*******************************************************************************************/
    //// MAIN LOOP
    /*******************************************************************************************/
    /*******************************************************************************************/
    
	for(int loop=0; loop <= TIMESTEPS; loop++){
        
		stringstream plotfiles;
		plotfiles.str(std::string());
		
        /*********************************************************************************/
		// Convolution simulates grain growth
		/*********************************************************************************/
		
		for (it = distances.begin(); it !=distances.end(); it++){

		    
            if (DISCRETE_CONVOLUTION) (*it).discrete_convolution(dt, h, grid_blowup, fp);
			else	(*it).convolution(dt);
			
            // Output

			if ((loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS){
				filename.str(std::string());
				filename << "Convoluted_matrix" << (*it).get_id() << "_"<< loop << ".gnu";

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
        
        compared_dist=distances;
		for (it = distances.begin(), itc= compared_dist.begin(); itc != compared_dist.end(); it++, itc++){
			bool exist = false;
			exist = (*itc).comparison(distances, grid_blowup);
			if (exist == false) {
				cout << "now we delete domain: "<< (*itc).get_id() << endl << endl;;
				itc = compared_dist.erase(itc);
				itc--;
				it = distances.erase(it);
				it--;
			}
			else {
				if ( (loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS){
					filename.str(std::string());
					filename << "Compared_matrix" << (*itc).get_id() << "_"<< loop << ".gnu";
					
					if (SAFEFILES) {
						(*itc).save_matrix(filename.str().c_str());
						cout << filename.str() << endl << endl;
					}
				}
			}
		}
        
		
        
        /* Slope-Field solution attempt
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
         }
         //***
         */
        
        
        
		/*********************************************************************************/
		// Redistancing Step
		/*********************************************************************************/
		
		int length = compared_dist.size();
        for (i=0, itc= compared_dist.begin(); itc != compared_dist.end(); itc++, i++){
			(*itc).redistancing(h, grid_blowup, distances, borderSlopes, slopeField);
			if ( (loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS){
				filename.str(std::string());
				filename << "Redistanced_matrix" << (*itc).get_id() << "_"<< loop << ".gnu";
				
				if (SAFEFILES) {
					(*itc).save_matrix(filename.str().c_str());
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
		
		distances = compared_dist;
	}
    
    /*******************************************************************************************/
    // end of simulation
    /*******************************************************************************************/
    
    
    
    if (IMAGEOUT) {
        // make gif
        utils::PNGtoGIF("test.mp4");
    }
    
    /*********************************************************************************/
    /*********************************************************************************/
	
	con.draw_cells_gnuplot("particles.gnu");
	cout << "number of distanzmatrices: "<< distances.size() << endl;
	//utils::print_2dim_array( ID, m, m );
	
 	free (ID);
    
	return 0;
}
