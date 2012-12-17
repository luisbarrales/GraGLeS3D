/*----------------------------------------------------------------------------------------------------------------------*/	
	
/*		Programm zur Simulation der Kornvergrößerung----------------------------------------------------------------------*/

/*		Ziel: Parallelisierung von bekannten Methoden---------------------------------------------------------------------*/

/*		Aufbau: Der Algorithmus wird in vier Teile gegliedert, sodass einzelne 
		Module ausgetauscht werden können:
			1. Initialisierungsschritt (siehe "initialisation.h")
			2. Wachstumsschritt (surfacemotion.h)
			3. Korrekturschritt (correction.h)
			4. Redistanzierungsschritt (redistancing.h)
	
		Die Schritte 2 und 3 stellen den Kern des Algorithmus dar, sie werden je Zeitschritt einmal ausgeführt. 
		Das "Redistancing" ist aus Stabilitätsgründen notwendig - jedoch nicht in jedem Zeitschritt.

		Beschreibung der Kornoberflächen: Hierzu wird ein spezieller Level Set Ansatz gewählt. Wir identifizieren 
		die Oberfläche eines Korns über die Nullstellenmenge einer Hilfsfunktion \Phi[i,j] = dist(x_i,j , \Gamma). 
		Diese Funktion wird auf einem äquidistanten Gitter ausgewertet und trägt die Werte des kürzesten 
		vorzeichenbehafteten Abstandes zur Kornoberfläche. Die Steigung ist in Normalenrichtung zur Oberfläche in 
		jedem Punkt 1. Wir werten dies Funktion nur auf eine eps-Schlauch um die zu beschreibende Oberfläche aus und 
		können so mehrere/viele Körner in einem Gitter speichern.----------------------------------------------------------*/
	
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
	const int particles = int(PARTICLES);
	const double dt = 1./(double(M*M));

	double x,y,z,rx,ry,rz;
	int m = M, current_cell, cell_id;	
	int cell_order[particles]; // stores the order of computed cells (vl.pid())
	
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

	bool randbedingung = false; // f�r false ist der container halb offen?! d.h. gitterwert mit 1 werden keinem partikel zugeordnet
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

	
	
	/*********************************************************************************/
	// Randomly add particles into the container
	/*********************************************************************************/
	for(int i=0;i<particles;i++) {
		x=utils::rnd();
		y=utils::rnd();
	// 		cout << x <<y <<endl;
		z=0;
		con.put(i,x,y,z);
	}



/*********************************************************************************/
// 	Output the Voronoi cells to a file, in the gnuplot format
	con.draw_cells_gnuplot("random_points_v.gnu");

	

/*********************************************************************************/
// find cell information fpr each grid point
/*********************************************************************************/
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
		
		compared_dist.push_front(matrix(resized_m,resized_m, vl.pid()));
		// write the distancematrices to outputfile:
		filename.str(std::string());
		filename << "Distanzmatrix" << vl.pid() << ".gnu";
		cout << filename.str() << endl << endl;
		if (SAFEFILES)(*distances.begin()).save_matrix(filename.str().c_str());
	}
	while(vl.inc());	


	
	
/*******************************************************************************************/
//// Starting Simulation of "Timesteps" times the Convolution/Comparison/Redistancing Steps
/*******************************************************************************************/
/*******************************************************************************************/

	for(int loop=0; loop <= TIMESTEPS; loop++){
		stringstream plotfiles;
		plotfiles.str(std::string());
		/*********************************************************************************/	
		// Convolution simulates grain growth
		/*********************************************************************************/
		
		for (it = distances.begin(); it !=distances.end(); it++){	
			bool exist = false;
		    if (MODUS)	exist = (*it).discrete_convolution(dt, h, grid_blowup, fp);
			else	(*it).convolution(dt);
			// 	(*it).five_point_formula(dt, h);
			if ((loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS){
				filename.str(std::string());
				filename << "Convoluted_matrix" << (*it).get_id() << "_"<< loop << ".gnu";
				cout << filename.str() << endl << endl;	
				if (SAFEFILES)(*it).save_matrix(filename.str().c_str());
			}
		}
		
		compared_dist=distances;

		/*********************************************************************************/
		// Comparison Step: step 0.5 *(A_k(x) - max A_i(x) | i!=k)
		/*********************************************************************************/
		// Create a list for storing the new distances after comparison
	
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
					cout << filename.str() << endl << endl;
					if (SAFEFILES) (*itc).save_matrix(filename.str().c_str());
				}
			}
		}

		

		/*********************************************************************************/
		// Redistancing Step
		/*********************************************************************************/
		
		int length = compared_dist.size();
			for (i=0, itc= compared_dist.begin(); itc != compared_dist.end(); itc++, i++){
			(*itc).redistancing(h, grid_blowup);
			if ( (loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS){
				filename.str(std::string());
				filename << "Redistanced_matrix" << (*itc).get_id() << "_"<< loop << ".gnu";
				cout << filename.str() << endl << endl;
				if (SAFEFILES) (*itc).save_matrix(filename.str().c_str());
				plotfiles << " \""<<filename.str();
				plotfiles << "\" matrix w l";
				if(i!=(length-1)) plotfiles << ",";
			}	
		}
		if ( (loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS){
			filename.str(std::string());
			filename << "GrainNetwork" << "_"<< loop << ".gnu";
			utils::plotGnu(filename.str().c_str(), plotfiles.str().c_str());
		}
		
		distances = compared_dist;
	}    

/*******************************************************************************************/
// end of simulation
/*******************************************************************************************/



/*********************************************************************************/
/*********************************************************************************/
	
	con.draw_cells_gnuplot("particles.gnu");	
	cout << "number of distanzmatrices: "<< distances.size() << endl;
	//utils::print_2dim_array( ID, m, m );
	
 	free (ID);

	return 0;
}
