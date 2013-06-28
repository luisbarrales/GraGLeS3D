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

/*		Stand: JUNI 2013 --------------------------------------------------------------------------------------------------*/

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
	
	
/*************************************************************/
// Init

    char buffer2;
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
    std::list<matrix>::iterator it, itc, it_domain;
    vector<LSbox*> grains;
    vector<LSbox*>::iterator it2,it2c;
    
    voronoicell_neighbor c;
    
    bool randbedingung = false; // bei false ist der container halb offen?! d.h. gitterwert mit 1 werden keinem partikel zugeordnet
    if (randbedingung == false) m = int(M)-1;
    
    const double h = 1.0/double(m); // there are m-2*gridbloup grid points in each direction in the domain
    const int grid_blowup = int(((double)DELTA / h)+1); // number of grid points to extract the domain at each boundary
    
    int resized_m = m + (2*grid_blowup); 
    LSbox ***ID; // array to asign a cell id to each grid point
    
    ID = new LSbox**[2];
	ID[0] = new LSbox*[resized_m*resized_m];
	ID[1] = new LSbox*[resized_m*resized_m];
        
    double  (*fp)(double,int, int, int); // function pointer
    fp = &kernel;
    
    container con(0,1,0,1,0,1,5,5,5,randbedingung,randbedingung,randbedingung,2);
    c_loop_all vl(con);
	
	vector<int> nr_grains(TIMESTEPS+1);

/**********************************************************/

    
    //program options:
    cout << endl << "******* PROGRAM OPTIONS: *******" << endl << endl;
    cout << "Number of Grains: " << PARTICLES << endl;
    cout << "simulated Timesteps: " << TIMESTEPS << endl;
    cout << "Timestepwidth " << dt << endl;
    cout << "Number of Gridpoints: " << M << endl << endl;
    
    cout << endl << "******* start simulation: *******" << endl << endl;
	
    
/**********************************************************/
// Randomly add particles into the container
   
    for(int i=0;i<particles;i++) {
        x=utils::rnd();
        y=utils::rnd();
        z=0;
        con.put(i,x,y,z);
    }
    
/**********************************************************/


/**********************************************************/
// generate random or deterministic surface tension coefficients

	double *ST;
	ST = (double*) calloc (PARTICLES*PARTICLES,sizeof(double));
	const double MIN = 0.6;
	const double MAX = 1.0;

	for(int i=0; i < PARTICLES; i++){
		for(int j=0; j < PARTICLES; j++){			
			double zahl=(double)(rand() / (((double)RAND_MAX+1)/ (double)(MAX-MIN)))+MIN;
			ST[i+(PARTICLES*j)] = zahl;
			ST[j+(PARTICLES*i)] = zahl;
		}
	}
	
/**********************************************************/


/**********************************************************/
// 	Output the Voronoi cells to a file, in the gnuplot format
	
	con.draw_cells_gnuplot("random_points_v.gnu");
	
/**********************************************************/
    

/**********************************************************/	
// find cell information for each grid point

    int *gridIDs;
	gridIDs = new int [resized_m*resized_m]; //new int[resized_m*resized_m];
    for(int i=0; i < m; i++) for(int j= 0; j < m; j++){
        x=double(i*h); 
	y=double(j*h); // only point within the domain
        if(con.find_voronoi_cell(x,y,z,rx,ry,rz,cell_id)){
			cell_id= cell_id+1;
			part_pos[3*(cell_id-1)]=rx;
            part_pos[3*(cell_id-1)+1]=ry;
            part_pos[3*(cell_id-1)+2]=rz;
			gridIDs[(i+grid_blowup)*resized_m + j + grid_blowup]= cell_id;
        }
        else fprintf(stderr,"# find_voronoi_cell error for %g %g 0\n",x,y);
    }  
	if(DRAW_PARTICLES)	con.draw_cells_gnuplot("particles.gnu");

/**********************************************************/



/*********************************************************************************/
// Initialisation: Voronoizellen -> Box
    
	int i=0;
	// iteration over all cells in the container con:
	if(vl.start()) 
	  do {
	    // compute the current cell, taken out of the container
	    con.compute_cell(c,vl);
	    cell_order[particles-1-i]=(vl.pid()+1);
	    
	    // create a new Box for the current cell
		int box_id = vl.pid()+1;
	    LSbox* newBox = new LSbox(box_id, c, part_pos, grid_blowup, h);
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
		cout << "failed; creating new domain:" << endl;
			    domains.emplace_back(resized_m,resized_m, i,INTERIMVAL);  // = push_back
			    i++;
		domains.back().addBox(newBox);            
	    } else cout << "success" << endl;
		    
	    // calculate distances	    
	    newBox->distancefunction(c, gridIDs, part_pos, grid_blowup, h);        

	  } while(vl.inc());

    int j;
    for (it = domains.begin(), j = 0; it !=domains.end(); it++, j++){
        filename.str(std::string());
        filename << "Distanzmatrix";
		grains = (*it).getBoxList();
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


/*********************************************************************************/
// Determine the box' neighbors 

for (it = domains.begin(); it !=domains.end(); it++){
	vector<LSbox*> grains = it->getBoxList();
	vector<LSbox*> grainstoComp;
	vector<LSbox*> ::iterator itl, itlc;
	cout <<"determine Neighbors for domain: " << (*it).get_id() << endl;
	itc = it++;
	it--;
	for (itc; itc!=domains.end(); itc++){
		grainstoComp = itc->getBoxList();
		
		for (itl = grains.begin(); itl != grains.end(); itl++){
			
			for(itlc = grainstoComp.begin(); itlc != grainstoComp.end(); itlc++){
				
				if ((*itlc)->getID()!=(*itl)->getID()){ // wozu ist diese abfrage? --- Damit die Boxen nicht Nachbarn von sich selbst sind
					if((*itl)->checkIntersect(*itlc)){
						(*itl)->neighbors.push_back(*itlc);
						(*itlc)->neighbors.push_back(*itl);
					}
				}
				
			}
			
		}
		
	}
	
// 	for(it2c=grains.begin();it2c!=grains.end();it2c++){ (*it2c)->plot_box(false); }
}

/*********************************************************************************/

	
	
	
	
/*********************************************************************************/	
/*********************************************************************************/
/*********************************************************************************/
// MAIN LOOP

for(int loop=0; loop <= TIMESTEPS; loop++){

	stringstream plotfiles;
	plotfiles.str(std::string());  

/*********************************************************************************/
// Convolution simulates grain growth
	
	// ACHTUNG hier kopieren wir die ganze LISTE!
	domains_copy=domains;

	for (it = domains.begin(); it !=domains.end(); it++){	
		(*it).convolution(dt, ID);
	}

	// Output	
	for (it = domains.begin(); it !=domains.end(); it++){		
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
// 				cout << filename.str() << endl << endl;
			}
		}
	}

/*********************************************************************************/	


/*************************************************************************************/
// Comparison Step: step 0.5 *(A_k(x) - max A_i(x) | i!=k)
	
	if(DOMAINCOMPARISON==true){
		/*********************************************************************************/
		// Comparison Step with domains
		/*********************************************************************************/
		domains_copy = domains;		
		for (it = domains.begin(); it != domains.end(); it++){		  
			(*it).comparison(domains_copy, grid_blowup);
			
			if ((loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS || loop == PRINTNOW){
				
				vector<LSbox*>::iterator it2;
				vector<LSbox*> grains;				
				grains = (*it).getBoxList();
				
				filename.str(std::string());
				filename << "Comparedmatrix_"<< "T"<<loop;
				for (it2 = grains.begin(); it2 != grains.end(); it2++) {
				filename << "_"<<(*it2)->getID();
				}
				
				filename << ".gnu";
				if (SAFEFILES) {
				(*it).save_matrix(filename.str().c_str());
// 				cout << filename.str() << endl << endl;
				}		
			}
		}
	}
	
	else {
		/*********************************************************************************/
		// Comparison Step on boxes
		/*********************************************************************************/	
		
		vector<LSbox*> grains;	
		vector<LSbox*>::iterator it2;
		domains_copy = domains;
		char buffer1;
		for (it = domains.begin(); it != domains.end(); it++){		     
			grains = (*it).getBoxList();
						
			for (it_domain = domains_copy.begin(); it_domain != domains_copy.end(); it_domain++){
				for (it2 = grains.begin(); it2 != grains.end(); it2++){		
					if(it_domain == domains_copy.begin()){ (**it2).add_n2o(); } // copy them once for each grain in the first cycle
					(*it2)->comparison(*it_domain, loop);
				}
			} 			
			filename.str(std::string());
			filename << "Comparedmatrix_"<< "T"<<loop;
						
			for (it2 = grains.begin(); it2 != grains.end(); it2++){
				filename << "_"<<(*it2)->getID();
				(**it2).comparison_set_to_domain();
			}
			it->set_border_to_INTERIMVAL(grid_blowup); // cut the grains at der boundary of the virtual domain
			
			if ((loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS || loop == PRINTNOW){
				filename << ".gnu";
				if (SAFEFILES) {
				(*it).save_matrix(filename.str().c_str());
				cout << filename.str() << endl << endl;
				}
			}
		}
	}
	
/*************************************************************************************/


//************************************************************************************/
// Redistancing Step

	/*****************************************************************/
	// Nullstellenverfolgung:
	// Speichert Nullstellen als NNZ-OBjekt in jeder Box
	// Testet Boxen auf Konlikte
	// Verschiebt "Konfliktboxen" in andere Domain
		
	// Ist dieser Schritt in jedem Zeitschritt notwenig????
	/*****************************************************************/
	
	
	/*****************************************************************/
	// checking for existence ++ resizing the boxes ++ swaping grains
	
	vector<LSbox*> buffer;
	for (it = domains.begin(); it != domains.end(); it++) {
		bool exist=true;
		//check domain it for intersecting grains
		exist = (*it).grainCheck(h, grid_blowup, buffer); // h Gitterabstand
		if (!exist){
			cout << (*it).get_id() <<"domain leer" << endl;
// 			domains.erase(it); it--;
// 			cerr << "GrainCheck" ;
// 			cin >> buffer1;
		}
		
	}
	// check if buffer is empty
	while(!buffer.empty()){
		cout << "created a new domain" << endl;
		domains.emplace_back(resized_m,resized_m, i,INTERIMVAL);
		domains.back().grainCheck(h, grid_blowup, buffer);
	}
	
	/*****************************************************************/
	
	
	/*****************************************************************/
	// fast sweeping
	
	for (i=0, it = domains.begin(); it != domains.end(); it++, i++) {
		//Nullstellenverfolgung:
		//cout << "Rechne Redistancing auf Boxen der Domain: " << (*it).get_id() << endl << endl;
		// zugriff auf Boxen über die Domain "it"
		// Intern können verschiedenRedistancing Routinen verwendet werden

		(*it).redistancing_2(h, grid_blowup);		

		nr_grains[loop]+=(*it).get_nr_of_grains();
		
		if ( (loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS || loop == PRINTNOW){
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
				int length = domains.size();
				if(i!=(length-1)) plotfiles << ",";
			}
		}
		if ( (loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS || loop == PRINTNOW ){
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
		
	/*****************************************************************/	
		
	/************************************************************************************/
			
	cout << "Timestep: "<< loop << " complete. " ;
	cout << "Number of remaining grains: "<< nr_grains[loop] << endl << endl;
}



/*******************************************************************************************/
// end of simulation
// Endausgabe
/*******************************************************************************************/
	

	ofstream myfile;
	myfile.open ("kinetics.txt");
	for(int i=0; i< TIMESTEPS; i++)
		myfile << nr_grains[i] << "\t";
	myfile.close();

	utils::PNGtoGIF("test.mp4");
	cout << "number of distanzmatrices: "<< domains.size() << endl;
//	//utils::print_2dim_array( ID, m, m );
	
 	delete	[] gridIDs;  
	delete	[] ID[0];
	delete 	[] ID[1];
	delete 	[] ID;
	return 0;
}
