#include "grainhdl.h"


grainhdl::grainhdl(){}
grainhdl::~grainhdl(){
	delete mymath;
	delete zeroBox;
}



void grainhdl::setSimulationParameter(){
	mymath = new mathMethods();
	// 	readInit();
	Mode = MODE; // 2 fuer lesen;  fuer erzeugen der mikrostrukture
	ngrains = PARTICLES;
	if(Mode==1) realDomainSize= M-1;			
	if(Mode==2) realDomainSize= M;
	
	dt = 1.0/double(M*M);
	h = 1.0/double(realDomainSize);
	tubeRadius = sqrt(2)*2*h + 0.00001;
	grid_blowup = BORDER; 
	
	ngridpoints = realDomainSize + (2*grid_blowup); 

	zeroBox = new LSbox();
	
	switch (Mode) {
		case 1: { 			
			ST = new double [ngrains*ngrains];
			std::fill_n(ST,ngrains*ngrains,0);
			VOROMicrostructure();
			generateRandomEnergy();
			break;
		}
		case 2: {
			ST=NULL;
			readMicrostructurefromVertex();
			break;
		}		
	}		
	construct_boundary();
	//program options:
    cout << endl << "******* PROGRAM OPTIONS: *******" << endl << endl;
    cout << "Number of Grains: " << ngrains << endl;
    cout << "simulated Timesteps: " << TIMESTEPS << endl;
	 cout << "DELTA TUBE: " << DELTA << endl;
    cout << "Timestepwidth " << dt << endl;
    cout << "Number of Gridpoints: " << ngridpoints << endl << endl;
    
    cout << endl << "******* start simulation: *******" << endl << endl;
}



void grainhdl::VOROMicrostructure(){	
	
	stringstream filename, plotfiles;
	int current_cell, cell_id;
	double x,y,z,rx,ry,rz;
	int cell_order[ngrains];
	
	grains.resize(ngrains+1);
	
	vector<LSbox*> local_grains;	

	std::vector<LSbox*>::iterator itg;
	
	// stores the centroids of the cells ; access by (3*Id, 3*ID +1, 3*ID +2)
    part_pos = new double[3*ngrains];
	
	bool randbedingung = false; // bei false ist der container halb offen?! d.h. gitterwert mit 1 werden keinem partikel zugeordnet
    if (randbedingung == false) realDomainSize-=1;	
	
	voronoicell_neighbor c;
	container con(0,1,0,1,0,1,5,5,5,randbedingung,randbedingung,randbedingung,2);
    c_loop_all vl(con);
	
	/**********************************************************/
	// Randomly add particles into the container
	
		for(int i=0;i<ngrains;i++) {
			x=utils::rnd();
			y=utils::rnd();
			z=0;
			con.put(i,x,y,z);
		}
		
	/**********************************************************/
		
    for(int i=0; i < realDomainSize; i++) for(int j= 0; j < realDomainSize; j++){
	x=double(i*h); 
	y=double(j*h); // only point within the domain
    if(con.find_voronoi_cell(x,y,z,rx,ry,rz,cell_id)){
	  cell_id= cell_id++;
	  part_pos[3*(cell_id-1)]=rx;
	  part_pos[3*(cell_id-1)+1]=ry;
	  part_pos[3*(cell_id-1)+2]=rz;

	}
    else fprintf(stderr,"# find_voronoi_cell error for %g %g 0\n",x,y);
    }  
// 	con.draw_cells_gnuplot("particles.gnu");

	int i=0;
	// iteration over all cells in the container con:
	if(vl.start()) 
	do {
	  // compute the current cell, taken out of the container
		con.compute_cell(c,vl);
		cell_order[ngrains-1-i]=(vl.pid()+1);
		
		// create a new Box for the current cell
		int box_id = vl.pid()+1;
		LSbox* newBox = new LSbox(box_id, c, part_pos,this);
		
		grains[box_id]= newBox;
		newBox->distancefunction(c, part_pos);        

	} while(vl.inc());

	delete [] part_pos;
}


void grainhdl::construct_boundary(){
	int nvertex;
	double p1 = 0.0;
	double p2 = 1.0  ;
	double vertices[16] = {p1,p1,p1,p2,p1,p2,p2,p2,p2,p2,p2,p1,p2,p1,p1,p1};
	
	boundary = new LSbox(0, 4, vertices, 0, 0, 0, this);
	grains[0]= boundary;
	boundary->distancefunction(4, vertices); 
	boundary->shape_distance();
	// get the inverse distancefunction which slope 4!
	(*boundary).plot_box(true,2,"boundary");
}

void grainhdl::readMicrostructurefromVertex(){
	FILE * levelset;	
	levelset = fopen( "lsInput.dat", "r" );
// 	levelset = fopen( "lsInput_DRAG.dat", "r" );
// 	levelset = fopen( "lsInput_quadrat.dat", "r" );

	long id;
	int nvertex;
	double phi1, PHI, phi2, xr, yr, xl, yl;
	double* vertices;
	
	fscanf(levelset, "%d\n", &ngrains);
	cout << "ngrains : " << ngrains << endl;;
	grains.resize(ngrains+1);
	
	int i=0;
	for(int nn=0; nn< ngrains; nn++){
		
		fscanf(levelset, "%ld\t %d\t %lf\t %lf\t%lf\n", &id, &nvertex, &phi1, &PHI, &phi2);
		vertices = new double [nvertex * 4];
		cout << id << " || " << nvertex << " || " << phi1 << " || " << PHI << " || " << phi2<< endl;
		
		for(unsigned int j=0; j<nvertex; j++){
			fscanf(levelset, "%lf\t %lf\t %lf\t%lf\n", &xl, &yl, &xr, &yr);	
			cout << xl << " ||\t "<< yl << " ||\t "<< xr << " ||\t "<< yr<< " ||\t " << endl;
			int k = 4*j;
			vertices[k]   = xl;
			vertices[k+1] = yl;
			vertices[k+2] = xr;
			vertices[k+3] = yr;
		}
		
		LSbox* newBox = new LSbox(id, nvertex, vertices, phi1, PHI, phi2, this);
		grains[id]= newBox;
				
	    // calculate distances	    
	    newBox->distancefunction(nvertex, vertices); 
		
		delete [] vertices;
	}
	
	ST = new double [ngrains*ngrains];			//Create ST array and fill with zeros
	std::fill_n(ST,ngrains*ngrains,0);
	
// 	compute_Boundary_Energy();
	
	for(unsigned int i=0; i<ngrains; i++){
		double buffer;
		fscanf(levelset, "%lf\t", &buffer);		
		for(unsigned int j=0; j<ngrains; j++){
			while(j < i) { 
				fscanf(levelset, "%lf\t", &buffer);
				j++;
			}
			fscanf(levelset, "%lf\t", &buffer);
			ST[j+(ngrains*i)]= (double) buffer;
			ST[i+(ngrains*j)] = ST[j+(ngrains*i)];
// 			cout << "buffer " << buffer <<endl ;
// 			fwrite(/*levelset*/, "test");
		}
		fscanf(levelset, "\n");
	} 
	fclose(levelset);

	for(unsigned int i=0; i<ngrains; i++){
		for(unsigned int j=0; j<ngrains; j++){
			cout << ST[i+(ngrains*j)] << "  \t";
		}
		cout << endl;
	}

}
 

void grainhdl::compute_Boundary_Energy(){
	double energy;
	double gamma_hagb = 0.6;
	double theta_ref = 15.0* PI / 180.;
	double theta_mis;
	for(int i= 0; i<ngrains; i++) {
		for(int j=0; j <=i; j++){	
			if(i==j) ST[j+(ngrains*i)] = 1.0;
			else{
				theta_mis = grains[i+1]->mis_ori(grains[j+1]);
				if (theta_mis <= theta_ref)	energy = gamma_hagb * ( theta_mis / theta_ref) * (1.0 - log( theta_mis / theta_ref));
				else energy = gamma_hagb;
				//richtiger Logarithmus??????
				ST[i+(ngrains*j)] = energy;
				ST[j+(ngrains*i)] = energy;
			}			
		}
	}
}
 
 
 
void grainhdl::generateRandomEnergy(){	
	const double MIN = 0.6;
	const double MAX = 1.5;
	for(int i=0; i < ngrains; i++){
		for(int j=0; j <=i; j++){			
			double zahl=(double)(rand() / (((double)RAND_MAX+1)/ (double)(MAX-MIN)))+MIN;
// 			if (i==6 && j ==5) { 
// 				ST[i+(PARTICLES*j)] = 0.1;
// 				ST[j+(PARTICLES*i)] = 0.1;
// 			}
			ST[i+(ngrains*j)] = zahl;
			ST[j+(ngrains*i)] = zahl;
			if(i==j) ST[j+(ngrains*i)] = 1.0;
		}
	} 
}
 
 
void grainhdl::convolution(){
	std::vector<LSbox*>::iterator it;
	for (it = ++grains.begin(); it !=grains.end(); it++){	
		(*it)->convolution();
	}
}
 


void grainhdl::comparison_box(){
	stringstream filename;
	vector<LSbox*>::iterator it;
	for (it = ++grains.begin(); it != grains.end(); it++){	
		(*it)->comparison();
	}
}


void grainhdl::level_set(){
	vector<LSbox*>::iterator it;
	for (it = ++grains.begin(); it != grains.end(); it++) {
		(*it)->find_contour();
	}
}


void grainhdl::redistancing(){
 
	std::vector<LSbox*>::iterator it;
	for (it = ++grains.begin(); it != grains.end(); it++) {
		(*it)->redist_box();
	}
}


void grainhdl::save_texture(){
	FILE* myfile;
	stringstream filename;
	double total_energy= 0.0;
	int numberGrains=0;
	filename << "Texture" << "_"<< loop << ".ori";	
	myfile = fopen(filename.str().c_str(), "w");
	double buffer = 0.24;
	fprintf(myfile, "%d\n", ngrains);
	vector<LSbox*> :: iterator it;
	
	for(it = ++grains.begin(); it != grains.end(); it++){
		double euler[3];
		(*mymath).quaternion2Euler( (*it)->quaternion, euler );
// 		printf( "%lf\t%lf\t%lf\t%lf\n", (*it)->quaternion[0], (*it)->quaternion[1], (*it)->quaternion[2], (*it)->quaternion[3]);
// 		printf( "%lf\t%lf\t%lf\t%lf\t%lf\n", euler[0], euler[1], euler[2], (*it)->volume, buffer);
		fprintf(myfile, "%lf\t%lf\t%lf\t%lf\t%lf\n", euler[0], euler[1], euler[2], (*it)->volume, buffer);
		total_energy += (*it)->energy;
		numberGrains+=1;
	}
    totalenergy.push_back(0.5*total_energy);
	nr_grains.push_back(numberGrains);
	cout << "Timestep " << loop << " complete:" << endl;
	cout << "Number of grains remaining in the Network :" << nr_grains.back()<< endl;
	cout << "Amount of free Energy in the Network :" << totalenergy.back()<< endl << endl;
	fclose(myfile);
}
 
  
 
 
void grainhdl::run_sim(){
	find_neighbors();
// 	determineIDs();
	for(loop=0; loop <= TIMESTEPS; loop++){		
		convolution();
		updateSecondOrderNeighbors();
		comparison_box();
		level_set();
		redistancing();
		if ( (loop % int(ANALYSESTEP)) == 0 || loop == TIMESTEPS ) {
			saveAllContourEnergies();
			save_texture();
		}
	}
	utils::CreateMakeGif();
	cout << "Simulation complete." << endl;
}  

// void grainhdl::determineIDs(){
//   	std::vector<LSbox*>::iterator it;
// 	for (it = ++grains.begin(); it !=grains.end(); it++){
// 		(*it)->determineIDs();
// 	}
// }

/*
void grainhdl::plot_contour(){
	stringstream filename;
	filename << "GrainNetwork_T" << loop << ".gnu";
	
	vector<LSbox*>::iterator it_gr; 
	for( it_gr= grains.begin(); it_gr!= grains.end(); it_gr++){
	  (*it_gr)->plot_box_contour(loop);
	}
	//change to openmp reduce
	stringstream dateinamen;
	int len=0;
	for( it_gr= grains.begin(); it_gr!= grains.end(); it_gr++){
	  if (len !=0) dateinamen << ", ";
	  dateinamen << "\"TempBox_"<< (*it_gr)->get_id() <<"_T"<< loop << ".gnu\"";
	  len+=15;
	}
	
	utils::plotContour(filename.str().c_str(), dateinamen.str().c_str(), len);
	
// 	datei << s.str();
// 	datei.close();
}
 
 */
 
  
 
void grainhdl::save_sim(){
// 	(*my_weights).plot_weightmap(ngridpoints, ID, ST, zeroBox);		
	ofstream myfile;
	myfile.open ("nr_of_rem_grains.txt");
	for(int i=0; i< nr_grains.size(); i++)
		myfile << nr_grains[i] << "\t";
	myfile.close();

	myfile.open ("energy.txt");
	for(int i=0; i<totalenergy.size(); i++)
		myfile << totalenergy[i] << "\t";
	myfile.close();
	
	if (SAVEIMAGE)utils::PNGtoGIF("test.mp4");
	//cout << "number of distanzmatrices: "<< domains.size() << endl;
}

void grainhdl::updateSecondOrderNeighbors(){
	std::vector<LSbox*>::iterator it,itc;
	for (it = ++grains.begin(); it !=grains.end(); it++){
		(*it)->add_n2o();
// 		(*it)->plot_box(false,1, "nothing");
	}
}

void grainhdl::find_neighbors(){
	std::vector<LSbox*>::iterator it,itc;
	
	for (it = ++grains.begin(); it !=grains.end(); it++)
		for (itc = ++grains.begin(); itc !=grains.end(); itc++)
			if(*it!=*itc) 
				if ((*it)->checkIntersect(*itc))
					(*it)->neighbors.push_back(*itc);

}

void grainhdl::saveAllContourEnergies(){
  	stringstream filename;
	filename<< "EnergyDistribution_T"<< loop << ".gnu";  
	ofstream dateiname;
	dateiname.open(filename.str());
	std::vector<LSbox*>::iterator it;	

	// use "#" for comments in gnuplot files for documentation purpose 
	dateiname << "#******* PROGRAM OPTIONS: *******" << endl << endl;
	dateiname << "#Number of Grains: " << ngrains << endl;
	dateiname << "#simulated Timesteps: " << TIMESTEPS << endl;
	dateiname << "#DELTA TUBE: " << DELTA << endl;
	dateiname << "#Timestepwidth " << dt << endl;
	dateiname << "#Number of Gridpoints: " << ngridpoints << endl << endl;
	
	dateiname << "set palette rgbformulae 33,13,10"<< endl;
//TODO uncomment for energy distribution	
// 	dateiname << "set cbrange[0:0.6]"<< endl;
// 	dateiname << "set cbtics 0.1" << endl;
	
	dateiname << "set cbrange[0:0.6]"<< endl;
 	dateiname << "set cbtics 1" << endl;
	
	dateiname << "set title \"Energy Distribution at Timestep " << loop <<"\""<<endl;
	dateiname << "plot \"-\" w l palette" << endl;
	for (it = ++grains.begin(); it !=grains.end(); it++)
		(*it)->plot_box_contour(loop, &dateiname, true);
	dateiname.close();
}

void grainhdl::saveAllContourLines(){	
	stringstream filename;
	filename<< "NetworkAtTime_"<< loop << ".gnu";  
	ofstream dateiname;
	dateiname.open(filename.str());
	std::vector<LSbox*>::iterator it;	
	for (it = ++grains.begin(); it !=grains.end(); it++)
		(*it)->plot_box_contour(loop, &dateiname);
	dateiname.close();
}


 
void grainhdl::clear_mem() {
	if (ST!=NULL) {delete  [] ST; }
}



