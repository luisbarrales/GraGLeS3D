#include "grainhdl.h"


grainhdl::grainhdl(){}
grainhdl::~grainhdl(){}



void grainhdl::setSimulationParameter(){
	
// 	readInit();
	Mode = MODE; // 2 für lesen;  für erzeugen der mikrostrukture
	ngrains = PARTICLES;
	if(Mode==1) realDomainSize= M-1;			
	if(Mode==2) realDomainSize= M;
	
	dt = 1.0/double(M*M);
	h = 1.0/double(realDomainSize);

	grid_blowup = BORDER; 
    totalenergy = new double[TIMESTEPS];
		
	ngridpoints = realDomainSize + (2*grid_blowup); 

	LSbox* zeroBox = new LSbox();
	
	switch (Mode) {
		case 1: { 			
			ST = new double [ngrains*ngrains];
			std::fill_n(ST,ngrains*ngrains,0);
			VOROMicrostructure();
			generateRandomEnergy();
			break;
		}
		case 2: {			
			readMicrostructurefromVertex();
			break;
		}		
	}		
	construct_boundary();
	//program options:
    cout << endl << "******* PROGRAM OPTIONS: *******" << endl << endl;
    cout << "Number of Grains: " << ngrains << endl;
    cout << "simulated Timesteps: " << TIMESTEPS << endl;
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
	// get the invers distancefunction
	
// 	(*boundary).save_box("boundary.gnu");
}

void grainhdl::readMicrostructurefromVertex(){
	FILE * levelset;	
// 	levelset = fopen( "lsInput.dat", "r" );
	levelset = fopen( "lsInput_DRAG.dat", "r" );
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
// 	if (((loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS )&& SAVECONV) save_conv_step();
}
 
 
// void grainhdl::save_conv_step(){
// 	std::list<domainCl>::iterator it;
// 	stringstream filename;
// 	for (it = domains.begin(); it !=domains.end(); it++){		
// 		filename.str(std::string());
// 		filename << "Convolutedmatrix_";
// 		vector<LSbox*> grains = (*it).getBoxList();
// 		vector<LSbox*>::iterator it2;
// 		filename << "T"<<loop;
// 		for (it2 = grains.begin(); it2 != grains.end(); it2++) {
// 			filename << "_"<<(*it2)->getID();
// 		}
// 		filename << ".gnu";
// 		(*it).save_domainCl(filename.str().c_str());
// 		cout << filename.str() << endl << endl;		
// 	}	
// }



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
// 	stringstream plotfiles;
// 	stringstream filename;
// 	vector<LSbox*>::iterator it2;
	std::vector<LSbox*>::iterator it;
	int i;
	nr_grains.push_back(0);
	
	for (it = ++grains.begin(); it != grains.end(); it++) {
		(*it)->redist_box();
		nr_grains[loop]+=1;
		
// 		if ( ((loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS ) && SAVEREDIST ){
// 				filename.str(std::string());
// 				filename << "Redistanced_matrix_";
// 				filename << it->get_id() << "_";
// 				vector<LSbox*> grains = (*it).getBoxList();
// 				filename << "T"<<loop;
// 
// 				for (it2 = grains.begin(); it2 != grains.end(); it2++) {
// 						filename <<"_"<< (*it2)->getID();
// 				}
// 				
// 				filename << ".gnu";
// 				
// 				(*it).save_domainCl(filename.str().c_str());
// 				cout << filename.str() << endl << endl;
// 				
// 				plotfiles << " \""<<filename.str();
// 				plotfiles << "\" matrix w l";
// 				int length = domains.size();
// 				if(i!=(length-1)) plotfiles << ",";
// 			}
	}
	
// 	if ( ((loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS )&& IMAGEOUT && SAVEREDIST){
// 		filename.str(std::string());
// 		filename << "GrainNetwork" << "_"<< loop << ".gnu";
// 		utils::plotGnu(filename.str().c_str(), plotfiles.str().c_str(), plotfiles.str().size());
// 		if(SAVEIMAGE){
// 			int imgnum = (loop/PRINTSTEP);
// 			filename.str(std::string());
// 			filename << "GrainNetwork";
// 			if (imgnum < 10000) filename << "0";
// 			if (imgnum < 1000) filename << "0";
// 			if (imgnum < 100) filename << "0";
// 			if (imgnum < 10) filename << "0";
// 			filename << imgnum << ".png";
// 			utils::plotGnuPNG(filename.str().c_str(), plotfiles.str().c_str(), plotfiles.str().size());
// 		}
// 	}
}

void grainhdl::save_texture(){
	FILE* myfile;
	stringstream filename;
	double total_energy=0;
	filename << "Texture" << "_"<< loop << ".ori";	
	myfile = fopen(filename.str().c_str(), "w");
	double buffer = 0.24;
	fprintf(myfile, "%lf\n", ngrains);
	vector<LSbox*> :: iterator it;
	for(it = ++grains.begin(); it != grains.end(); it++){
		fprintf(myfile, "%lf\t%lf\t%lf\t%lf\t%lf\n", (*it)->phi1, (*it)->PHI, (*it)->phi2, (*it)->volume, buffer);
		total_energy += (*it)->energy;
	}
    if (loop != TIMESTEPS) totalenergy[int (loop/ANALYSESTEP)]=total_energy;
    else totalenergy[TIMESTEPS-1] = 0.5 *total_energy;
	fclose(myfile);
}
 
 
 
 
 
 
 
void grainhdl::run_sim(){
	find_neighbors();
	for(loop=0; loop <= TIMESTEPS; loop++){		
		convolution();
		comparison_box();	
		level_set();
		char buffer;
		cin >> buffer;
		redistancing();
		
// 		if ( (loop % int(ANALYSESTEP)) == 0 || loop == TIMESTEPS ) {
// 			cout << "Grain Volumes after Timestep " << loop << endl;
// 			save_texture();
// 		}
	}
}  
 
 
 
  
 
void grainhdl::save_sim(){
// 	(*my_weights).plot_weightmap(ngridpoints, ID, ST, zeroBox);		
	ofstream myfile;
	myfile.open ("nr_of_rem_grains.txt");
	for(int i=0; i< TIMESTEPS; i++)
		myfile << nr_grains[i] << "\t";
	myfile.close();

	myfile.open ("energy.txt");
	for(int i=0; i<TIMESTEPS; i++)
		myfile << totalenergy[i] << "\t";
	myfile.close();
	
	if (SAVEIMAGE)utils::PNGtoGIF("test.mp4");
	//cout << "number of distanzmatrices: "<< domains.size() << endl;
}


void grainhdl::find_neighbors(){
	std::vector<LSbox*>::iterator it,itc;
	
	for (it = ++grains.begin(); it !=grains.end(); it++)
		for (itc = ++grains.begin(); itc !=grains.end(); itc++)
			if(*it!=*itc) 
				if ((*it)->checkIntersect(*itc))
					(*it)->neighbors.push_back(*itc);

}

 
void grainhdl::clear_mem() {
	delete  [] ST;
}



