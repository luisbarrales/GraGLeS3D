#include "grainhdl.h"


grainhdl::grainhdl(){}
grainhdl::~grainhdl(){}



void grainhdl::setSimulationParameter(){
	
// 	readInit();
	ngrains= PARTICLES;
	realDomainSize= M-1;
	
	dt = 1.0/double(M*M);
	h = 1.0/double(realDomainSize);
	
	compare_mod =1; //1 f�r box, 2 f�r domain	
	Mode = 1; // 2 f�r lesen;  f�r erzeugen der mikrostrukture
	
	 
	grid_blowup = 2*int(((double)DELTA / h)+1); //
	ngridpoints = realDomainSize + (2*grid_blowup); 
	
	my_weights = new weightmap();
	
	ID = new LSbox**[3];
	ID[0] = new LSbox*[ngridpoints*ngridpoints];
	ID[1] = new LSbox*[ngridpoints*ngridpoints];
	ID[2] = new LSbox*[ngridpoints*ngridpoints];
	LSbox* zeroBox = new LSbox();
	std::fill_n(ID[0],ngridpoints*ngridpoints,zeroBox);
	std::fill_n(ID[1],ngridpoints*ngridpoints,zeroBox);
	std::fill_n(ID[2],ngridpoints*ngridpoints,zeroBox);
	
	zeroBox = new LSbox();
	
	ST = new double [ngrains*ngrains];
	std::fill_n(ST,ngrains*ngrains,0);
	
	switch (Mode) {
		case 1: { 
			VOROMicrostructure();
			generateRandomEnergy();
			break;
		}
		case 2: readMicrostructurefromVertex(); break;
		
	}		
	
	//program options:
    cout << endl << "******* PROGRAM OPTIONS: *******" << endl << endl;
    cout << "Number of Grains: " << PARTICLES << endl;
    cout << "simulated Timesteps: " << TIMESTEPS << endl;
    cout << "Timestepwidth " << dt << endl;
    cout << "Number of Gridpoints: " << ngridpoints << endl << endl;
    
    cout << endl << "******* start simulation: *******" << endl << endl;
}





void grainhdl::readInit(){
	FILE * Init;
	Init = fopen( "Init.dat", "r" );
	fscanf(Init, "%ld\n", ngrains);
	
	
	fclose(Init);
}





void grainhdl::VOROMicrostructure(){	
	
	stringstream filename, plotfiles;
	int current_cell, cell_id;
	double x,y,z,rx,ry,rz;
	int cell_order[ngrains];
	vector<LSbox*> grains;	
	
	std::list<domainCl>::iterator it;
	std::vector<LSbox*>::iterator itg;
	
	// stores the centroids of the cells ; access by (3*Id, 3*ID +1, 3*ID +2)
    part_pos = new double[3*ngrains];
	
	bool randbedingung = false; // bei false ist der container halb offen?! d.h. gitterwert mit 1 werden keinem partikel zugeordnet
    if (randbedingung == false) realDomainSize-=1;	
	
	voronoicell_neighbor c;
	container con(0,1,0,1,0,1,5,5,5,randbedingung,randbedingung,randbedingung,2);
    c_loop_all vl(con);
	
	gridIDs = new int [ ngridpoints* ngridpoints]; //new int[ngridpoints*ngridpoints];
	std::fill_n(gridIDs, ngridpoints* ngridpoints, 0); 
	
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
			cell_id= cell_id+1;
			part_pos[3*(cell_id-1)]=rx;
            part_pos[3*(cell_id-1)+1]=ry;
            part_pos[3*(cell_id-1)+2]=rz;
			gridIDs[(i+grid_blowup)*ngridpoints + j + grid_blowup]= cell_id;
        }
        else fprintf(stderr,"# find_voronoi_cell error for %g %g 0\n",x,y);
    }  
	if(DRAW_PARTICLES)	con.draw_cells_gnuplot("particles.gnu");

	int i=0;
	// iteration over all cells in the container con:
	if(vl.start()) 
	  do {
	    // compute the current cell, taken out of the container
	    con.compute_cell(c,vl);
	    cell_order[ngrains-1-i]=(vl.pid()+1);
	    
	    // create a new Box for the current cell
		int box_id = vl.pid()+1;
	    LSbox* newBox = new LSbox(box_id, c, part_pos, grid_blowup, h, this);
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
			domains.emplace_back(ngridpoints,ngridpoints, i,INTERIMVAL, this);  
			i++;
			domains.back().addBox(newBox);            
	    } 
	    else cout << "success" << endl;
		    
	    // calculate distances	    
	    newBox->distancefunction(c, gridIDs, part_pos, grid_blowup, h);        

	  } while(vl.inc());

    int j;
    for (it = domains.begin(), j = 0; it !=domains.end(); it++, j++){
        filename.str(std::string());
        filename << "Distanzmatrix";
		grains = (*it).getBoxList();
        for (itg = grains.begin(); itg != grains.end(); itg++) {
            filename << "_" <<(*itg)->getID();
        }
        filename << ".gnu";    
		
        if (SAFEFILES){
        	cout << filename.str() << endl << endl;        
			(*it).save_domainCl(filename.str().c_str());
		}
    }
		
	delete [] part_pos;
}




void grainhdl::readMicrostructurefromVertex(){
	FILE * levelset;	
	levelset = fopen( "lsInput.dat", "r" );
	
	std::list<domainCl>::iterator it;
	
	unsigned int id;
	int nvertex;
	float phi1, PHI, phi2, xr, yr, xl, yl;
	float* vertices;
	
	fscanf(levelset, "%ld\n", &ngrains);
	cout << "ngrains : " << ngrains << endl;;
	
	
	for(int i=0; i<ngrains; i++){
		
		fscanf(levelset, "%ld\t %d\t %f\t %f\t%f\n", &id, &nvertex, &phi1, &PHI, &phi2);
		vertices = new float [nvertex * 4];
		
		for(unsigned int j=0; j<nvertex; j++){
			fscanf(levelset, "%f\t %f\t %f\t%f\n", &xl, &yl, &xr, &yr);	
			fprintf(levelset, "%f\t %f\t %f\t%f\n", xl, yl, xr, yr);	
			int k = 4*j;
			vertices[k]   = xl;
			vertices[k+1] = yl;
			vertices[k+2] = xr;
			vertices[k+3] = yr;
		}
		
		LSbox* newBox = new LSbox(id, nvertex, vertices, phi1, PHI, phi2, grid_blowup, h, this);
		
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
			domains.emplace_back(ngridpoints,ngridpoints, i,INTERIMVAL, this);  
			i++;
			domains.back().addBox(newBox);            
	    } 
	    else cout << "success" << endl;
				
	    // calculate distances	    
	    newBox->distancefunction(nvertex, vertices, grid_blowup, h); 
		delete [] vertices;
	}
	
	for(unsigned int i=0; i<ngrains; i++){
		for(unsigned int j=i; j<ngrains; j++){
			fscanf(levelset, "%f\t", &ST[j+(ngrains*i)]);
			ST[i+(ngrains*j)] = ST[j+(ngrains*i)];
		}
		fscanf(levelset, "\n");
	} 
	
	fclose(levelset);
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
	std::list<domainCl>::iterator it, itc;
	domains_copy=domains;
	for (it = domains.begin(), itc=domains_copy.begin(); it !=domains.end(); it++, itc++){	
		(*it).convolution(dt,ST,ID,(*itc), zeroBox, grid_blowup, my_weights);
	}
	if ((loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS || loop == PRINTNOW) save_conv_step();
}
 
 
void grainhdl::save_conv_step(){
	std::list<domainCl>::iterator it;
	stringstream filename;
	for (it = domains.begin(); it !=domains.end(); it++){		
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
			(*it).save_domainCl(filename.str().c_str());
			cout << filename.str() << endl << endl;
		}
	}	
}


void grainhdl::comparison_domain(){
	std::list<domainCl>::iterator it;
	stringstream filename;
	vector<LSbox*>::iterator it2;
	vector<LSbox*> grains;	
	
	domains_copy = domains;		
	
	for (it = domains.begin(); it != domains.end(); it++){		  
		(*it).comparison(domains_copy, grid_blowup);
		
		if ((loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS || loop == PRINTNOW){
			
			grains = (*it).getBoxList();			
			filename.str(std::string());
			filename << "Comparedmatrix_"<< "T"<<loop;
			for (it2 = grains.begin(); it2 != grains.end(); it2++) {
				filename << "_"<<(*it2)->getID();
			}
			
			filename << ".gnu";
			if (SAFEFILES) {
			(*it).save_domainCl(filename.str().c_str());
			}		
		}
	}
}


void grainhdl::comparison_box(){
	stringstream filename;
	std::list<domainCl>::iterator it, it_domain;
	vector<LSbox*>::iterator itLS;
//     vector<LSbox*>::iterator itLSc;
	vector<LSbox*> grains;

	domains_copy = domains;
	
	for (it = domains.begin(); it != domains.end(); it++){		     
		grains = (*it).getBoxList();
					
		for (it_domain = domains_copy.begin(); it_domain != domains_copy.end(); it_domain++){
			for (itLS = grains.begin(); itLS != grains.end(); itLS++){	
				if( (**itLS).get_status() == true ){
					if(it_domain == domains_copy.begin() )
						{ (**itLS).add_n2o(); } // copy them once for each grain in the first cycle
					(*itLS)->comparison(*it_domain, loop ); 
				}
			}
		} 			
		filename.str(std::string());
		filename << "Comparedmatrix_"<< "T"<<loop;
// 			cout << "set to dpomain" << endl;
		for (itLS = grains.begin(); itLS != grains.end(); itLS++){
			filename << "_"<<(*itLS)->getID();
			if( (**itLS).get_status() == true ){
			(**itLS).comparison_set_to_domain(ID, ngridpoints, grid_blowup);
			}
		}
		it->set_border_to_INTERIMVAL(grid_blowup); // cut the grains at der boundary of the virtual domain
// 			cout << "huhuh" << endl;	
		if ((loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS || loop == PRINTNOW){
			filename << ".gnu";
			if (SAFEFILES) {
			(*it).save_domainCl(filename.str().c_str());
			cout << filename.str() << endl << endl;
			}
		}
	}
}


void grainhdl::swap_grains(){
	std::list<domainCl>::iterator it;
	vector<LSbox*> buffer;
	for (it = domains.begin(); it != domains.end(); it++) {
		bool exist=true;
		//check domain it for intersecting grains
		exist = (*it).grainCheck(h, grid_blowup, buffer, loop); // h Gitterabstand
		if (!exist){
			cout << (*it).get_id() <<"domain leer" << endl;
			// 			domains.erase(it); it--;
			// 			cerr << "GrainCheck" ;
			// 			cin >> buffer1;
		}		
	}
	int i=domains.size()-1;
	// check if buffer is empty
	while(!buffer.empty()){
		i++;
		cout << "created a new domain" << endl;
	// 		cin >> buffer1;
		domains.emplace_back(ngridpoints,ngridpoints, i,INTERIMVAL, this);
		domains.back().grainCheck(h, grid_blowup, buffer, loop);
	}
}


void grainhdl::redistancing(){
	stringstream plotfiles;
	stringstream filename;
	vector<LSbox*>::iterator it2;
	std::list<domainCl>::iterator it;
	int i;
	nr_grains.push_back(0);
	
	for (i=0, it = domains.begin(); it != domains.end(); it++, i++) {
		//Nullstellenverfolgung:
		//cout << "Rechne Redistancing auf Boxen der Domain: " << (*it).get_id() << endl << endl;
		// zugriff auf Boxen �ber die Domain "it"
		// Intern k�nnen verschiedenRedistancing Routinen verwendet werden

		(*it).redistancing_2(h, grid_blowup);
		
		nr_grains[loop]+=(*it).get_nr_of_grains();
		
		if ( (loop % int(PRINTSTEP)) == 0 || loop == TIMESTEPS || loop == PRINTNOW){
				filename.str(std::string());
				filename << "Redistanced_matrix_";
				filename << it->get_id() << "_";
				vector<LSbox*> grains = (*it).getBoxList();
				filename << "T"<<loop;

				for (it2 = grains.begin(); it2 != grains.end(); it2++) {
						filename <<"_"<< (*it2)->getID();
				}
				
				filename << ".gnu";
				if (SAFEFILES) {
					(*it).save_domainCl(filename.str().c_str());
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
			if (imgnum < 10000) filename << "0";
			if (imgnum < 1000) filename << "0";
			if (imgnum < 100) filename << "0";
			if (imgnum < 10) filename << "0";
			filename << imgnum << ".png";
			utils::plotGnuPNG(filename.str().c_str(), plotfiles.str().c_str());
		}
	}
}
 
 
void grainhdl::run_sim(){
	find_neighbors();
	for(loop=0; loop <= TIMESTEPS; loop++){		
		convolution();
// 		domains_copy.clear();
		switch (compare_mod){
			case 1: comparison_box(); break;
			case 2: comparison_domain(); break;
		}	
		swap_grains();
// 		domains_copy.clear();
		redistancing();
	}
}  
 
 
void grainhdl::save_sim(){
	(*my_weights).plot_weightmap(ngridpoints, ID, ST, zeroBox);
	
	for(int i=0; i < realDomainSize; i++) for(int j= 0; j < realDomainSize; j++){
		gridIDs[(i+grid_blowup)*ngridpoints + j + grid_blowup] = ID[1][(i+grid_blowup)*ngridpoints + j + grid_blowup]->get_id();
	}	
	utils::save_2dim_array( gridIDs, ngridpoints, ngridpoints, "ID_Feld" );	
	
	ofstream myfile;
	myfile.open ("kinetics.txt");
	for(int i=0; i< TIMESTEPS; i++)
		myfile << nr_grains[i] << "\t";
	myfile.close();

	utils::PNGtoGIF("test.mp4");
	cout << "number of distanzmatrices: "<< domains.size() << endl;
}


void grainhdl::find_neighbors(){
	std::list<domainCl>::iterator it,itc;
	
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
							(*itl)	->neighbors.push_back(*itlc);
							(*itlc)	->neighbors.push_back(*itl);
						}
					}
					
				}
				
			}
			
		}
	}
}

 
void grainhdl::clear_mem() {
	delete  [] ST;
 	delete	[] gridIDs;  
	delete	[] ID[0];
	delete 	[] ID[1];
	delete 	[] ID[2];
	delete 	[] ID; 
	delete my_weights;
}