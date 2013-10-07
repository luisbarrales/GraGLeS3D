#include "grainhdl.h"


grainhdl::grainhdl(){}
grainhdl::~grainhdl(){}



void grainhdl::setSimulationParameter(){
	
// 	readInit();
	ngrains = PARTICLES;
	realDomainSize= M-1;
	
	dt = 1.0/double(M*M);
	h = 1.0/double(realDomainSize);
	
	compare_mod =1; //1 für box, 2 für domain	
	Mode = 2; // 2 für lesen;  für erzeugen der mikrostrukture
	
	 
	grid_blowup = 2*int(((double)DELTA / h)+1); //
	ngridpoints = realDomainSize + (2*grid_blowup); 
	
	my_weights = new weightmap(this);
	
	ID = new LSbox**[3];
	ID[0] = new LSbox*[ngridpoints*ngridpoints];
	ID[1] = new LSbox*[ngridpoints*ngridpoints];
	ID[2] = new LSbox*[ngridpoints*ngridpoints];
	LSbox* zeroBox = new LSbox();
	std::fill_n(ID[0],ngridpoints*ngridpoints,zeroBox);
	std::fill_n(ID[1],ngridpoints*ngridpoints,zeroBox);
	std::fill_n(ID[2],ngridpoints*ngridpoints,zeroBox);
	
	zeroBox = new LSbox();
	
	
	switch (Mode) {
		case 1: { 
			ST = new double [ngrains*ngrains];
			std::fill_n(ST,ngrains*ngrains,0);
    
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
	fscanf(Init, "%d\n", &ngrains);
	
	
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
	
	gridIDs = new int [ngridpoints* ngridpoints]; //new int[ngridpoints*ngridpoints];
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
	  cell_id= cell_id++;
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
	
	long id;
	int nvertex;
	float phi1, PHI, phi2, xr, yr, xl, yl;
	float* vertices;
	
	fscanf(levelset, "%d\n", &ngrains);
	cout << "ngrains : " << ngrains << endl;;
	
	int i=0;
	for(int nn=0; nn< ngrains; nn++){
		
		fscanf(levelset, "%ld\t %d\t %f\t %f\t%f\n", &id, &nvertex, &phi1, &PHI, &phi2);
		vertices = new float [nvertex * 4];
		cout << id << " || " << nvertex << " || " << phi1 << " || " << PHI << " || " << phi2<< endl;
		
		for(unsigned int j=0; j<nvertex; j++){
			fscanf(levelset, "%f\t %f\t %f\t%f\n", &xl, &yl, &xr, &yr);	
			cout << xl << " ||\t "<< yl << " ||\t "<< xr << " ||\t "<< yr<< " ||\t " << endl;
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
	
	ST = new double [ngrains*ngrains];			//Create ST array and fill with zeros
	std::fill_n(ST,ngrains*ngrains,0);

	for(unsigned int i=0; i<ngrains; i++){
		float buffer;
		fscanf(levelset, "%f\t", &buffer);		
		for(unsigned int j=0; j<ngrains; j++){
			while(j < i) { 
				fscanf(levelset, "%f\t", &buffer);
				j++;
			}
			fscanf(levelset, "%f\t", &buffer);
			ST[j+(ngrains*i)]= (double) buffer;
			ST[i+(ngrains*j)] = ST[j+(ngrains*i)];
// 			cout << "buffer " << buffer <<endl ;
// 			fwrite(/*levelset*/, "test");
		}
		fscanf(levelset, "\n");
	} 
	for(unsigned int i=0; i<ngrains; i++){
		for(unsigned int j=0; j<ngrains; j++){
			cout << ST[i+(ngrains*j)] << "  \t";
		}
		cout << endl;
	}
	fclose(levelset);
	
	// Ausgabe der Distanzmatrizen
	
	vector<LSbox*> grains;
	vector<LSbox*>::iterator itg;;	
	stringstream filename;
	int j;
	for (it = domains.begin(), j = 0; it !=domains.end(); it++, j++){
	  filename.str(std::string());
	  filename <<  
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
		// zugriff auf Boxen über die Domain "it"
		// Intern können verschiedenRedistancing Routinen verwendet werden

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






//=============================================================================
//
//     CONREC is a contouring subroutine for rectangularily spaced data.
//
//     It emits calls to a line drawing subroutine supplied by the user
//     which draws a contour map corresponding to real*4data on a randomly
//     spaced rectangular grid. The coordinates emitted are in the same
//     units given in the x() and y() arrays.
//
//     Any number of contour levels may be specified but they must be
//     in order of increasing value.
//
//     As this code is ported from FORTRAN-77, please be very careful of the
//     various indices like ilb,iub,jlb and jub, remeber that C/C++ indices
//     starts from zero (0)
//
//=============================================================================


int grainhdl::conrec() {
	// d               ! matrix of data to contour
	double **d = (*domain);
	
	// nc              ! number of contour levels
	// z               ! contour levels in increasing order
	
	int nc = 1;
	double *z = new double(1); z[0]=0.0;		
		
		
	std::list<domainCl>::iterator it;
	vector<LSbox*> ::iterator itl;
	
	for (it = domains.begin(); it !=domains.end(); it++){
		vector<LSbox*> grains = it->getBoxList();	
		for (itl = grains.begin(); itl != grains.end(); itl++){	
			// ilb,iub,jlb,jub ! index bounds of data matrix
			int ilb = ymin;
			int iub = ymax;
			int jlb = xmin;
			int jub = xmax;
			
			// x               ! data matrix column coordinates
			// y               ! data matrix row coordinates
			double *x = new double(ymax-ymin);
			double *y = new double(ymax-ymin);			
			for (int i=ymin; i< ymax; i++) x[i]=i*h;
			for (int j=xmin; j< xmax; i++) y[j]=j*h;
			
			
			
			int m1,m2,m3,case_value;
			double dmin,dmax,x1,x2,y1,y2;
			register int i,j,k,m;
			double h[5];
			int sh[5];
			double xh[5],yh[5];
			//===========================================================================
			// The indexing of im and jm should be noted as it has to start from zero
			// unlike the fortran counter part
			//===========================================================================
			int im[4] = {0,1,1,0},jm[4]={0,0,1,1};
			//===========================================================================
			// Note that castab is arranged differently from the FORTRAN code because
			// Fortran and C/C++ arrays are transposed of each other, in this case
			// it is more tricky as castab is in 3 dimension
			//===========================================================================
			int castab[3][3][3] =
			{
				{
				{0,0,8},{0,2,5},{7,6,9}
				},
				{
				{0,3,4},{1,3,1},{4,3,0}
				},
				{
				{9,6,7},{5,2,0},{8,0,0}
				}
			};
			for (j=(jub-1);j>=jlb;j--) {
				for (i=ilb;i<=iub-1;i++) {
					double temp1,temp2;
					temp1 = min(d[i][j],d[i][j+1]);
					temp2 = min(d[i+1][j],d[i+1][j+1]);
					dmin = min(temp1,temp2);
					temp1 = max(d[i][j],d[i][j+1]);
					temp2 = max(d[i+1][j],d[i+1][j+1]);
					dmax = max(temp1,temp2);
					if (dmax>=z[0]&&dmin<=z[nc-1]) {
						for (k=0;k<nc;k++) {
							if (z[k]>=dmin&&z[k]<=dmax) {
								for (m=4;m>=0;m--) {
									if (m>0) {
									//=============================================================
									// The indexing of im and jm should be noted as it has to
									// start from zero
									//=============================================================
									h[m] = d[i+im[m-1]][j+jm[m-1]]-z[k];
									xh[m] = x[i+im[m-1]];
									yh[m] = y[j+jm[m-1]];
									} else {
									h[0] = 0.25*(h[1]+h[2]+h[3]+h[4]);
									xh[0]=0.5*(x[i]+x[i+1]);
									yh[0]=0.5*(y[j]+y[j+1]);
									}
									if (h[m]>0.0) {
										sh[m] = 1;
									} 
									else if (h[m]<0.0) {
										sh[m] = -1;
									} 
									else sh[m] = 0;
								}
								//=================================================================
								//
								// Note: at this stage the relative heights of the corners and the
								// centre are in the h array, and the corresponding coordinates are
								// in the xh and yh arrays. The centre of the box is indexed by 0
								// and the 4 corners by 1 to 4 as shown below.
								// Each triangle is then indexed by the parameter m, and the 3
								// vertices of each triangle are indexed by parameters m1,m2,and
								// m3.
								// It is assumed that the centre of the box is always vertex 2
								// though this isimportant only when all 3 vertices lie exactly on
								// the same contour level, in which case only the side of the box
								// is drawn.
								//
								//
								//      vertex 4 +-------------------+ vertex 3
								//               | \               / |
								//               |   \    m-3    /   |
								//               |     \       /     |
								//               |       \   /       |
								//               |  m=2    X   m=2   |       the centre is vertex 0
								//               |       /   \       |
								//               |     /       \     |
								//               |   /    m=1    \   |
								//               | /               \ |
								//      vertex 1 +-------------------+ vertex 2
								//
								//
								//
								//               Scan each triangle in the box
								//
								//=================================================================
								for (m=1;m<=4;m++) {
									m1 = m;
									m2 = 0;
									if (m!=4) m3 = m+1;
									else m3 = 1;
									case_value = castab[sh[m1]+1][sh[m2]+1][sh[m3]+1];
									if (case_value!=0) {
										switch (case_value) {
											//===========================================================
											//     Case 1 - Line between vertices 1 and 2
											//===========================================================
											case 1:
											x1=xh[m1];
											y1=yh[m1];
											x2=xh[m2];
											y2=yh[m2];
											break;
											//===========================================================
											//     Case 2 - Line between vertices 2 and 3
											//===========================================================
											case 2:
											x1=xh[m2];
											y1=yh[m2];
											x2=xh[m3];
											y2=yh[m3];
											break;
											//===========================================================
											//     Case 3 - Line between vertices 3 and 1
											//===========================================================
											case 3:
											x1=xh[m3];
											y1=yh[m3];
											x2=xh[m1];
											y2=yh[m1];
											break;
											//===========================================================
											//     Case 4 - Line between vertex 1 and side 2-3
											//===========================================================
											case 4:
											x1=xh[m1];
											y1=yh[m1];
											x2=xsect(m2,m3);
											y2=ysect(m2,m3);
											break;
											//===========================================================
											//     Case 5 - Line between vertex 2 and side 3-1
											//===========================================================
											case 5:
											x1=xh[m2];
											y1=yh[m2];
											x2=xsect(m3,m1);
											y2=ysect(m3,m1);
											break;
											//===========================================================
											//     Case 6 - Line between vertex 3 and side 1-2
											//===========================================================
											case 6:
											x1=xh[m3];
											y1=yh[m3];
											x2=xsect(m1,m2);
											y2=ysect(m1,m2);
											break;
											//===========================================================
											//     Case 7 - Line between sides 1-2 and 2-3
											//===========================================================
											case 7:
											x1=xsect(m1,m2);
											y1=ysect(m1,m2);
											x2=xsect(m2,m3);
											y2=ysect(m2,m3);
											break;
											//===========================================================
											//     Case 8 - Line between sides 2-3 and 3-1
											//===========================================================
											case 8:
											x1=xsect(m2,m3);
											y1=ysect(m2,m3);
											x2=xsect(m3,m1);
											y2=ysect(m3,m1);
											break;
											//===========================================================
											//     Case 9 - Line between sides 3-1 and 1-2
											//===========================================================
											case 9:
											x1=xsect(m3,m1);
											y1=ysect(m3,m1);
											x2=xsect(m1,m2);
											y2=ysect(m1,m2);
											break;
											default:
											break;
										}
										//=============================================================
										// Put your processing code here and comment out the printf
										//=============================================================
										printf("%f %f %f %f %f\n",x1,y1,x2,y2,z[k]);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return 0;
}





// void grainhdl::compute_grain_vol(){
// 	Rf_initEmbeddedR(argc, argv);
// 	grain_vol();
// 	Rf_endEmbeddedR(0);
// }
// 
// static double grain_vol() {
//     SEXP x,y,data,level, A;
//     int errorOccurred;
// 
//     // create and evaluate 'library(splines)'
//     PROTECT(e = lang2(install("library"), mkString("sos")));
//     R_tryEval(e, R_GlobalEnv, &errorOccurred);
//     if (errorOccurred) {
//         // handle error
//     }
//     UNPROTECT(1);
// 
// //     // 'options(FALSE)' ...
// //     PROTECT(e = lang2(install("options"), ScalarLogical(0)));
// //     // ... modified to 'options(example.ask=FALSE)' (this is obscure)
// //     SET_TAG(CDR(e), install("example.ask"));
// //     R_tryEval(e, R_GlobalEnv, NULL);
// //     UNPROTECT(1);
// 
// 	data= as.matrix(read.table("C:/Users/miessen.IMM/Documents/Redistanced_matrix_2_T4_3.gnu", header=FALSE))
// 	x <- 1:nrow(data)
// 	y <- 1:ncol(data)
// 	contour(x, y, data, levels=c(0));  
// 	clines <- contourLines(x, y, data, levels=c(0))
// 	x <- clines[[1]][["x"]]
// 	y <- clines[[1]][["y"]]
// 	level <- clines[[1]][["level"]]
// 	level
// 	A = 0.5* abs( sum( x[1:(length(x)-1)]*y[2:length(x)] - y[1:(length(x)-1)]*x[2:length(x)] ) )
// 	return A;   
// }