#include "box.h"



LSbox::LSbox() {}

LSbox::LSbox(int id, int xmin, int xmax, int ymin, int ymax, double phi1, double PHI, double phi2): id(id), xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), phi1(phi1), PHI(PHI), phi2(phi2) {}

LSbox::LSbox(int aID, voro::voronoicell_neighbor& c, double *part_pos, int grid_blowup, double h, grainhdl* owner) : id(aID), phi1(0), PHI(0), phi2(0), nvertices(0), handler(owner) {
    
    // determine size of grain
	xmax = 0; xmin = handler->get_ngridpoints(); 
	ymax = 0; ymin = xmin;
	
    vektor x1(2), x2(2);
    vector<double> vv;
    exist = true;
	distance = NULL;
	c.vertices(part_pos[3*(id-1)],part_pos[3*(id-1)+1],part_pos[3*(id-1)+2],vv);
    for(int ii=0;ii<c.p;ii++) {
        for(int jj=0;jj<c.nu[ii];jj++) {
            
            int k=c.ed[ii][jj];
            x1[0]=vv[3*ii];x1[1]=vv[3*ii+1];
            x2[0]=vv[3*k]; x2[1]=vv[3*k+1];
            
			//	for convention: 
			//	x[i][j]:
			//	i = Zeilenindex(y-direction)   
			// 	j = Spaltenindex(x-direction)
			
			// check for "Zeilen" Minima/Maxima
			if (x1[0]/h < ymin) ymin = x1[0]/h;
            if (x2[0]/h < ymin) ymin = x2[0]/h;
            
            if (x1[0]/h > ymax) ymax = (x1[0]/h);
            if (x2[0]/h > ymax) ymax = (x2[0]/h);
			
			// check for "Spalten" Minima/Maxima
            if (x1[1]/h < xmin) xmin = x1[1]/h;
            if (x2[1]/h < xmin) xmin = x2[1]/h;
            
            if (x1[1]/h > xmax) xmax = (x1[1]/h);
            if (x2[1]/h > xmax) xmax = (x2[1]/h);
          
            
        }
    }
    
	xmax += 2*grid_blowup;
	ymax += 2*grid_blowup;
        
   cout << "made a new box: xmin="<<xmin<< " xmax="<<xmax <<" ymin="<<ymin << " ymax="<<ymax<<endl;
}


LSbox::LSbox(int id, int nvertex, double* vertices, double phi1, double PHI, double phi2, int grid_blowup, double h, grainhdl* owner) : id(id), phi1(phi1), PHI(PHI), phi2(phi2), nvertices(nvertex), handler(owner){
    
    // determine size of grain
    xmax = 0; xmin = handler->get_ngridpoints(); 
	ymax = 0; ymin = xmin;
	    
    vektor x1(2), x2(2);
    exist = true;
	distance = NULL;
    for (unsigned int k=0; k < nvertex; k++){       
		x1[0]=vertices[(4*k)+1]; x1[1]=vertices[4*k];
		x2[0]=vertices[(4*k)+3]; x2[1]=vertices[(4*k)+2];
		
		//	for convention: 
		//	x[i][j]:
		//	i = Zeilenindex(y-direction)   
		// 	j = Spaltenindex(x-direction)
		
		// check for "Zeilen" Minima/Maxima
		if (x1[0]/h < ymin) ymin = x1[0]/h;
		if (x2[0]/h < ymin) ymin = x2[0]/h;
		
		if (x1[0]/h > ymax) ymax = (x1[0]/h);
		if (x2[0]/h > ymax) ymax = (x2[0]/h);
		
		// check for "Spalten" Minima/Maxima
		if (x1[1]/h < xmin) xmin = x1[1]/h;
		if (x2[1]/h < xmin) xmin = x2[1]/h;
		
		if (x1[1]/h > xmax) xmax = (x1[1]/h);
		if (x2[1]/h > xmax) xmax = (x2[1]/h);
    }
	xmax += 2*grid_blowup;
	ymax += 2*grid_blowup;
        
   cout << "made a new box: xmin="<<xmin<< " xmax="<<xmax <<" ymin="<<ymin << " ymax="<<ymax<<endl;

}





LSbox::~LSbox() {
}

void LSbox::setDomain(domainCl* aDomain) {
    domain = aDomain;
}

int LSbox::getID() {
    return id;
}
LSbox LSbox::distancefunction(int nvertex, double* vertices, int grid_blowup, double h){
// 	plot_box(false);
	int i,j,k;
	double d, dmin,lambda;
	int m=handler->get_ngridpoints();
	vektor u(2), a(2), p(2), x1(2), x2(2);

	for (i=ymin;i<ymax;i++){ // ¸ber gitter iterieren
	  for (j=xmin;j<xmax;j++){
            dmin=1000.;
            p[0]=(i-grid_blowup)*h; p[1]=(j-grid_blowup)*h;            
            
            for(int k=0; k < nvertices; k++) {                
				x1[0]=vertices[(4*k)+1]; x1[1]=vertices[4*k];
				x2[0]=vertices[(4*k)+3]; x2[1]=vertices[(4*k)+2];				
				if (x1 != x2){
					a = x1;
					u = x2-x1;
					lambda=((p-a)*u)/(u*u); 
					
					if(lambda <= 0.) 				d = (p-x1).laenge();
					if((0. < lambda) && (lambda < 1.)) 		d = (p-(a+(u*lambda))).laenge();
					if(lambda >= 1.) 				d = (p-x2).laenge();
// 					if(((grid_blowup < i) && (i < (m- grid_blowup))) && ((grid_blowup < j) && (j < (m- grid_blowup)))) {
// 						d=abs(d);
// 					}
// 					else d= abs(d);
					d= abs(d);
					if(abs(d)< abs(dmin)) dmin=d;
				}
            }
			// 			(*domain)[i][j]= dmin;
			if (abs(dmin) < DELTA) (*domain)[i][j]= dmin;
            else (*domain)[i][j]= DELTA * utils::sgn(dmin);
        }
	}
	int count = 0;
	for (i=xmin;i<xmax;i++){ // ¸ber gitter iterieren
		j=ymin;
		count = 0;
		while( j<ymax  && count < 1) {
			(*domain)[j][i] = - abs((*domain)[j][i]);
			if ( -((*domain)[j][i]) <=  h ) count++;
			j++;
		} 		
		j=ymax-1;
		count =0;
		while( j>=ymin && count < 1) {
			(*domain)[j][i] = - abs((*domain)[j][i]);	
			if ( -((*domain)[j][i]) <= h ) count++;
			j--;
		} 
	}

	for (j=ymin;j<ymax;j++){ // ¸ber gitter iterieren
		i=xmin;
		count = 0;
		while( i<xmax  && count < 1 ) {
			(*domain)[j][i] = - abs((*domain)[j][i]);
			if ( -((*domain)[j][i]) <=  h ) count++;
			i++;
		} 
		i=xmax-1;
		count =0;
		while( i>=xmin   && count < 1  ) {			
			(*domain)[j][i] = - abs((*domain)[j][i]);	
			if ( -((*domain)[j][i]) <=  h ) count++;
			i--;
		} 
	}
	return(*this);
}


LSbox LSbox::distancefunction(voro::voronoicell_neighbor& c, int *gridIDs, double *part_pos, int grid_blowup, double h ){
	int i,j,k;
	double d, dmin,lambda;
	int m=domain->get_m();
	int n=domain->get_n();
	vektor u(2), a(2), p(2), x1(2), x2(2);
	vector<double> vv;
	c.vertices (part_pos[3*(id-1)],part_pos[3*(id-1)+1],part_pos[3*(id-1)+2],vv);
	double domain_vertices[] = {0.,0.,1.,0.,1.,1.,0.,1.,0.,0.}; // array of vertices to loop over
	distance=NULL;
	for (i=ymin;i<ymax;i++){ // ¸ber gitter iterieren
	  for (j=xmin;j<xmax;j++){
            dmin=1000.;
            p[0]=(i-grid_blowup)*h; p[1]=(j-grid_blowup)*h;            
            
            for(int ii = 0;ii < c.p; ii++) {
                for(int jj = 0;jj < c.nu[ii]; jj++) {
                    
                    k=c.ed[ii][jj];
                    
                    x1[0] = vv[3*ii]; x1[1] = vv[3*ii+1];
                    x2[0] = vv[3*k];  x2[1] = vv[3*k+1];
                    
                    if (x1 != x2){
                        a = x1;
                        u = x2-x1;
                        lambda=((p-a)*u)/(u*u); 
                        
                        if(lambda <= 0.) 				d = (p-x1).laenge();
                        if((0. < lambda) && (lambda < 1.)) 		d = (p-(a+(u*lambda))).laenge();
                        if(lambda >= 1.) 				d = (p-x2).laenge();
						if((id == gridIDs[i*m +j]) && ((grid_blowup < i) && (i < (m- grid_blowup))) && ((grid_blowup < j) && (j < (m- grid_blowup)))) d=abs(d);
						else d= -abs(d);
			// 			cout << "ID klappt" << endl;
			// 			char buffer;
			//                         cin >> buffer ;
						if(abs(d)< abs(dmin)) dmin=d;
                    }
                }
            }
// 			(*domain)[i][j]= dmin;
            if (dmin < DELTA) (*domain)[i][j]= dmin;
            else (*domain)[i][j]= DELTA * utils::sgn(dmin);
        }
	}
	return(*this);
}

void LSbox::setZeros(double h, int grid_blowup, int loop) {
    zeros.clear();
    int m = handler->get_ngridpoints();
    stringstream s;
    int first_i, first_j;
    int current_i, current_j;
    int next_i, next_j;
    char direction = 1; // x+
    // directions 0 = y-  //  2 = y+  //  3 x-  (y- = up // y + down; (0,0)left upper corner)
	double pointx; 
	double pointy;
	exist = false;
	
	int dist = ymax - ymin;
	int i = ymin+ int(dist/2);
    // look for zero in row y
    for (int j = xmin; j < xmax-1; j++) {
        if ((*domain)[i][j] * (*domain)[i][j+1] <= 0) {
            first_i = i; 	first_j = j; 
            current_i = i; 	current_j = j;
	    next_i =i; 		next_j = j+1;
	    exist= true;
	    break;
        }
	}
	if (!exist) {
		cout << "search in y-direction" << endl;
		int dist = xmax - xmin;
		int j = xmin+ int(dist/2);
		for (int i = ymin; i < ymax-1; i++) {	
			if ((*domain)[i][j] * (*domain)[i+1][j] <= 0) {
				first_i = i; 	first_j = j; 
				current_i = i; 	current_j = j;
				next_i =i+1; 	next_j = j;
				exist= true;
				direction = 2;
				cout << "boundary found"<< endl;
				break;
			}
		}
	}
	if (!exist) {
		cout << "no boundary found in box "<<  id << endl;
		exist = false;
		cout << "set exist status to false: " << exist << endl;
		return;
	}

    // begin zero-tracking and interpolations
    bool newZero = true;
	int sgn = -1; //(1 = left turn; -1 right turn)  

	// reste the min and:
	xmax = 0; xmin = m; ymax = 0; ymin = m;
	SPoint point;
	vector<SPoint> points;
	
// 	 begin search
	
    while (newZero) {		
	  current_i = next_i;
	  current_j = next_j;
	  
// check for size change
	  if (current_j-grid_blowup < xmin) 		
		xmin = current_j - grid_blowup;
	  else if (current_j > xmax-grid_blowup) 	
		xmax = current_j + grid_blowup;
	  if (current_i < ymin+grid_blowup) 		
		ymin = current_i - grid_blowup;
	  else if (current_i > ymax-grid_blowup) 	
		ymax = current_i + grid_blowup;
	  
	  // change search directions
	  //(1 = left turn; -1 right turn)  
	  sgn = utils::sgn((*domain)[current_i][current_j]);            
	  if (sgn == 0) { 
		if (direction == 0) 	  {next_i = current_i-1;next_j = current_j;}
		else if (direction == 2)  {next_i = current_i+1;next_j = current_j;}
		else if (direction == 1)  {next_j = current_j+1;next_i = current_i;}
		else if (direction == 3)  {next_j = current_j-1;next_i = current_i;}
		if ((*domain)[next_i][next_j]>0) (*domain)[current_i][current_j]= 0.000001;
		else	(*domain)[current_i][current_j]=-0.000001;
		  sgn = utils::sgn((*domain)[current_i][current_j]); 
	  } 
	  // search next zero
	  // turn  right
	  direction = (direction+sgn+4) %4; 
	  
	  bool foundnext = false;
	  for (int i = 0; i < 3; i++) {
		  if (direction == 0) 	 	{next_i = current_i-1;next_j = current_j;}
		  else if (direction == 2)  {next_i = current_i+1;next_j = current_j;}
		  else if (direction == 1)  {next_j = current_j+1;next_i = current_i;}
		  else if (direction == 3)  {next_j = current_j-1;next_i = current_i;}
		  
//            cerr << current_i  <<"  "<< current_j <<"  "<< next_i <<"  "<< next_j <<endl ;
		  if ((*domain)[current_i][current_j] * (*domain)[next_i][next_j] <= 0) {
			  foundnext = true;
			  if(((*domain)[current_i][current_j] * (*domain)[next_i][next_j] != 0.0))
			  {
				double slope = (*domain)[next_i][next_j]-(*domain)[current_i][current_j];
				slope = -1.0 *slope;
				if (direction == 1) {	 
				  point.x= current_j + (*domain)[current_i][current_j]/slope;
				  point.y = current_i;
				}
				else if (direction == 3){	 
				  point.x = current_j - (*domain)[current_i][current_j]/slope;
				  point.y = current_i;
				}
				else if (direction == 0) {	 				
				  point.y =  current_i - (*domain)[current_i][current_j]/slope;
				  point.x = current_j;
				}
				else if (direction == 2){	
				  point.y =  current_i + (*domain)[current_i][current_j]/slope;
				  point.x = current_j;
				}		 
			  }
			  else {
				cerr << "levelset on gridpoint  " << current_i << "\t" << current_j<<"\t" << (*domain)[current_i][current_j]<< endl;
				
				if ((*domain)[current_i][current_j] == 0.0) { point.y =  current_i;  point.x = current_j; }
				else if ((*domain)[next_i][next_j]  == 0.0) { point.y =  next_i;  point.x = next_j; }
			  }
			  points.emplace_back(point);			 
			  break; //springen aus der for-schleife, falls nullstelle gefunden wird
		  }        		  
		  
		  current_j = next_j; 
		  current_i = next_i;
		  // turn further if no zero found
		  direction = (direction+sgn+4)%4; 
	  }
	  if (!foundnext) {
		  break; //springen aus der while-schleife, falls keine nullstelle mehr gefunden wird
	  }
	  
	  // check if completed round
	  if (current_j == first_j && current_i == first_i) {
		  newZero = false;	//springen aus der while-schleife, wir haben eine geschlosse kurve gefunden
	  }
    }
    
    // compute Volume and Energy
    if ( (loop % int(ANALYSESTEP)) == 0 || loop == TIMESTEPS ) {
		vector<SPoint>::iterator volumeit=points.begin();
		energy = 0;
		volume = 0;
		double px, py;
		double gamma_hagb = 0.6;
		double theta_ref = 15.0* PI / 180.;
		double theta_mis;
		
		px= (*volumeit).x;
		py= (*volumeit).y;
		volumeit++;
		for (; volumeit!= points.end(); volumeit++){
			s << (*volumeit).x << "\t" << (*volumeit).y<<endl;      
			volume += (py+(*volumeit).y)*(px-(*volumeit).x);
			px= (*volumeit).x;
			py= (*volumeit).y;
			
			theta_mis=mis_ori( handler->ID[1][(int(py+0.5)*m) + int(px+0.5)] );
			if (theta_mis <= theta_ref)	energy += h* gamma_hagb * ( theta_mis / theta_ref) * (1.0 - log( theta_mis / theta_ref));
				else energy += h* gamma_hagb;
		}		
		if (xmin < 0) xmin = 0;
		if (xmax > m) xmax = m;
		if (ymin < 0) ymin = 0;
		if (ymax > m) ymax = m;
		
		stringstream dateiname;
		dateiname << "testpoints_" << id << ".gnu";
		ofstream datei;
		datei.open(dateiname.str());
		datei << s.str();
		datei.close();
		cerr<< "Volume of " << id << "= " << abs(volume)*0.5<< endl;
	}
    return;
}

void LSbox::copy_distances(){
 	//double *buffer = new double[(ymax-ymin)*(xmax-xmin)];
	//distance = buffer;
	distance = new double [(ymax-ymin)*(xmax-xmin)];    
	for (int i = ymin; i < ymax; i++){
		for (int j = xmin; j < xmax; j++){		
			distance[(i-ymin)*(xmax-xmin)+(j-xmin)]=(*domain)[i][j];
			(*domain)[i][j] = INTERIMVAL;
		}
	}
}

void LSbox::copy_distances_to_domain(){
	for (int i = ymin; i < ymax; i++){
		for (int j = xmin; j < xmax; j++){
			(*domain)[i][j]=distance[(i-ymin)*(xmax-xmin)+(j-xmin)];
		}
	}
	delete [] distance;
	distance = NULL;
}

void LSbox::comparison_set_to_domain(LSbox ***ID, int grid_blowup){
	int m = (*handler).get_ngridpoints();
	double h = handler->get_h();
	for (int i = ymin; i < ymax; i++){
		for (int j = xmin; j < xmax; j++){
			if ((i <= grid_blowup) || (m-grid_blowup <= i) || (j <= grid_blowup) || (m-grid_blowup <= j)) {
				(*domain)[i][j] = -DELTA;
			}
			if( abs(distance[(i-ymin)*(xmax-xmin)+(j-xmin)]) < (0.7* DELTA) && (abs((*domain)[i][j]) < ( 0.7 * DELTA)) ) {
// 				 update only in a tube around the n boundary - numerical stability!s
				(*domain)[i][j] = 0.5 * ((*domain)[i][j]-distance[(i-ymin)*(xmax-xmin)+(j-xmin)]);
			}
	
			if ((*domain)[i][j]> 0){
				ID[0][(i*m) + j] = this;
				ID[1][(i*m) + j] = IDLocal[0][(i-ymin)*(xmax-xmin)+(j-xmin)];
				ID[2][(i*m) + j] = IDLocal[1][(i-ymin)*(xmax-xmin)+(j-xmin)];
			}
			else {
	// 			  we are otside the cureent grain, but we has to assure that every gridpoint is updated!!
				ID[0][(i*m) + j] = IDLocal[0][(i-ymin)*(xmax-xmin)+(j-xmin)];
				ID[1][(i*m) + j] = this;
				ID[2][(i*m) + j] = IDLocal[1][(i-ymin)*(xmax-xmin)+(j-xmin)];
			}			
		}
	}
// 	utils::print_2dim_array(distance,ymax-ymin,xmax-xmin);
 	delete [] distance;
	delete [] IDLocal[0];
	delete [] IDLocal[1];
 	distance = NULL;
}


bool LSbox::checkIntersect(LSbox* box2) {    
    if ((xmin > box2->xmax || xmax < box2->xmin || ymin > box2->ymax || ymax < box2->ymin)) return false;
    return true;
}

void LSbox::comparison(const domainCl &domain_copy, int loop){
	domainCl *temp = new domainCl(ymax-ymin, xmax-xmin, 0, -1.0);
	std::vector<LSbox*>::iterator it;
	std::vector<LSbox*>::iterator it_nn;

	if(distance == NULL) {
		IDLocal[0]	=	new LSbox*[(xmax-xmin)*(ymax-ymin)];
		IDLocal[1]	=	new LSbox*[(xmax-xmin)*(ymax-ymin)];
		distance 	= 	new double [(ymax-ymin)*(xmax-xmin)];
		std::fill_n(distance,(ymax-ymin)*(xmax-xmin), -1.0); // Neccesary for the comparison
		std::fill_n(IDLocal[0],(ymax-ymin)*(xmax-xmin), this);
		std::fill_n(IDLocal[1],(ymax-ymin)*(xmax-xmin), this);
	}
	for(it_nn = neighbors_2order.begin(); it_nn != neighbors_2order.end();){		
		if((domain_copy.get_id() == (*(**it_nn).domain).get_id()) ){
			if( ((**it_nn).get_status() == true )) {
			if (checkIntersect(*it_nn)){
				neighbors.push_back(*it_nn);
				int x_min_new, x_max_new, y_min_new, y_max_new;
				
				if(xmin < (**it_nn).xmin) x_min_new = (**it_nn).xmin;
					else x_min_new = xmin;
				
				if(xmax > (**it_nn).xmax) x_max_new = (**it_nn).xmax;
					else x_max_new = xmax;
								
				if(ymin < (**it_nn).ymin) y_min_new = (**it_nn).ymin;
					else y_min_new = ymin;
					
				if(ymax > (**it_nn).ymax) y_max_new = (**it_nn).ymax;
					else y_max_new = ymax;
					
				for (int i = y_min_new; i < y_max_new; i++){
					for (int j = x_min_new; j < x_max_new; j++){
						if( domain_copy[i][j]< 0.7*DELTA){
							if( distance[(i-ymin)*(xmax-xmin)+(j-xmin)] < domain_copy[i][j] ){ 
								distance[(i-ymin)*(xmax-xmin)+(j-xmin)] = domain_copy[i][j];
								if( IDLocal[1][(i-ymin)*(xmax-xmin)+(j-xmin)]!= this ) {
									// we just have found 2 neighbour
									(*temp)[i-ymin][j-xmin] = distance[(i-ymin)*(xmax-xmin)+(j-xmin)];
								}
								IDLocal[1][(i-ymin)*(xmax-xmin)+(j-xmin)] = IDLocal[0][(i-ymin)*(xmax-xmin)+(j-xmin)];							
								IDLocal[0][(i-ymin)*(xmax-xmin)+(j-xmin)] = *it_nn;
							}
							else if( (*temp)[i-ymin][j-xmin] <= domain_copy[i][j] ){ // falls 2 k�rner gleichweit entfernt - geht ohne <= eins verloren!
								(*temp)[i-ymin][j-xmin] = domain_copy[i][j]; 
								IDLocal[1][(i-ymin)*(xmax-xmin)+(j-xmin)] = *it_nn;								  
							}
						}
					}
				}
			}		
		}
		neighbors_2order.erase(it_nn);
	}
		else it_nn++;
		
	}
	// checke schnitt zum randkorn:
	checkIntersect_zero_grain(temp);
	delete temp;
}


void LSbox::checkIntersect_zero_grain(domainCl* temp){
	domainCl* boundary = handler->boundary;
	int grid_blowup = handler->get_grid_blowup();
	vector<LSbox*> mid_in = boundary->getBoxList();
	int m = handler->get_ngridpoints();
	if (!(xmin > mid_in[0]->xmin && xmax < mid_in[0]->xmax && ymin > mid_in[0]->ymin && ymax < mid_in[0]->ymax)){
// 	if (checkIntersect(mid_in[0])){
		for (int i = ymin; i < ymax; i++){
			for (int j = xmin; j < xmax; j++){	
				if ((i <= 2* grid_blowup) || (m-2*grid_blowup <= i) || (j <= 2*grid_blowup) || (m-2*grid_blowup <= j)){
					if(distance[(i-ymin)*(xmax-xmin)+(j-xmin)] < (*boundary)[i][j]){ 
						distance[(i-ymin)*(xmax-xmin)+(j-xmin)] = (*boundary)[i][j];						  
					}
				}
			}
		}
	}
}


void LSbox::add_n2o(){
	neighbors_2order = neighbors;
	vector<LSbox*> ::iterator it_com, it_n2o, it;
	bool just_in;						
	for(it = neighbors.begin(); it != neighbors.end(); it++){		
		for( it_n2o = (**it).neighbors.begin(); it_n2o != (**it).neighbors.end(); it_n2o++){
			just_in = false;
			for(it_com = neighbors_2order.begin(); it_com != neighbors_2order.end(); it_com++){
				if(*it_n2o == *it_com || this == *it_n2o) just_in = true;
			}
			if(just_in == false) neighbors_2order.push_back(*it_n2o);
		}
	}
	neighbors.clear();
}


void LSbox::sweeping (double h, int start_i, int start_j, int direction){
	// directions 0 = y-  // 1= x+ //  2 = y+  //  3 x-  (y- = up // y + down; (0,0)left upper corner)
	int signk=0, signl=0, k=start_i, l=start_j;
	switch (direction) {
		case 0 :   signk = -1; break; 
		case 1 :   signl =  1; break; 
		case 2 :   signk =  1; break; 
		case 3 :   signl = -1; break; 
	}
	int sgn = 1;
	bool zero_found= false;
	do {
		if( (*domain)[k][l] * (*domain)[k+signk][l+signl] < 0 ) zero_found = true;
		
		if(zero_found){
			double candidate = (*domain)[k][l] + (sgn * h);
			k += signk;
			l += signl;
			if ((*domain)[k][l] == INTERIMVAL) 
				(*domain)[k][l] = candidate;
			else if (abs((*domain)[k][l]) < h) {
				// found value close to zero
				// possible candidate, if no zero will be crossed
				candidate = (*domain)[k][l] + (sgn * h);
				// check next
				k += signk;
				l += signl;
				// if also zero -> change sweepvalue
				if ( abs((*domain)[k][l]) < h ) {
					sgn = utils::sgn((*domain)[k][l]);
					candidate += (2 * sgn * h);
				}
			}
			
			if (abs(candidate) < abs((*domain)[k][l])) 
				(*domain)[k][l] = candidate;
		}
		else { k += signk;	l += signl;}
	} 
	while ( ( ymin<=(k+signk)) && ( (k+(signk))<ymax ) && ( xmin<(l+(signl)) ) && ( (l+(signl))<xmax ) );
}


void LSbox::redist_box(double h, int grid_blowup ) {
// 	plot_box(false);
	int m=ymax-ymin;
	int n=xmax-xmin;
	int ii, jj;
	domainCl *temp = new domainCl(m,n,id,-1.0);
	double slope = 1;
	double candidate, i_slope, zero;
	// x-direction forward
	for (int i = ymin; i < ymax; i++) {
		for (int j = xmin; j < xmax-1; j++) {
			ii = i-ymin; jj = j-xmin;
			//check for sign change
			if ((*domain)[i][j] * (*domain)[i][j+1] <= 0.0) {
				// interpolate
				i_slope  = ((*domain)[i][j+1] - (*domain)[i][j]) / h;
				zero = -(*domain)[i][j] / i_slope;
				if ( abs((*temp)[ii][jj]) > abs(zero)) (*temp)[ii][jj] = -zero * utils::sgn(i_slope);
			}
				// calculate new distance candidate and assign if appropriate
			candidate = (*temp)[ii][jj] + (utils::sgn((*domain)[i][j+1]) * h);
			if (abs(candidate) < abs((*temp)[ii][jj+1])) (*temp)[ii][jj+1] = candidate;
		}
	}
	
	// x-direction backward
	for (int i = ymin; i < ymax; i++) {
		for (int j = xmax-1; j > xmin; j--) {
			ii = i-ymin; jj = j-xmin;
			if ((*domain)[i][j] * (*domain)[i][j-1] <= 0.0) {
				// interpolate
				i_slope  = ((*domain)[i][j-1] - (*domain)[i][j]) / h;
				zero = -(*domain)[i][j] / i_slope;
				if ( abs((*temp)[ii][jj]) > abs(zero)) (*temp)[ii][jj] = -zero * utils::sgn(i_slope);
			}
			// calculate new distance candidate and assign if appropriate
			candidate = (*temp)[ii][jj] + (utils::sgn((*domain)[i][j-1]) * h); // replace with the "a"-slope stuff...
			if (abs(candidate) < abs((*temp)[ii][jj-1])) (*temp)[ii][jj-1] = candidate;
		}
// 		(*temp)[ymin][jj-1] = (*temp)[ymin][jj] -h;
	}
// 	
	// y-direction forward
	for (int j = xmin; j < xmax; j++) {		
		for (int i = ymin; i < ymax-1; i++) {	
			ii = i-ymin; jj = j-xmin;
			// check for sign change
			if ((*domain)[i][j] * (*domain)[i+1][j] <= 0.0) {  
				// interpolate
				i_slope  = ((*domain)[i+1][j] - (*domain)[i][j]) / h;
				zero = -(*domain)[i][j] / i_slope;
				if ( abs((*temp)[ii][jj]) > abs(zero)) (*temp)[ii][jj] = -zero * utils::sgn(i_slope);
			}
			// calculate new distance candidate and assign if appropriate
			candidate = (*temp)[ii][jj] + (utils::sgn((*domain)[i+1][j]) * h);
			if (abs(candidate) < abs((*temp)[ii+1][jj])) (*temp)[ii+1][jj] = candidate;
		}

	}
	
	// y-direction backward
	for (int j = xmin; j < xmax; j++) {
		for (int i = ymax-1; i > ymin; i--) {	
			ii = i-ymin;	jj = j-xmin;
				if ((*domain)[i][j] * (*domain)[i-1][j] <= 0.0) {  
					// interpolate
					i_slope  = ((*domain)[i-1][j] - (*domain)[i][j]) / h;
					zero = -(*domain)[i][j] / i_slope;
					if ( abs((*temp)[ii][jj]) > abs(zero)) (*temp)[ii][jj] = -zero * utils::sgn(i_slope);
				}
			// calculate new distance candidate and assign if appropriate
			candidate = (*temp)[ii][jj] + (utils::sgn((*domain)[i-1][j]) * h); // replace with the "a"-slope stuff...
			if (abs(candidate) < abs((*temp)[ii-1][jj])) (*temp)[ii-1][jj] = candidate;
		}
// 		(*temp)[ii-1][xmin] = (*temp)[ii][xmin] - h;
	}
	
	for (int i = ymin; i < ymax; i++){
		for (int j = xmin; j < xmax; j++){
			ii = i-ymin;	jj = j-xmin;
// 			(*domain)[i][j] = (*temp)[ii][jj];
			if (abs((*temp)[i-ymin][j-xmin])< DELTA) (*domain)[i][j] = (*temp)[i-ymin][j-xmin];
			else (*domain)[i][j] = DELTA * utils::sgn((*temp)[i-ymin][j-xmin]);
// 			if (((i <= grid_blowup) && ((m-grid_blowup <= j) || (j <= grid_blowup)) ) || ( (m-grid_blowup <= i) && ((m-grid_blowup <= j) || (j <= grid_blowup)))){
// 				(*domain)[i][j]= (*temp)[i-ymin][j-xmin];
// 			}
		}
	}
	delete temp;
	
}

double LSbox::curvature(int x, int y, double h){	
//     Returns the curvature in the point x,y.
    double phix,phiy, phixx, phiyy,phixy;
    phix	=	((*domain)[y][x+1]-(*domain)[y][x-1])/(2*h);
    phiy	=	((*domain)[y+1][x]-(*domain)[y-1][x])/(2*h);
    phixx	=	((*domain)[y][x-1]-(2*(*domain)[y][x])+(*domain)[y][x+1])/(h*h);
    phiyy	=	((*domain)[y-1][x]-(2*(*domain)[y][x])+(*domain)[y+1][x])/(h*h);
    phixy	=	((*domain)[y-1][x-1]+(*domain)[y+1][x+1]-(*domain)[y-1][x+1]+(*domain)[y+1][x-1])/(4*h*h);	
	double kappa = ((phixx*phiy*phiy)-(2*phix*phiy*phixy)+(phiyy*phix*phix))/sqrt( ((phix*phix)+(phiy*phiy)) * ((phix*phix)+(phiy*phiy)) * ((phix*phix)+(phiy*phiy)) + h);
// 	double N = ((phiyy*phiy*phiy)+(2*phix*phiy*phixy)+(phixx*phix*phix))/sqrt( ((phix*phix)+(phiy*phiy)) * ((phix*phix)+(phiy*phiy)) * ((phix*phix)+(phiy*phiy)) );
	double laplace = (phixx + phiyy);// sqrt((phix*phix) + (phiy*phiy) + h);
//     return(kappa);
	return (laplace);
}      

void LSbox::euler_forward(double dt, double h){
	double kappa, v_n;
	for (int i = ymin+1; i < ymax-1; i++){
		for (int j = xmin+1; j < xmax-1; j++){
			kappa = curvature(i,j,h);
			v_n = kappa;
			(*domain)[i][j] = (dt*v_n) + (*domain)[i][j];
		}
	}
}


void LSbox::plot_box(bool distanceplot){
	cout <<" \nGrain  Info: " << endl;
	cout << " ID :" <<id << endl;
    cout << " xmin, xmax, ymin, ymax :" << xmin << " || "<< xmax << " || " << ymin << " || " << ymax << endl;
    if (distance !=NULL && distanceplot==true) utils::print_2dim_array(distance,ymax-ymin,xmax-xmin);
		else cout << " no distance values in storage!" << endl;
    cout << " Box is in Domain: " << domain->get_id() << endl;
	
	if(neighbors.empty()!=true){
		cout << " List of Neighbors : ";
		vector<LSbox*> :: iterator it;
		for (it = neighbors.begin(); it != neighbors.end(); it++){
			cout << (*it)->getID() << " || ";
		}
		cout << endl;
	} else cout << " neighbors unknown " << endl;
	
	if(neighbors_2order.empty()!=true){
		cout << " List of 2order Neighbors : ";
		vector<LSbox*> :: iterator it;
		for (it = neighbors_2order.begin(); it != neighbors_2order.end(); it++){
			cout << (*it)->getID() << " || ";
		}
		cout << endl;
	} else cout << " neighbors_2order unknown " << endl;
	
}


double LSbox::mis_ori(LSbox* grain_2){
	return misorientationCubic(phi1,PHI,phi2,grain_2->get_phi1(), grain_2->get_PHI(), grain_2->get_phi2());
}






