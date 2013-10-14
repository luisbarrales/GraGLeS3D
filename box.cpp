#include "box.h"

LSbox::LSbox() {}

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


LSbox::LSbox(int id, int nvertex, float* vertices, float phi1, float PHI, float phi2, int grid_blowup, double h, grainhdl* owner) : id(id), phi1(phi1), PHI(PHI), phi2(phi2), nvertices(nvertex), handler(owner){
    
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
LSbox LSbox::distancefunction(int nvertex, float* vertices, int grid_blowup, double h){
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
					if(((grid_blowup < i) && (i < (m- grid_blowup))) && ((grid_blowup < j) && (j < (m- grid_blowup)))) {
						d=abs(d);
// 						char buffer;
// 						cout << " ich wei� nicht ob ich im korn bin";
// 						cin >> buffer;
					}
					else d= abs(d);
		// 			cout << "ID klappt" << endl;
		// 			char buffer;
		//                         cin >> buffer ;
					if(abs(d)< abs(dmin)) dmin=d;
				}
            }

            (*domain)[i][j]= dmin;
        }
	}
	
	for (i=xmin;i<xmax;i++){ // ¸ber gitter iterieren
		j=ymin;
		while( j<ymax && (*domain)[j][i] >= h  ) {
			(*domain)[j][i] = - (*domain)[j][i];	  
			j++;
		} if ( j<ymax ) (*domain)[j][i] = - (*domain)[j][i];
		
		j=ymax-1;
		while( j>=ymin &&  (*domain)[j][i] >= h ) {
			(*domain)[j][i] = - (*domain)[j][i];	  
			j--;
		} if ( j>=ymin ) (*domain)[j][i] = - (*domain)[j][i];
		
	}
	for (j=ymin;j<ymax;j++){ // ¸ber gitter iterieren
		i=xmin;
		while( i<xmax && (*domain)[j][i] >= h   ) {
			(*domain)[j][i] = - (*domain)[j][i];	  
			i++;
		} if ( i<xmax ) (*domain)[j][i] = - (*domain)[j][i];
		
		i=xmax-1;
		while( i>=xmin   && (*domain)[j][i] >= h ) {
			(*domain)[j][i] = - (*domain)[j][i];	  
			i--;
		} if ( i>=xmin ) (*domain)[j][i] = -(*domain)[j][i];
		
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

            (*domain)[i][j]= dmin;
        }
	}
	return(*this);
}

void LSbox::setZeros(double h, int grid_blowup, int loop) {
	if(loop == PRINTNOW) plot_box(false);
    // clear current vector
    zeros.clear();
    
    int first_i, first_j;
    int current_i, current_j;
    int next_i, next_j;
    char direction = 1; // x+
    // directions 0 = y-  //  2 = y+  //  3 x-  (y- = up // y + down; (0,0)left upper corner)

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
				next_i =i+1; 		next_j = j;
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
	xmax = 0; xmin = M; ymax = 0; ymin = M;
    while (newZero) {
		
		double val1 = (*domain)[current_i][current_j];
        double val2 = (*domain)[next_i][next_j];
		
		// interpolate current zero
        double i_slope = (val2 - val1) / h;
        double zero = -val1 / i_slope;
		double bufferVal = -zero; //* i_slope
		int sweep_direction = direction;
// 		if (val1 > 0) sweep_direction=  (sweep_direction+2)%4;
//         zeros.emplace_back(current_i, current_j, bufferVal, sweep_direction);
        // add to zero-list		
		// add left and right point to zerolist, with the direction towards the zero!!
//         bufferVal = (-zero+h);
// 		   zeros.emplace_back(next_i, next_j, bufferVal, (sweep_direction+2)% 4 );
		
		// iterate forward
		current_i = next_i;
		current_j = next_j;
	
        // check for size change
        if (current_j-grid_blowup < xmin) 		xmin = current_j - grid_blowup;
        else if (current_j > xmax-grid_blowup) 	xmax = current_j + grid_blowup;
        if (current_i < ymin+grid_blowup) 		ymin = current_i - grid_blowup;
        else if (current_i > ymax-grid_blowup) 	ymax = current_i + grid_blowup;
		
		if (direction % 2 == 0) {
            if (direction == 0) next_i = current_i-1;
            if (direction == 2) next_i = current_i+1;
            next_j = current_j;
        } else {
            if (direction == 1) next_j = current_j+1;
            if (direction == 3) next_j = current_j-1;
            next_i = current_i;
        }
        // change search directions
		//(1 = left turn; -1 right turn)  
		sgn *= -1;            
		
        // search next zero
		// turn  right
        direction = (direction+sgn+4) %4; 
 
        bool foundnext = false;
        for (int i = 0; i < 3; i++) {
			// who is next?
            if (direction % 2 == 0) {
                if (direction == 0) next_i = current_i-1;
                if (direction == 2) next_i = current_i+1;
                next_j = current_j;
            } else {
                if (direction == 1) next_j = current_j+1;
                if (direction == 3) next_j = current_j-1;
                next_i = current_i;
            }
            
            if ((*domain)[current_i][current_j] * (*domain)[next_i][next_j] < 0) {
                foundnext = true;
                break;
            }            
            
            current_j = next_j; current_i = next_i;
			// turn further if no zero found
			direction = (direction+sgn+4)%4; 
        }
        if (!foundnext) {
            break;
        }
        
        // check if completed round
        if (current_j == first_j && current_i == first_i) {
            newZero = false;
        }
    }
    
    if (xmin < 0) xmin = 0;
    if (xmax > (*domain).get_m()) xmax = (*domain).get_m();
	if (ymin < 0) ymin = 0;
	if (ymax > (*domain).get_n()) ymax = (*domain).get_n();

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

void LSbox::comparison_set_to_domain(LSbox ***ID, int resized_m, int grid_blowup){
	for (int i = ymin; i < ymax; i++){
		for (int j = xmin; j < xmax; j++){
  			(*domain)[i][j]=0.5*((*domain)[i][j]-distance[(i-ymin)*(xmax-xmin)+(j-xmin)]);
			if ((*domain)[i][j]- distance[(i-ymin)*(xmax-xmin)+(j-xmin)] > 0){
			  ID[0][(i*resized_m) + j] = this;
			  ID[1][(i*resized_m) + j] = IDLocal[0][(i-ymin)*(xmax-xmin)+(j-xmin)];
			  ID[2][(i*resized_m) + j] = IDLocal[1][(i-ymin)*(xmax-xmin)+(j-xmin)];
			}
			else {
// 			  we are otside the cureent grain, but we has to assure that every gridpoint is updated!!
			  ID[0][(i*resized_m) + j] = IDLocal[0][(i-ymin)*(xmax-xmin)+(j-xmin)];
			  ID[1][(i*resized_m) + j] = this;
			  ID[2][(i*resized_m) + j] = IDLocal[1][(i-ymin)*(xmax-xmin)+(j-xmin)];
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
	//LSbox*
	domainCl *temp = new domainCl(ymax-ymin, xmax-xmin, 0, 10.0);
	std::vector<LSbox*>::iterator it;
	std::vector<LSbox*>::iterator it_nn;
	double max;
	if(distance == NULL) {
	  IDLocal[0]	=	new LSbox*[(xmax-xmin)*(ymax-ymin)];
	  IDLocal[1]	=	new LSbox*[(xmax-xmin)*(ymax-ymin)];
	  distance 	= 	new double [(ymax-ymin)*(xmax-xmin)];
	  std::fill_n(distance,(ymax-ymin)*(xmax-xmin), INTERIMVAL); // Neccesary for the comparison
	  std::fill_n(IDLocal[0],(ymax-ymin)*(xmax-xmin), this);
	  std::fill_n(IDLocal[1],(ymax-ymin)*(xmax-xmin), this);
	}
	for(it_nn = neighbors_2order.begin(); it_nn != neighbors_2order.end();){		
		if((domain_copy.get_id() == (*(**it_nn).domain).get_id()) ){
			if( ((**it_nn).get_status() == true )) {
			if (checkIntersect(*it_nn)){
			  if ((*(**it_nn).domain).get_id()==domain->get_id()) cout << "GRAINS IN THE SAME DOMAIN !!!!!!!!!!!!!!!!! " << (**it_nn).get_id()<<","<<get_id()<<endl << endl;
  
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
						if(distance[(i-ymin)*(xmax-xmin)+(j-xmin)] < domain_copy[i][j]){ 
							distance[(i-ymin)*(xmax-xmin)+(j-xmin)] = domain_copy[i][j];
							IDLocal[1][(i-ymin)*(xmax-xmin)+(j-xmin)] = IDLocal[0][(i-ymin)*(xmax-xmin)+(j-xmin)];
							IDLocal[0][(i-ymin)*(xmax-xmin)+(j-xmin)] = *it_nn;
						}
						else if( distance[(i-ymin)*(xmax-xmin)+(j-xmin)] < (*temp)[i-ymin][j-xmin]){
							(*temp)[i-ymin][j-xmin] = distance[(i-ymin)*(xmax-xmin)+(j-xmin)];
							IDLocal[1][(i-ymin)*(xmax-xmin)+(j-xmin)] = *it_nn;						
						}
					}
				}
			}
			
		
		}
		neighbors_2order.erase(it_nn);
	}
		else it_nn++;
	}
	
	delete temp;
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




/*void LSbox::free_memory_distance(){
	free (distance);
	distance=NULL;
}*/
	

/*********************************************************************/
/*************** Redistancing Helperfunctions ************************/
/*********************************************************************/


// void LSbox::sweep(pointVal zero, double h){
// 	int sweep_direction = (zero.direction + 2) % 4;
// // 	cout << "direction: " << sweep_direction << endl;
// 	int i=zero.y, l=0, signl = 0; 
// 	int j=zero.x, k=0, signk = 0; 
// 	
// 	if (sweep_direction == -1) {
// 		cout << "start at: "<< i    << " || " << j     << endl;
// 		cout << "min: "     << ymin << " || " << xmin  << endl;
// 		cout << "max: "     << ymax << " || " << xmax  << endl;
// 		cout << "sign: " 	<< signk<< " || " << signl << endl << endl;
// 		return;
// 	}
// 	int sgn = utils::sgn(zero.val);
// 
// 	
// 	// direction always looks into the grain
// 	// to get the opposite for ah negativ boundary value:
// 	// int sweep_direction = 1; // 0 y+  2 y-  1 x+  3 x-  (1 is firstDir)
// 	
// 	switch (sweep_direction) {
// 		case 0 :   signk = 1; break; 
// 		B
// 		case 1 :   signl = 1; break; 
// 		case 2 :   signk = -1; break; 
// 		case 3 :   signl = -1; break; 
// 		
// 	}
// 	if (sgn == 0) sgn = utils::sgn((*domain)[i+signk][j+signl]);
// 	
// 	bool sign_change = false;
// 	
// 	while (  (ymin < (i+k)) && ((i+k) < ymax-1) && (xmin < (j+l)) && ((j+l) < xmax-1) && sign_change == false ){
// 		double candidate = (*domain)[i+k][j+l] + (sgn * h);
// 		k += (signk * 1);
// 		l += (signl * 1);
// 		if ((*domain)[i+k][j+l] == INTERIMVAL )
// 			(*domain)[i+k][j+l] = candidate;
// 		else if (abs(candidate) < abs((*domain)[i+k][j+l])) 
// 			(*domain)[i+k][j+l] = candidate;	
// 		
// 	}
// 	
// }
void LSbox::maximum(const domainCl &A, const domainCl &B){
// 	assert(A.n == B.n);
// 	assert(A.m == B.m);
// 	for (int i = 0; i < m; i++)
// 		for (int j = 0; j < n; j++) {
// 			if (A[i][j] > B[i][j]) (*x[i])[j] = A[i][j];
// 			else (*x[i])[j] = B[i][j];
// 		}
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


void LSbox::redistancing(double h, int grid_blowup /*,std::list<domainCl> distances, double** borderSlopes, double** slopeField*/) {
	
	cout << "Berechne Redist fuer Box: "<< id << endl;	
// 	write zeros to domain:
	vector<pointVal> :: iterator k;
	for (k = zeros.begin(); k != zeros.end(); k++){
		if ( abs((*domain)[(*k).y][(*k).x]) > abs((*k).val)) (*domain)[(*k).y][(*k).x]= (*k).val;	
	}
	for (int j = xmin; j < xmax; j++){
		// sweep_y_forward
		//Startvalues for sweep (ymin, j)
		sweeping(h, ymin, j, 2);
// 		sweep_y_backward
		sweeping(h, ymax-1, j, 0);
	}		
	cout << "y sweeps complete" << endl;
	
	for (int i = ymin; i < ymax; i++){
		sweeping(h, i, xmin, 1);
		sweeping(h, i, xmax-1, 3);
	}
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


void LSbox::find_LevelSet(){
	
	vector<double*> levelset;
	
	int nc = 1;
	double z[1]; z[0]=0.0;		
	// x               ! data matrix column coordinates
	// y               ! data matrix row coordinates
	double x[xmax-xmin];
	double dx =  (handler->get_h());
	double y[ymax-ymin];	
	int gb = (handler)->get_grid_blowup();
	for (int i=ymin; i< ymax; i++) y[i-ymin]=i* dx - gb*dx;
	for (int j=xmin; j< xmax; j++) x[j-xmin]=j* dx - gb*dx;	
	
// 	if(id==2)  utils:: print_2dim_array(x,1,xmax-xmin);
	
	
	int m1,m2,m3,case_value;
	double dmin,dmax,x1,x2,y1,y2;
	int i,j,k,m;
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
	// Fortran and C/C++ arrays are transposed of each #22 0x00000000004287cd in 
	// grainhdl::conrec (this=0x7fffffff9f08) at /home/cm654063/LevelSetProject/grainhdlother, in this case
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
	
	int jub  = xmax;
	int jlb  = xmin;
	int iub  = ymax;
	int ilb  = ymin;
	
	for (j=(jub-2);j>jlb;j--) {
		for (i=ilb+1;i<iub-1;i++) {	
			double temp1,temp2;
			temp1 = min((*domain)[i][j],(*domain)[i][j+1]);
			temp2 = min((*domain)[i+1][j],(*domain)[i+1][j+1]);
			dmin = min(temp1,temp2);
			temp1 = max((*domain)[i][j],(*domain)[i][j+1]);
			temp2 = max((*domain)[i+1][j],(*domain)[i+1][j+1]);
			dmax = max(temp1,temp2);
// 			cout<< i<<" || " <<j << endl;
			
			if ( dmax >= z[0] && dmin <= z[nc-1] ) {
				for (k=0;k<nc;k++) {
					if (z[k]>=dmin&&z[k]<=dmax) {
						for (m=4;m>=0;m--) {
							if (m>0) {
								//=============================================================
								// The indexing of im and jm should be noted as it has to
								// start from zero
								//=============================================================
								h[m] = (*domain)[i+im[m-1]][j+jm[m-1]]-z[k];
								xh[m] = x[i+im[m-1]];
								yh[m] = y[j+jm[m-1]];							
							} 
							else {
								h[0] = 0.25*(h[1]+h[2]+h[3]+h[4]);
								xh[0]=0.5*(x[i]+x[i+1]);
								yh[0]=0.5*(y[j]+y[j+1]);							
							}
							if (h[m]>0.0) {	sh[m] = 1;} 
							else if (h[m] < 0.0) {	sh[m] = -1; } 
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
								double line[4] = {x1,y1,x2,y2};
								// 	if(id==2)  {utils:: print_2dim_array(line,1,4);
								printf("%d\t%d\t%f\t%f\t%f\t%f\n", i,j ,x1, y1, x2, y2);
								levelset.push_back(line);
							}
						}
					}
				}
			}
		}
		domainCl temp(2* levelset.size(),2);
		for(int i=0; i< levelset.size(); i++){
			
			temp[2*i][0]= levelset[i][0];
			temp[2*i][1]= levelset[i][1];
			temp[(2*i)+1][0]= levelset[i][2];
			temp[(2*i)+1][1]= levelset[i][3];
		}
		temp.save_domainCl("grain.gnu");
	}
	// compute cell volume
	
	compute_volume(levelset);
	char buffer;
	cin >> buffer;
}

void LSbox::compute_volume(vector<double*> &levelset){
	float sum=.0;
	cout << " size: "  << (levelset).size() << endl;
	int N = (levelset).size();
	
	for( int i=0;i< N ;i++ ){
		sum += ( ((levelset)[i][0] + (levelset)[i][2] ) - ((levelset)[i][3] - (levelset)[i][1] )); 
		//Gau�sche Trapezformel
    }
    volume = 0.5* abs(sum);	
	
	cout <<"id: " << id <<"   vol: "<< volume << endl;
}
