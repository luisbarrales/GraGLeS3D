#include "box.h"

LSbox::LSbox() {
    
}

LSbox::LSbox(int aID, voro::voronoicell_neighbor& c, double *part_pos, int grid_blowup, double h):id(aID) {
    
    // determine size of grain
    xmax = 0; xmin = M; ymax = 0; ymin = M;
    vektor x1(2), x2(2);
    vector<double> vv;
	c.vertices(part_pos[3*id],part_pos[3*id+1],part_pos[3*id+2],vv);
    distance = (double*) calloc ((ymax-ymin)*(xmax-xmin),sizeof(double));
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
            
            if (x1[0]/h > ymax) ymax = (x1[0]/h) + 1;
            if (x2[0]/h > ymax) ymax = (x2[0]/h) + 1;
			
			// check for "Spalten" Minima/Maxima
            if (x1[1]/h < xmin) xmin = x1[1]/h;
            if (x2[1]/h < xmin) xmin = x2[1]/h;
            
            if (x1[1]/h > xmax) xmax = (x1[1]/h) + 1;
            if (x2[1]/h > xmax) xmax = (x2[1]/h) + 1;
          
            
        }
    }
    
    xmax += 2*grid_blowup;
	ymax += 2*grid_blowup;
        
    cout << "made a new box: xmin="<<xmin<< " xmax="<<xmax <<" ymin="<<ymin << " ymax="<<ymax<<endl;
       
}

LSbox::~LSbox() {
}

void LSbox::setDomain(matrix* aDomain) {
    domain = aDomain;
}

int LSbox::getID() {
    return id;
}


LSbox LSbox::distancefunction(voro::voronoicell_neighbor& c, int *ID_mat, double *part_pos, int grid_blowup, double h){
    int i,j,k;
	
		
	double d, dmin,lambda;
	int m=domain->get_m();
	int n=domain->get_n();

	vektor u(2), a(2), p(2), x1(2), x2(2);
    vector<double> vv;
	c.vertices(part_pos[3*id],part_pos[3*id+1],part_pos[3*id+2],vv);
	double domain_vertices[] = {0.,0.,1.,0.,1.,1.,0.,1.,0.,0.}; // array of vertices to loop over
	    
	for (i=ymin;i<ymax;i++){ // ¸ber gitter iterieren
        for (j=xmin;j< xmax;j++){
            dmin=1000.;
            p[0]=(i-grid_blowup)*h; p[1]=(j-grid_blowup)*h;            
            
            for(int ii=0;ii<c.p;ii++) {
                for(int jj=0;jj<c.nu[ii];jj++) {
                    
                    k=c.ed[ii][jj];
                    
                    x1[0]=vv[3*ii];x1[1]=vv[3*ii+1];
                    x2[0]=vv[3*k]; x2[1]=vv[3*k+1];
                    
                    if (x1!=x2){
                        a=x1;
                        u = x2-x1;
                        lambda=((p-a)*u)/(u*u);
                        
                        if(lambda <= 0.) 					    d= (p-x1).laenge();
                        if((0. < lambda) && (lambda < 1.)) 		d= (p-(a+(u*lambda))).laenge();
                        if(lambda >= 1.) 					    d= (p-x2).laenge();
                        
                        if(id==ID_mat[i*m +j] && ((grid_blowup < i) && (i < (m- grid_blowup))) && ((grid_blowup < j) && (j < (m- grid_blowup)))) d=abs(d);
                        else d= -abs(d);
                        
                        if(abs(d)< abs(dmin)) dmin=d;
                    }
                }
            }

            (*domain)[i][j]= dmin;
        }
	}
	return(*this);
}

void LSbox::setZeros(double h, int grid_blowup) {

    // clear current vector
    zeros.clear();
    
    int first_i, first_j;
    int current_i, current_j;
	int next_i, next_j;
    char direction = 1; // x+
    // directions 0 = y-  //  2 = y+  //  3 x-  (y- = up // y + down; (0,0)left upper corner)
    int dist = ymax - ymin;

    
    int i = ymin+ int(dist/2);
    // look for zero in row y
    for (int j = xmin; j < xmax; j++) {
        if ((*domain)[i][j] * (*domain)[i][j+1] <= 0) {
            first_i = i; 	first_j = j; 
            current_i = i; 	current_j = j;
			next_i =i; 		next_j = j+1;
			break;
        }
    }
// 	cout << "neue Abmessungen : " << endl;
// 	cout << xmin << " || " << xmax << endl;
// 	cout << ymin << " || " << ymax << endl << endl;
//     cout << "first zero: " << first_i << " || " << first_j << endl;

    // begin zero-tracking and interpolations
    bool newZero = true;
	int sgn = -1; //(1 = left turn; -1 right turn)  

    while (newZero) {
		
		double val1 = (*domain)[current_i][current_j];
        double val2 = (*domain)[next_i][next_j];
		
		// interpolate current zero
        double i_slope = (val2 - val1) / h;
        double zero = -val1 / i_slope;
		double bufferVal = -zero; //* i_slope
		int sweep_direction = direction;
		if (val1 > 0) sweep_direction=  (sweep_direction+2)%4;
        zeros.emplace_back(current_i, current_j, bufferVal, sweep_direction);
        // add to zero-list		
		// add left and right point to zerolist, with the direction towards the zero!!
        bufferVal = (-zero+h);
        zeros.emplace_back(next_i, next_j, bufferVal, (sweep_direction+2)% 4 );
		
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
        if (current_j == first_j && current_i == first_i && direction == 1) {
            newZero = false;
        }
    }
    
    if (xmin < 0) xmin = 0;
    if (xmax > (*domain).get_m()) xmax = (*domain).get_m();
	if (ymin < 0) ymin = 0;
	if (ymax > (*domain).get_n()) ymax = (*domain).get_n();

// 	cout << "neue Abmessungen : " << endl;
// 	cout << xmin << " || " << xmax << endl;
// 	cout << ymin << " || " << ymax << endl << endl;

}

void LSbox::copy_distances(){
	for (int j = xmin; j < xmax; j++){
		for (int i = ymin; i < ymax; i++){
			distance[i*(xmax-xmin)+j]=(*domain)[i][j];
			(*domain)[i][j] = INTERIMVAL;
		}
	}
}

void LSbox::copy_distances_to_domain(){
	for (int j = xmin; j < xmax; j++){
		for (int i = ymin; i < ymax; i++){
			(*domain)[i][j]=distance[i*(xmax-xmin)+j];
		}
	}
}


bool LSbox::checkIntersect(LSbox* box2) {    
    //if (xmin <= box2->xmax && xmax >= box2->xmin && ymin <= box2->ymax && ymax >= box2->ymin) return true;
    if ((xmin > box2->xmax || xmax < box2->xmin || ymin > box2->ymax || ymax < box2->ymin)) return false;
    return true;
}




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

void LSbox::sweeping (double h, int start_i, int start_j, int direction){
	// directions 0 = y-  // 1= x+ //  2 = y+  //  3 x-  (y- = up // y + down; (0,0)left upper corner)
	int signk=0, signl=0, k=start_i, l=start_j;
	switch (direction) {
		case 0 :   signk = -1; break; 
		case 1 :   signl = 1; break; 
		case 2 :   signk = 1; break; 
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


void LSbox::redistancing(double h, int grid_blowup /*,std::list<matrix> distances, double** borderSlopes, double** slopeField*/) {
	
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
