#include "box.h"

LSbox::LSbox() {
    
}

LSbox::LSbox(int aID, voro::voronoicell_neighbor& c, double *part_pos, int grid_blowup, double h):id(aID) {
    
    // determine size of grain
    xmax = 0; xmin = 1000; ymax = 0; ymin = 1000;
    vektor x1(2), x2(2);
    vector<double> vv;
	c.vertices(part_pos[3*id],part_pos[3*id+1],part_pos[3*id+2],vv);
    
    for(int ii=0;ii<c.p;ii++) {
        for(int jj=0;jj<c.nu[ii];jj++) {
            
            int k=c.ed[ii][jj];
            x1[0]=vv[3*ii];x1[1]=vv[3*ii+1];
            x2[0]=vv[3*k]; x2[1]=vv[3*k+1];
            
            if (x1[0]/h < xmin) xmin = x1[0]/h;
            if (x2[0]/h < xmin) xmin = x2[0]/h;
            if (x1[1]/h < ymin) ymin = x1[1]/h;
            if (x2[1]/h < ymin) ymin = x2[1]/h;
            if (x1[0]/h > xmax) xmax = x1[0]/h;
            if (x2[0]/h > xmax) xmax = x2[0]/h;
            if (x1[1]/h > ymax) ymax = x1[1]/h;
            if (x2[1]/h > ymax) ymax = x2[1]/h;
            
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
	
    
	for (i=xmin;i<xmax;i++){ // ¸ber gitter iterieren
        for (j=ymin;j< ymax;j++){
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
    bool zero_found = false;
	cout << "grain: " << id << endl << "Boxabmessung: " << endl;;
	
	cout << xmin << " || " << xmax << endl;
	cout << ymin << " || " << ymax << endl;
    
    int firstx, firsty;
    int currentx=0, currenty=0;
    char direction = 1; // 0 y+  2 y-  1 x+  3 x-  (1 is firstDir)
    
    int dist = ymax - ymin;
    int y = ymin+ int(dist/2);
    // look for zero in row y
    for (int j = xmin; j < xmax; j++) {
// 		cout << y <<" || "<< j << " || value: " << (*domain)[y][j] <<endl;
        if ((*domain)[y][j] * (*domain)[y][j+1] <= 0) {
            firstx = j; firsty = y;
            currentx = j; currenty = y;
			cout << "found first zero: " << currentx << " || " << currenty << endl;
            zero_found= true;
			break;
        }
    }
    
    if (!zero_found) {
		cout << "error: no grain in box: " << id << endl << endl;
		return;
	}
    // begin zero-tracking and interpolations
    bool newZero = true;

    cout << "suche neue boxgroesse" << endl;
    while (newZero) {
        // interpolate current zero
        int nextx, nexty;
        if (direction % 2 == 0) {
            if (direction == 0) nexty = currenty+1;
            else nexty = currenty-1;
            nextx = currentx;
        } else {
            if (direction == 1) nextx = currentx+1;
            else nextx = currentx-1;
            nexty = currenty;
        }

        double val1 = (*domain)[currenty][currentx];
        double val2 = (*domain)[nexty][nextx];

        double i_slope = (val2 - val1) / h;
        double zero = -val1 / i_slope;
        

        
        // add to zero-list
        double bufferVal = -zero * i_slope;
        zeros.emplace_back(currentx, currenty, bufferVal);
        bufferVal = (h-zero) * i_slope;
        zeros.emplace_back(nextx, nexty, bufferVal);
        
// 		cout << currentx << " - " << currenty << endl;
		
        // check for size change
        if (currentx-grid_blowup < xmin) xmin = currentx-grid_blowup;
        else if (currentx > xmax-grid_blowup) xmax = currentx+grid_blowup;
        if (currenty < ymin+grid_blowup) ymin = currenty-grid_blowup;
        else if (currenty > ymax-grid_blowup) ymax = currenty+grid_blowup;
        
        // find next zero
        direction = (direction-1)%4; //left turn
 
        bool foundnext = false;
        for (int i = 0; i < 3; i++) {
            
            int nextx, nexty;
            if (direction % 2 == 0) {
                if (direction == 0) nexty = currenty+1;
                else nexty = currenty-1;
                nextx = currentx;
            } else {
                if (direction == 1) nextx = currentx+1;
                else nextx = currentx-1;
                nexty = currenty;
            }
            
            if ((*domain)[currenty][currentx] * (*domain)[nexty][nextx] < 0) {
                foundnext = true;
                break;
            }
            
            direction = (direction+1)%4; // right turn
            currentx = nextx; currenty = nexty;
        }
        if (!foundnext) {
            break;
        }
        
        // check if completed round
        if (currentx == firstx && currenty == firsty && direction == 1) {
            newZero = false;
        }
    }
    
    if (xmin < 0) xmin = 0;
    if (xmax > (*domain).get_m()) xmax = (*domain).get_m();
	if (ymin < 0) ymin = 0;
	if (ymax > (*domain).get_n()) ymax = (*domain).get_n();
	cout <<"neue Abmessungen : " << endl;
	cout << xmin << " || " << xmax << endl;
	cout << ymin << " || " << ymax << endl << endl;
}

bool LSbox::checkIntersect(LSbox* box2) {
    
    //if (xmin <= box2->xmax && xmax >= box2->xmin && ymin <= box2->ymax && ymax >= box2->ymin) return true;
    if ((xmin > box2->xmax || xmax < box2->xmin || ymin > box2->ymax || ymax < box2->ymin)) return false;
    
    return true;
}

void LSbox::redistancing(double h, int grid_blowup /*,std::list<matrix> distances, double** borderSlopes, double** slopeField*/) {
	int m = ymax - ymin;
    int n = xmax - xmin;
	matrix *temp = new matrix(m,n,id,-INTERIMVAL);

    double limiter = INTERIMVAL;
    double slope = 1;
	
	cout << "Berechne Redist fuer Box: "<< id << endl;
	
	// write zeros to domain:
	vector<pointVal> :: iterator k;
	for (k = zeros.begin(); k != zeros.end(); k++){
		(*temp)[(*k).y][(*k).x]= (*k).val;	
	}
	(*temp).redistancing(h, grid_blowup); // ruft rististancing aus Matrixklasse auf
// 	(*temp).redistancing_2(h, grid_blowup);

	// copy temp to domain
	int ii,jj,i,j;
    for (ii=0, i = xmin; i < xmax; i++, ii++) {
		for (jj=0, j = ymin; j < ymax; j++, jj++) {
			(*domain)[i][j]= (*temp)[ii][jj];			
		}
    }
    cout << "Box in Domain kopiert -- success"<< endl;
	delete temp;
}
