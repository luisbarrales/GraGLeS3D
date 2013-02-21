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
                        if((0. < lambda) && (lambda < 1.)) 	d= (p-(a+(u*lambda))).laenge();
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

void LSbox::setZeros(double h) {
    
    // clear current vector
    zeros.clear();
    
    int firstx, firsty;
    int currentx, currenty;
    char direction = 1; // 0 y+  2 y-  1 x+  3 x-  (1 is firstDir)
    
    // find first zero
    bool success = false;
    int dist = ymax - ymin;
    int k = 2;
    while (!success) {
        for (int i = 0; i < k; i+=2) {
            int y = ((1+i)*dist) / k;
            // look for zero in row y
            for (int j = xmin; j < xmax; j++) {
                if ((*domain)[y][j] * (*domain)[y][j+1] <= 0) {
                    firstx = j; firsty = y;
                    currentx = j; currenty = y;
                    success = true;
                    break;
                }
            }
        }
        k*=2;
    }
    
    // begin zero-tracking and interpolations
    bool newZero = true;
    
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
        
        // check for size change
        if (currentx < xmin) xmin = currentx;
        else if (currentx > xmax) xmax = currentx;
        if (currenty < ymin) ymin = currenty;
        else if (currenty > ymax) ymax = currenty;
        
        // find next zero
        direction = (direction-1)%4; //left turn
        
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
            
            if ((*domain)[currenty][currentx] * (*domain)[nexty][nextx] < 0) break;
            
            direction = (direction+1)%4; // right turn
            currentx = nextx; currenty = nexty;
        }
        
        // check if completed round
        if (currentx == firstx && currenty == firsty && direction == 1) newZero = false;
    }

}

bool LSbox::checkIntersect(LSbox* box2) {
    
    //if (xmin <= box2->xmax && xmax >= box2->xmin && ymin <= box2->ymax && ymax >= box2->ymin) return true;
    if ((xmin > box2->xmax || xmax < box2->xmin || ymin > box2->ymax || ymax < box2->ymin)) return false;
    
    return true;
}
