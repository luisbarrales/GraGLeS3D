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
    
	cout << "grain: " << id << endl;
    
    int firstx, firsty;
    int currentx=0, currenty=0;
    char direction = 1; // 0 y+  2 y-  1 x+  3 x-  (1 is firstDir)
    
    int dist = ymax - ymin;
    int y = ymin+ (dist/2);
    cout << dist << "  " << y << endl;
    // look for zero in row y
    for (int j = xmin; j < xmax; j++) {
        if ((*domain)[y][j] * (*domain)[y][j+1] <= 0) {
            firstx = j; firsty = y;
            currentx = j; currenty = y;
            cout << "found" << endl;
            break;
        }
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
}

bool LSbox::checkIntersect(LSbox* box2) {
    
    //if (xmin <= box2->xmax && xmax >= box2->xmin && ymin <= box2->ymax && ymax >= box2->ymin) return true;
    if ((xmin > box2->xmax || xmax < box2->xmin || ymin > box2->ymax || ymax < box2->ymin)) return false;
    
    return true;
}

void LSbox::redistancing(/*double h, int grid_blowup, std::list<matrix> distances, double** borderSlopes, double** slopeField*/) {

//     matrix *temp = new matrix(m,n,id);

    double limiter = INTERIMVAL;
    double slope = 1;
	
	   // THIS IS THE VERSION USING SIGN CHANGES TO GET THE SLOPES
    //slope=1;
    // x-direction forward
	cout << "berechne box: "<< id << endl;
	
    for (int i = xmin; i < xmax; i++) {
		for (int j = ymin; j < ymax; j++) {
			(*domain)[i][j]= (double)id;			
		}
//             if (j==0) (*temp)[i][j] = -limiter;
//             (*temp)[i][j+1] = limiter * utils::sgn((*this)[i][j+1]); // set temp to limiter initially
//             
//             // check for sign change
//             if ((*this)[i][j] * (*this)[i][j+1] < 0.0) {
//                 // find grain with minimal distance to [i][j]
//                 int rightID = (*this).id;
//                 int leftID = minimumInPoint(distances, i, j, rightID);
//                 slope = borderSlopes[leftID][rightID];
//                 
//                 if (slope == 0) slope = 1;
//                 
//                 // interpolate
//                 double i_slope  = ((*this)[i][j+1] - (*this)[i][j]) / h;
//                 double zero = -(*this)[i][j] / i_slope;
//                 if ( abs((*temp)[i][j]) > abs(-zero)) (*temp)[i][j] = -zero * utils::sgn(i_slope);
// 			}
//             // calculate new distance candidate and assign if appropriate
// 			double candidate = (*temp)[i][j] + (utils::sgn((*this)[i][j+1]) * h * slope); 
// 			if (abs(candidate) < abs((*temp)[i][j+1])) (*temp)[i][j+1] = candidate;
//         }
    }
    
//     // y-direction forward
//     for (int j = ymin; j <= ymax; j++) {
//         slope = 1;
//         for (int i = xmin; i <= xmax-1; i++) {
//             
//             // check for sign change
//             if ((*this)[i][j] * (*this)[i+1][j] < 0.0) {
//                 // find grain with minimal distance to [i][j]
//                 int bottomID = (*this).id;
//                 int topID = minimumInPoint(distances, i, j, bottomID);
//                 slope = borderSlopes[topID][bottomID];
//                 
//                 if (slope == 0) slope = 1;
//                 
//                 // interpolate
//                 double i_slope  = ((*this)[i+1][j] - (*this)[i][j]) / h;
//                 double zero = -(*this)[i][j] / i_slope;
//                 if ( abs((*temp)[i][j]) > abs(-zero)) (*temp)[i][j] = -zero * utils::sgn(i_slope);
// 			}
//             // calculate new distance candidate and assign if appropriate
// 			double candidate = (*temp)[i][j] + (utils::sgn((*this)[i+1][j]) * h * slope);
// 			if (abs(candidate) < abs((*temp)[i+1][j])) (*temp)[i+1][j] = candidate;
//         }
//     }
//     
//     // x-direction backward
//     for (int i = xmin; i <= xmax; i++) {
//         slope = 1;
//         for (int j = ymax-1; j >= ymin; j--) {
//             
//             // check for sign change
//             if ((*this)[i][j] * (*this)[i][j-1] < 0.0) {
//                 // find grain with minimal distance to [i][j]
//                 int leftID = (*this).id;
//                 int rightID = minimumInPoint(distances, i, j, leftID);
//                 slope = borderSlopes[leftID][rightID];
//                 
//                 if (slope == 0) slope = 1;
//             }
//             
//             
//             
//             // calculate new distance candidate and assign if appropriate
// 			double candidate = (*temp)[i][j] + (utils::sgn((*this)[i][j-1]) * h * slope); // replace with the "a"-slope stuff...
// 			if (abs(candidate) < abs((*temp)[i][j-1])) (*temp)[i][j-1] = candidate;
//         }
//     }
//     
//     
//     // y-direction backward
//     for (int j = ymin; j <= ymax; j++) {
//         slope = 1;
//         for (int i = xmax-1; i >= xmin; i--) {
//             
//             // check for sign change
//             if ((*this)[i][j] * (*this)[i-1][j] < 0.0) {
//                 // find grain with minimal distance to [i][j]
//                 int topID = (*this).id;
//                 int bottomID = minimumInPoint(distances, i, j, topID);
//                 slope = borderSlopes[topID][bottomID];
//                 
//                 if (slope == 0) slope = 1;
//             }
//             
//             // calculate new distance candidate and assign if appropriate
// 			double candidate = (*temp)[i][j] + (utils::sgn((*this)[i-1][j]) * h * slope); // replace with the "a"-slope stuff...
// 			if (abs(candidate) < abs((*temp)[i-1][j])) (*temp)[i-1][j] = candidate;
//         }
//     }
//     
//     
 }
