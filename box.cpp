#include "box.h"



LSbox::LSbox() {}

LSbox::LSbox(int id, int xmin, int xmax, int ymin, int ymax, double phi1, double PHI, double phi2): id(id), xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), phi1(phi1), PHI(PHI), phi2(phi2) {
  
  
  	IDLocal.resize((xmax-xmin)*(ymax-ymin));	
  
}

LSbox::LSbox(int aID, voro::voronoicell_neighbor& c, double *part_pos, grainhdl* owner) : id(aID), phi1(0), PHI(0), phi2(0), nvertices(0), handler(owner) {
    
	int grid_blowup = owner->get_grid_blowup(); 
	double h = owner->get_h();
    // determine size of grain
    xmax = 0; xmin = handler->get_ngridpoints(); 
    ymax = 0; ymin = xmin;
	
    vektor x1(2), x2(2);
    vector<double> vv;
    exist = true;
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
		
	IDLocal.resize((xmax-xmin)*(ymax-ymin));	
	distance_current.resize(xmax-xmin * ymax-ymin);
	distance_new.resize(xmax-xmin * ymax-ymin);
	cout << "made a new box: xmin="<<xmin<< " xmax="<<xmax <<" ymin="<<ymin << " ymax="<<ymax<<endl;
}


LSbox::LSbox(int id, int nvertex, double* vertices, double phi1, double PHI, double phi2, grainhdl* owner) : id(id), phi1(phi1), PHI(PHI), phi2(phi2), nvertices(nvertex), handler(owner){
    
	int grid_blowup = owner->get_grid_blowup(); 
	double h = owner->get_h();
    // determine size of grain
    xmax = 0; xmin = handler->get_ngridpoints(); 
	ymax = 0; ymin = xmin;
	    
    vektor x1(2), x2(2);
    exist = true;

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
	IDLocal.resize((xmax-xmin)*(ymax-ymin));
    
	distance_current.resize(xmax-xmin * ymax-ymin);
	distance_new.resize(xmax-xmin * ymax-ymin);
	
	cout << "made a new box: xmin="<<xmin<< " xmax="<<xmax <<" ymin="<<ymin << " ymax="<<ymax<<endl;

}





LSbox::~LSbox() {
}


int LSbox::getID() {
    return id;
}

void LSbox::distancefunction(int nvertex, double* vertices ){
// 	plot_box(false);
	int grid_blowup = owner->get_grid_blowup(); 
	double h = owner->get_h();
	int i,j,k;
	double d, dmin,lambda;
	vektor u(2), a(2), p(2), x1(2), x2(2);

	for (i=ymin;i<ymax;i++){ // ¸ber gitter iterieren
	  for (j=xmin;j<xmax;j++){
            dmin = 1000.;
            p[0] = (i-grid_blowup)*h; p[1] = (j-grid_blowup)*h;            
            
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
			if (abs(dmin) < DELTA) distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]= dmin;
			else distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]= DELTA * utils::sgn(dmin);
        }
	}
	int count = 0;
	for (i=xmin;i<xmax;i++){ // ¸ber gitter iterieren
		j=ymin;
		count = 0;
		while( j<ymax  && count < 1) {
			distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)] = - abs(distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]]);
			if ( -(distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]) <=  h ) count++;
			j++;
		} 		
		j=ymax-1;
		count =0;
		while( j>=ymin && count < 1) {
			distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)] = - abs(distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]);	
			if ( -(distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]) <= h ) count++;
			j--;
		} 
	}

	for (j=ymin;j<ymax;j++){ // ¸ber gitter iterieren
		i=xmin;
		count = 0;
		while( i<xmax  && count < 1 ) {
			distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)] = - abs(distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]);
			if ( -(distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]) <=  h ) count++;
			i++;
		} 
		i=xmax-1;
		count =0;
		while( i>=xmin   && count < 1  ) {			
			distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)] = - abs(distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]);	
			if ( -(distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]) <=  h ) count++;
			i--;
		} 
	}
}


void LSbox::distancefunction(voro::voronoicell_neighbor& c, double *part_pos){
	int grid_blowup = owner->get_grid_blowup(); 
	double h = owner->get_h();
	
	int i,j,k;
	double d, dmin,lambda;

	vektor u(2), a(2), p(2), x1(2), x2(2);
	vector<double> vv;
	c.vertices (part_pos[3*(id-1)],part_pos[3*(id-1)+1],part_pos[3*(id-1)+2],vv);
	double domain_vertices[] = {0.,0.,1.,0.,1.,1.,0.,1.,0.,0.}; // array of vertices to loop over
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
                        if(lambda >= 1.) d = (p-x2).laenge();
						if(abs(d)< abs(dmin)) dmin=d;
                    }
                }
            }
// 			(*domain)[i][j]= dmin;
            if (abs(dmin) < DELTA) distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]= dmin;
            else distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]= DELTA * utils::sgn(dmin);
        }
	}
		int count = 0;
	for (i=xmin;i<xmax;i++){ // ¸ber gitter iterieren
		j=ymin;
		count = 0;
		while( j<ymax  && count < 1) {
			distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)] = - abs(distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]]);
			if ( -(distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]) <=  h ) count++;
			j++;
		} 		
		j=ymax-1;
		count =0;
		while( j>=ymin && count < 1) {
			distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)] = - abs(distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]);	
			if ( -(distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]) <= h ) count++;
			j--;
		} 
	}

	for (j=ymin;j<ymax;j++){ // ¸ber gitter iterieren
		i=xmin;
		count = 0;
		while( i<xmax  && count < 1 ) {
			distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)] = - abs(distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]);
			if ( -(distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]) <=  h ) count++;
			i++;
		} 
		i=xmax-1;
		count =0;
		while( i>=xmin   && count < 1  ) {			
			distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)] = - abs(distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]);	
			if ( -(distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]) <=  h ) count++;
			i--;
		} 
	}
}

void domainCl::convolution(const double dt, double *ST, LSbox ***ID, domainCl &ref, LSbox* zeroBox, int grid_blowup, weightmap* my_weights){
	
	int n= grainhdl->get_ngridpoints();
	
	fftTemp = (fftw_complex*) fftw_malloc(n*(floor(n/2)+1)*sizeof(fftw_complex));

	makeFFTPlans(val, fftTemp,&fwdPlan,&bwdPlan);
	conv_generator(val,fftTemp,fwdPlan,bwdPlan,dt);

	fftw_destroy_plan(fwdPlan);
	fftw_destroy_plan(bwdPlan);

	fftw_free (fftTemp);
	/*********************************************************************************/
	// Velocity Corrector Step: 
	/*********************************************************************************/
	// hier soll energycorrection gerechnet werden.
	// in der domainCl steht die urspr�nglich distanzfunktion, in dem arry die gefaltete
	if(!ISOTROPIC){
		double rad =  DELTA* 0.7; // radius in dem ein drag wirkt
		double weight;
// 		int* rep = new int[3];
	  vector<LSbox*>::iterator it;
	  int intersec_xmin, intersec_xmax, intersec_ymin, intersec_ymax;
	  
	  if (old_xmin < xmin) 		
		intersec_xmin = xmin;
	  else 	intersec_xmin = old_xmin;
	  
	  if (old_ymin < ymin) 		
		intersec_ymin = ymin;
	  else 	intersec_ymin = old_ymin;
	  
	  if (old_xmax > xmax) 	
		intersec_xmax= xmax;
	  else  intersec_xmax = old_xmax;
	  
	  if (old_ymax > ymax) 	
		intersec_ymax= ymax;
	  else  intersec_ymax = old_ymax;
	  
	  
	  for (int i = intersec_ymin; i < intersec_ymax; i++){
		for (int j = intersec_xmin; j < intersec_xmax; j++) {
// 		  if(  ID[0][i*m +j] == zeroBox ) continue; ///////// unneccessary
// 		  if (!( ID[0][i*m +j] != ID[1][i*m +j] && ID[1][i*m +j] != ID[2][i*m +j] && ID[0][i*m +j] != ID[2][i*m +j])) continue;
// 		    if ( rad < abs(ref[i][j]) ) continue;
				
				ID[(i-old_xmin)*(old_xmax-old_xmin) + (j-old_ymin)] //== (**it).get_id() || ID[1][i*m +j]->get_id() == (**it).get_id() || ID[2][i*m +j]->get_id() == (**it).get_id() )
// 				{	
					
					weight = local_weights.load_weights(ID[(i-old_xmin)*(old_xmax-old_xmin) + (j-old_ymin)]);
// 					weight = ( 1-abs(rad - abs(ref[i][j])) ) * weight;		nur sinnvoll um einen drag zu simulieren			
				
					(*this)[i][j] = ref[i][j] + (((*this)[i][j] -ref[i][j]) * weight);
					(**it).(*(IDLocal[i])[j]).clear(); // removes all ID's add that that gridpoint -> neccessary for update in the comparison
// 				}
				else
				{
					
// 					weight = (*my_weights).load_weights(m,ST)
					cout << "ID not found! " << (**it).get_id() << endl;
					cout << (*this)[i][j] << "   DELTA: " << DELTA << "  h=  "<< owner->get_h() <<endl;
// 					cout << ID[0][i*m +j]->get_id()  << "  "<< ID[1][i*m +j]->get_id() << "  "<< ID[2][i*m +j]->get_id() <<endl;
// 					cout << ID[0][i*m +j]->domain->entry(i,j)  << "  "<< ID[1][i*m +j]->domain->entry(i,j) << "  "<< ID[2][i*m +j]->domain->entry(i,j) <<endl;
					(**it).plot_box(true);
// 					ID[0][i*m +j]->plot_box(true);
					char buffer;
// 					owner->my_weights->plot_weightmap(n,ID, ST, zeroBox);
					cin >> buffer;
				}						  
				}
			}
		}
	}
}



void LSbox::find_contour() {
	int grid_blowup = owner->get_grid_blowup(); 
	double h = owner->get_h();
	int loop = owner->loop;
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
        if (distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]  * distance_current[(i-ymin)*(xmax-xmin)+(j-xmin+1)]  <= 0) {
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
			if (distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]  * distance_current[(i-ymin+1)*(xmax-xmin)+(j-xmin)]  <= 0) {
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
	old_xmin = xmin; 
	old_xmax = xmax; 
	old_ymin = ymin; 
	old_ymax = ymax;
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
		if (direction == 0) 	   {next_i = current_i-1;next_j = current_j;}
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
		  if ( distance_current[(current_i-ymin)*(xmax-xmin)+(current_j-xmin)] * distance_current[(next_i-ymin)*(xmax-xmin)+(next_j-xmin)] <= 0) {
			  foundnext = true;
			  if( distance_current[(current_i-ymin)*(xmax-xmin)+(current_j-xmin)] * distance_current[(next_i-ymin)*(xmax-xmin)+(next_j-xmin)] != 0.0 )
			  {
				double slope =  distance_current[(next_i-ymin)*(xmax-xmin)+(next_j-xmin)] - distance_current[(current_i-ymin)*(xmax-xmin)+(current_j-xmin)];
				slope = -1.0 *slope;
				if (direction == 1) {	 
				  point.x= current_j + distance_current[(current_i-ymin)*(xmax-xmin)+(current_j-xmin)]/slope;
				  point.y = current_i;
				}
				else if (direction == 3){	 
				  point.x = current_j - distance_current[(current_i-ymin)*(xmax-xmin)+(current_j-xmin)]/slope;
				  point.y = current_i;
				}
				else if (direction == 0) {	 				
				  point.y =  current_i - distance_current[(current_i-ymin)*(xmax-xmin)+(current_j-xmin)]/slope;
				  point.x = current_j;
				}
				else if (direction == 2){	
				  point.y =  current_i + distance_current[(current_i-ymin)*(xmax-xmin)+(current_j-xmin)]/slope;
				  point.x = current_j;
				}		 
			  }
			  else {
				cerr << "levelset on gridpoint  " << current_i << "\t" << current_j<<"\t" << (*domain)[current_i][current_j]<< endl;
				
				if (distance_current[(current_i-ymin)*(xmax-xmin)+(current_j-xmin)] == 0.0) { point.y =  current_i;  point.x = current_j; }
				else if (distance_current[(next_i-ymin)*(xmax-xmin)+(next_j-xmin)]  == 0.0) { point.y =  next_i;  point.x = next_j; }
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
    
	diff_xmin = diff_xmin - xmin; 
	diff_xmax = diff_xmax - xmax; 
	diff_ymin = diff_ymin - ymin; 
	diff_ymax = diff_ymax - ymax;
    
    
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
		if (xmin < 0) {cout <<"undefined box size for xmin: "<< xmin << endl; abort();}//xmin = 0;
		if (xmax > m) {cout <<"undefined box size for xmax: "<< xmax << endl; abort();}//xmax = m;
		if (ymin < 0) {cout <<"undefined box size for ymin: "<< ymin << endl; abort();}//ymin = 0;
		if (ymax > m) {cout <<"undefined box size for ymax: "<< xmax << endl; abort();}// ymax = m;
		
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



void LSbox::set_comparison(){
	int grid_blowup = (*handler).get_grid_blowup();
	int m = (*handler).get_ngridpoints();
	double h = handler->get_h();
	LSbox* zero = handler->zeroBox;
	for (int i = ymin; i < ymax; i++){
		for (int j = xmin; j < xmax; j++){
			if ((i <= grid_blowup) || (m-grid_blowup <= i) || (j <= grid_blowup) || (m-grid_blowup <= j)) {
				distance_new[(i-ymin)*(xmax-xmin)+(j-xmin)]  = -DELTA;
			}
			if( abs(distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]) < (0.9* DELTA) && (abs(distance_new[(i-ymin)*(xmax-xmin)+(j-xmin)]) < ( 0.9 * DELTA)) ) {
// 				 update only in a tube around the n boundary - numerical stability!s
				distance_new[(i-ymin)*(xmax-xmin)+(j-xmin)] = 0.5 * (distance_new[(i-ymin)*(xmax-xmin)+(j-xmin)] - distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)] );
			}
// 	perhaps better to shift this line to the convolution function:
	neighbors_old = neighbors;
// 	update the old neighborlist - for read access by the other grains in the next timestep
	delete [] distance_2neighbor;
}


bool LSbox::checkIntersect(LSbox* box2) {    
    if ((xmin > box2->xmax || xmax < box2->xmin || ymin > box2->ymax || ymax < box2->ymin)) return false;
    return true;
}

void LSbox::comparison(){
	int loop = owner->loop;
	std::vector<LSbox*>::iterator it_nn;

	for(it_nn = neighbors_2order.begin(); it_nn != neighbors_2order.end();){		
		
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
// 						after the Convolution the updated distancefunction is in the distance_new arry of each box. so we have to compare with this array. 
// 						the nearest value we save for comparison in the distance_current array of the current grain.
						if( abs((**it_nn).distance_new[(i-ymin)*(xmax-xmin)+(j-xmin)]) < 0.9*DELTA ){
// 							we are in the 0.9*DELTA-tube of grain (**it_nn) -> so we update our arrays and neighborlist
							if( distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)] < (**it_nn).distance_new[(i-ymin)*(xmax-xmin)+(j-xmin)] ){ 	
								if( IDLocal[(i-ymin)*(xmax-xmin)+(j-xmin)].empty() ) {
// 									falls noch kein nachbar vorhanden:
									distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]  = (**it.nn).distance_new[(i-ymin)*(xmax-xmin)+(j-xmin)];	
									IDLocal[(i-ymin)*(xmax-xmin)+(j-xmin)].push_back(*it_nn);									
								}
								else {
// 									new next neighbor found:
									distance_2neighbor[(i-ymin)*(xmax-xmin)+(j-xmin)] = distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)];
									distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]  = (**it_nn).distance_new[(i-ymin)*(xmax-xmin)+(j-xmin)];
									IDLocal[(i-ymin)*(xmax-xmin)+(j-xmin)].insert( IDLocal[(i-ymin)*(xmax-xmin)+(j-xmin)].begin(), *it_nn);										}
							}
							else if(  (**it_nn).distance_new[(i-ymin)*(xmax-xmin)+(j-xmin)] > distance_2neighbor[(i-ymin)*(xmax-xmin)+(j-xmin)] ){ 
// 								Kandidat ist n�her dran, als 2ter nachbar oder gleich
								distance_2neighbor[(i-ymin)*(xmax-xmin)+(j-xmin)] = (*(distance[i-ymin]))[j-xmin]; 
								IDLocal[(i-ymin)*(xmax-xmin)+(j-xmin)].insert( ++IDLocal[(i-ymin)*(xmax-xmin)+(j-xmin)].begin() , *it_nn);							  
							}
							else { 
								// probably there are more than 3 grains nearer than DELTA to this gridpoint
								IDLocal[(i-ymin)*(xmax-xmin)+(j-xmin)].push_back(*it_nn);						  
							}
						}
					}
				}
			}		
		}
		neighbors_2order.erase(it_nn);
		else it_nn++;
		
	}
	  // checke schnitt zum randkorn:
	checkIntersect_zero_grain();
// 	be careful for parralisation!!!!!
	set_comparison();
// 	write the compared values to the distance_new array
	
}


void LSbox::checkIntersect_zero_grain(){
	LSbox* boundary = handler->boundary;
	int grid_blowup = handler->get_grid_blowup();
	int m = handler->get_ngridpoints();
	if (!(xmin > boundary->xmin && xmax < boundary->xmax && ymin > boundary->ymin && ymax < boundary->ymax)){
// 	if (checkIntersect(mid_in[0])){
		for (int i = ymin; i < ymax; i++){
			for (int j = xmin; j < xmax; j++){	
				if ((i <= 2* grid_blowup) || (m-2*grid_blowup <= i) || (j <= 2*grid_blowup) || (m-2*grid_blowup <= j)){
					if(distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]  < (*boundary)[(i-ymin)*(xmax-xmin)+(j-xmin)]){ 
						distance_current[(i-ymin)*(xmax-xmin)+(j-xmin)]  = (*boundary)[(i-ymin)*(xmax-xmin)+(j-xmin)];						  
					}
				}
			}
		}
	}
}


void LSbox::add_n2o(){
	neighbors_2order = neighbors;
	neighbor.clear();
	vector<LSbox*> ::iterator it_com, it_n2o, it;
	bool just_in;						
	for(it = neighbors_old.begin(); it != neighbors_old.end(); it++){		
		for( it_n2o = (**it).neighbors_old.begin(); it_n2o != (**it).neighbors_old.end(); it_n2o++){
			just_in = false;
			for(it_com = neighbors_2order.begin(); it_com != neighbors_2order.end(); it_com++){
				if(*it_n2o == *it_com || this == *it_n2o) just_in = true;
			}
			if(just_in == false) neighbors_2order.push_back(*it_n2o);
		}
	}
}


void LSbox::redist_box() {
	int grid_blowup = owner->get_grid_blowup(); 
	double h = owner->get_h();
// 	plot_box(false);
	int m=ymax-ymin;
	int n=xmax-xmin;
	int ii, jj;
	double slope = 1;
	double candidate, i_slope, zero;
	
// 	resize the distance_current array. be careful because during this part of algorithm both arrays have not the same size!!
	
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
	copy_distances();
}



void LSbox::plot_box(bool distanceplot){
	cout <<" \nGrain  Info: " << endl;
	cout << " ID :" <<id << endl;
    cout << " xmin, xmax, ymin, ymax :" << xmin << " || "<< xmax << " || " << ymin << " || " << ymax << endl;
//     if (distanceplot==true) utils::print_2dim_array(distance,ymax-ymin,xmax-xmin);
// 		else cout << " no distance values in storage!" << endl;
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
	  
     if (distanceplot)
     {
       stringstream filename;
       filename<< "BoxDistance_"<< id << ".gnu";
       ofstream datei;
       datei.open(filename.str());
// 	for (int j = ymin; j < ymax; j++){
// 	    for (int i = xmin; i < xmax; i++){
// 		 datei << ::std::fixed << (*domain)[j][i] << "\t";
// 	      }
// 	    datei << endl;
// 	}
	
	for (int j = 0; j < handler->get_ngridpoints(); j++){
	    for (int i = 0; i < handler->get_ngridpoints(); i++){
		if( j >= ymin && j < ymax && i >=xmin && i < xmax) 
		datei << ::std::fixed << (*domain)[j][i] << "\t";
		else datei << ::std::fixed << -DELTA<< "\t";
	      }
	    datei << endl;
	}
	
	datei.close();
     }
}


double LSbox::mis_ori(LSbox* grain_2){
// 	here we could work direktly with quarternions
	return misorientationCubic(phi1,PHI,phi2,grain_2->get_phi1(), grain_2->get_PHI(), grain_2->get_phi2());
}






