#include "box.h"



LSbox::LSbox() :quaternion(NULL), inputDistance(NULL), outputDistance(NULL), local_weights(NULL){}

LSbox::LSbox(int id, int xmin, int xmax, int ymin, int ymax, grainhdl* owner) :
		id(id)
{

	handler = owner;
	quaternion = new double[4];
	inputDistance = new DimensionalBuffer<double>(xmin, ymin, xmax, ymax);
	outputDistance = new DimensionalBuffer<double>(xmin, ymin, xmax, ymax);
	
	IDLocal.resize((xmax - xmin) * (ymax - ymin));

	local_weights = new Weightmap(owner);

}

LSbox::LSbox(int aID, voro::voronoicell_neighbor& c, double *part_pos, grainhdl* owner) : id(aID), nvertices(0), handler(owner) {

	int grid_blowup = owner->get_grid_blowup(); 
	double h = owner->get_h();
    // determine size of grain
    quaternion = new double[4];
	
	(*(handler->mymath)).randomOriShoemakeQuat( quaternion );
	
    int xmax = 0; 
	int xmin = handler->get_ngridpoints(); 
	int ymax = 0; 
	int ymin = xmin;
	
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
	
	inputDistance = new DimensionalBuffer<double>(xmin, ymin, xmax, ymax);
	outputDistance = new DimensionalBuffer<double>(xmin, ymin, xmax, ymax);
	inputDistance->resizeToSquare(handler->get_ngridpoints());
	outputDistance->resizeToSquare(handler->get_ngridpoints());
	inputDistance->clearValues(0.0);
	outputDistance->clearValues(0.0);
	
	get_new_IDLocalSize();
	IDLocal.resize((xmaxId-xminId)*(ymaxId-yminId));
	
	local_weights=new Weightmap(owner);
	cout << "made a new box: xmin="<<xmin<< " xmax="<<xmax <<" ymin="<<ymin << " ymax="<<ymax<<endl;
}






LSbox::LSbox(int id, int nvertex, double* vertices, double phi1, double PHI, double phi2, grainhdl* owner) : id(id), nvertices(nvertex), handler(owner){
	quaternion = new double[4];

	double euler[3] = {phi1,PHI,phi2};
	(*(handler->mymath)).euler2quaternion( euler, quaternion );
	
	int grid_blowup = owner->get_grid_blowup(); 
	double h = owner->get_h();
    // determine size of grain
    int xmax = 0; 
	int xmin = handler->get_ngridpoints(); 
	int ymax = 0; 
	int ymin = xmin;
	    
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


	inputDistance = new DimensionalBuffer<double>(xmin, ymin, xmax, ymax);
	outputDistance = new DimensionalBuffer<double>(xmin, ymin, xmax, ymax);
	inputDistance->resizeToSquare(handler->get_ngridpoints());
	outputDistance->resizeToSquare(handler->get_ngridpoints());
	inputDistance->clearValues(0.0);
	outputDistance->clearValues(0.0);
	
	get_new_IDLocalSize();
	IDLocal.resize((xmaxId-xminId)*(ymaxId-yminId));	
	
	cout << "made a new box: xmin="<<xmin<< " xmax="<<xmax <<" ymin="<<ymin << " ymax="<<ymax<<endl;
	local_weights=new Weightmap(owner);
}

LSbox::~LSbox() {
	if(quaternion!=NULL){
		delete [] quaternion;
		delete inputDistance;
		delete outputDistance;
		delete local_weights;
	}
}

int LSbox::getID() {
    return id;
}

void LSbox::distancefunction(int nvertex, double* vertices){
// 	plot_box(false);
	int grid_blowup = handler->get_grid_blowup(); 
	double h = handler->get_h();
	int i,j,k;
	double d, dmin,lambda;
	vektor u(2), a(2), p(2), x1(2), x2(2);

	for (i=outputDistance->getMinY();i<outputDistance->getMaxY();i++){ // ¸ber gitter iterieren
	  for (j=outputDistance->getMinX();j<outputDistance->getMaxX();j++){
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
					d= abs(d);
					if(abs(d)< abs(dmin)) dmin=d;
				}
            }
			if (abs(dmin) < DELTA){
				outputDistance->setValueAt(i, j, dmin);
			}
			else{
				outputDistance->setValueAt(i, j, DELTA * utils::sgn(dmin));
			}
        }
	}
	
	int count = 0;
    for (j=outputDistance->getMinX();j<outputDistance->getMaxX();j++){ // ¸ber gitter iterieren
	    i=outputDistance->getMinY();
	    count = 0;
	    while( i<outputDistance->getMaxY() && count < 1) {
	    	outputDistance->setValueAt(i, j, -abs(outputDistance->getValueAt(i,j)));
		    if ( -outputDistance->getValueAt(i,j) <=  h )
		    	count++;
		    i++;
	    } 		
	    i=outputDistance->getMaxY()-1;
	    count =0;
	    while( i>=outputDistance->getMinY() && count < 1) {
	    	outputDistance->setValueAt(i, j, -abs(outputDistance->getValueAt(i,j)));
		    if ( -outputDistance->getValueAt(i,j) <= h )
		    	count++;
		    i--;
	    } 
    }

    for (i=outputDistance->getMinY();i<outputDistance->getMaxY();i++){
	    j=outputDistance->getMinX();
	    count = 0;
	    while( j<outputDistance->getMaxX() && count < 1 ) {
	    	outputDistance->setValueAt(i, j, -abs(outputDistance->getValueAt(i,j)));
		    if ( -outputDistance->getValueAt(i,j) <=  h )
		    	count++;
		    j++;
	    } 
	    j=outputDistance->getMaxX()-1;
	    count =0;
	    while( j>=outputDistance->getMinX() && count < 1  ) {
	    	outputDistance->setValueAt(i, j, -abs(outputDistance->getValueAt(i,j)));
		    if ( -outputDistance->getValueAt(i,j)<=  h )
		    	count++;
		    j--;
	    } 
    }
}

void LSbox::distancefunction(voro::voronoicell_neighbor& c, double *part_pos){
    int grid_blowup = handler->get_grid_blowup(); 
    double h = handler->get_h();
    
    int i,j,k;
    double d, dmin,lambda;

    vektor u(2), a(2), p(2), x1(2), x2(2);
    vector<double> vv;
    c.vertices (part_pos[3*(id-1)],part_pos[3*(id-1)+1],part_pos[3*(id-1)+2],vv);
    double domain_vertices[] = {0.,0.,1.,0.,1.,1.,0.,1.,0.,0.}; // array of vertices to loop over
    for (i=outputDistance->getMinY();i<outputDistance->getMaxY();i++){ // ¸ber gitter iterieren
    	  for (j=outputDistance->getMinX();j<outputDistance->getMaxX();j++){
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
						if(lambda <= 0) 					d = (p-x1).laenge();
						if((0 < lambda) && (lambda < 1)) 	d = (p-(a+(u*lambda))).laenge();
						if(lambda >= 1) 					d = (p-x2).laenge();
						if(abs(d)< abs(dmin)) 				dmin=d;
					}
				}
			}
			if (abs(dmin) < DELTA)
			{
				outputDistance->setValueAt(i, j, dmin);
			}
			else
			{
				outputDistance->setValueAt(i, j, DELTA * utils::sgn(dmin));
			}
		}
    }
	int count = 0;
    for (j=outputDistance->getMinX();j<outputDistance->getMaxX();j++){
	    i=outputDistance->getMinY();
	    count = 0;
	    while( i<outputDistance->getMaxY()  && count < 1) {
	    	outputDistance->setValueAt(i, j, -abs(outputDistance->getValueAt(i,j)));
		    if ( -outputDistance->getValueAt(i,j) <=  h )
		    	count++;
		    i++;
	    } 		
	    i=outputDistance->getMaxY()-1;
	    count =0;
	    while( i>=outputDistance->getMinY() && count < 1) {
	    	outputDistance->setValueAt(i, j, -abs(outputDistance->getValueAt(i,j)));
		    if ( -outputDistance->getValueAt(i,j) <= h )
		    	count++;
		    i--;
	    } 
    }

    for (i=outputDistance->getMinY();i<outputDistance->getMaxY();i++){
	    j=outputDistance->getMinX();
	    count = 0;
	    while( j<outputDistance->getMaxX()  && count < 1 ) {
	    	outputDistance->setValueAt(i, j, -abs(outputDistance->getValueAt(i,j)));
		    if ( -outputDistance->getValueAt(i,j) <=  h )
		    	count++;
		    j++;
	    } 
	    j=outputDistance->getMaxX()-1;
	    count =0;
	    while( j>=outputDistance->getMinX()   && count < 1  ) {
	    	outputDistance->setValueAt(i, j, -abs(outputDistance->getValueAt(i,j)));
		    if ( -outputDistance->getValueAt(i,j) <=  h )
		    	count++;
		    j--;
	    } 
    }
}


// Convolution und Helperfunctions 
/**************************************/
/**************************************/

void LSbox::convolution(){
	double h = handler->get_h();
	if(get_status() != true ) return;
	//  set references for the convolution step

	switch_in_and_out();
	double* ST = handler->ST;
	int n = outputDistance->getMaxX()-outputDistance->getMinX();
	int dt 	= handler->get_dt();
	
	fftw_complex *fftTemp;
	fftw_plan fwdPlan, bwdPlan;
	
	fftTemp = (fftw_complex*) fftw_malloc(n*(floor(n/2)+1)*sizeof(fftw_complex));
		
	makeFFTPlans(inputDistance->getRawData(),outputDistance->getRawData(), fftTemp, &fwdPlan, &bwdPlan);
	conv_generator(fftTemp,fwdPlan,bwdPlan);

	fftw_destroy_plan(fwdPlan);
	fftw_destroy_plan(bwdPlan);

	fftw_free (fftTemp);
	/*********************************************************************************/
	// Velocity Corrector Step: 
	/*********************************************************************************/
	// hier soll energycorrection gerechnet werden.
	// in der domainCl steht die urspr�nglich distanzfunktion, in dem arry die gefaltete
	
	if(!ISOTROPIC){	    
	    vector<LSbox*>::iterator it;
	    int intersec_xmin, intersec_xmax, intersec_ymin, intersec_ymax;
		double weight;
	    double val;
		double tubeRadius = sqrt(2)*h + 0.00001;
	    
	    if (xminId < outputDistance->getMinX())
		  intersec_xmin = outputDistance->getMinX();
	    else  intersec_xmin = xminId;
	    
	    if (yminId < outputDistance->getMinY())
		  intersec_ymin = outputDistance->getMinY();
	    else  intersec_ymin = yminId;
	    
	    if (xmaxId > outputDistance->getMaxX())
		  intersec_xmax= outputDistance->getMaxX();
	    else  intersec_xmax = xmaxId;
	    
	    if (ymaxId > outputDistance->getMaxY())
		  intersec_ymax= outputDistance->getMaxY();
	    else  intersec_ymax = ymaxId;
	    
	    for (int i = intersec_ymin; i < intersec_ymax; i++){
			for (int j = intersec_xmin; j < intersec_xmax; j++) {
				val = inputDistance->getValueAt(i,j);
				if(val <= tubeRadius ){
					if(IDLocal[(i-yminId)*(xmaxId-xminId) + (j-xminId)].size() != 2){ // || IDLocal[(i-yminId)*(xmaxId-xminId) + (j-xminId)].size() > 4 ) {
						weight=1;
						// for n < 3: here we are far away from a triple point - so we only take curvature into account.
						// for n >4: here we are at a multiple point, which only occurs because of unregular intialisation by voronoi cells - so we only take curvature into account.
					}
					else  weight = local_weights->loadWeights(IDLocal[(i-yminId)*(xmaxId-xminId) + (j-xminId)], this, handler->ST);
			
					outputDistance->setValueAt(i,j, val + (outputDistance->getValueAt(i,j) - val) * weight );
				}
			}
		}	   
	}
	
	IDLocal.clear();	
	get_new_IDLocalSize();
	IDLocal.resize((xmaxId-xminId)*(ymaxId-yminId));

// 	plot_box(true,1,"Convoluted_1");
// 	plot_box(true,2,"Convoluted_2");

	switch_in_and_out();

}
void LSbox::get_new_IDLocalSize(){
	xmaxId = outputDistance->getMaxX();
	xminId = outputDistance->getMinX();
	ymaxId = outputDistance->getMaxY();
	yminId = outputDistance->getMinY();
}



void LSbox::makeFFTPlans(double *in, double* out,fftw_complex *fftTemp, fftw_plan *fftplan1, fftw_plan *fftplan2)
{ /* creates plans for FFT and IFFT */
	int n = outputDistance->getMaxX() - outputDistance->getMinX();
	*fftplan1 = fftw_plan_dft_r2c_2d(n,n,in,fftTemp,FFTW_ESTIMATE);
	*fftplan2 = fftw_plan_dft_c2r_2d(n,n,fftTemp,out,FFTW_ESTIMATE);
}

void LSbox::conv_generator(fftw_complex *fftTemp, fftw_plan fftplan1, fftw_plan fftplan2)
{
	/* Function returns in u the updated value of u as described below..
	u -> (G_{dt})*u
	Assumptions:
	fftplan1 converts u to it's FT (in fftTemp), and
	fftplan2 converts the FT (after pointwise multiplication with G)
	(*outputDistance) to a real-valued level set function at u.
	Memory is already allocated in fftTemp
	(necessary to create the plans) */
	
	int n = outputDistance->getMaxX() - outputDistance->getMinX();
// 	int m = outputDistance->getMaxY() - outputDistance->getMinY();
// 	assert(m!=n);
	double dt = handler->get_dt();
	int n2 = floor(n/2) + 1;
	int nn = (*handler).get_ngridpoints();
	double nsq =  nn*nn; 
	double k = 2.0 * PI / n;
	double G;
	double coski;
	fftw_execute(fftplan1);
	for(int i=0;i<n2;i++) {
		coski=cos(k*i);
		for(int j=0;j<n;j++){
			// 	  G= exp((-2.0 * dt) * nsq * (2.0-cos(k*i)-cos(k*j)));			
			G = 2.0*(2.0 - coski - cos(k*j)) * nsq;
			G = 1.0/(1.0+(dt*G)) / (n*n);
			//        USE this line for Richardson-type extrapolation
			//       G = (4.0/pow(1+1.5*(dt)/40*G,40) - 1.0 / pow(1+3.0*(dt)/40*G,40)) / 3.0 / (double)(n*n);
			/* normalize G by n*n to pre-normalize convolution results */
			fftTemp[i+n2*j][0] = fftTemp[i+n2*j][0]*G;
			fftTemp[i+n2*j][1] = fftTemp[i+n2*j][1]*G;
		}
	}
	fftw_execute(fftplan2);
}

/**************************************/
/**************************************/


void LSbox::determineIDs(){
	DimensionalBuffer<double> distance_2neighbor(outputDistance->getMinX(), outputDistance->getMinY(),
										 	 	 outputDistance->getMaxX(), outputDistance->getMaxY());
	distance_2neighbor.clearValues(-1.0);
	inputDistance->clearValues(-1.0);
	int loop = 0;
	std::vector<LSbox*>::iterator it_nn;

	for(it_nn = neighbors.begin(); it_nn != neighbors.end(); it_nn++){		

			if (checkIntersect(*it_nn)){
				int x_min_new, x_max_new, y_min_new, y_max_new;
				
				if(outputDistance->getMinX() < (**it_nn).outputDistance->getMinX()) 
					  x_min_new = (**it_nn).outputDistance->getMinX();
				else x_min_new = outputDistance->getMinX();
				
				if(outputDistance->getMaxX() > (**it_nn).outputDistance->getMaxX()) 
					 x_max_new = (**it_nn).outputDistance->getMaxX();
				else x_max_new = outputDistance->getMaxX();
								
				if(outputDistance->getMinY() < (**it_nn).outputDistance->getMinY()) 
					y_min_new = (**it_nn).outputDistance->getMinY();
				else y_min_new = outputDistance->getMinY();
					
				if(outputDistance->getMaxY() > (**it_nn).outputDistance->getMaxY()) 
					 y_max_new = (**it_nn).outputDistance->getMaxY();
				else y_max_new = outputDistance->getMaxY();
					
// 				cout << "box: intersec_xmin="<<x_min_new<< " intersec_xmax="<<x_max_new <<" intersec_ymin="<<y_min_new << " intersec_ymax="<<y_max_new<<endl;
	
				for (int i = y_min_new; i < y_max_new; i++){
					for (int j = x_min_new; j < x_max_new; j++){					
						double dist = (**it_nn).outputDistance->getValueAt(i,j);
						if( dist > inputDistance->getValueAt(i,j) ){
							    if( !IDLocal[(i-yminId)*(xmaxId-xminId)+(j-xminId)].empty() ){ 
								      distance_2neighbor.setValueAt(i,j,inputDistance->getValueAt(i,j));
							    }
							    inputDistance->setValueAt(i, j, dist);
							    IDLocal[(i-yminId)*(xmaxId-xminId)+(j-xminId)].insert( IDLocal[(i-yminId)*(xmaxId-xminId)+(j-xminId)].begin(), *it_nn);	
						}
						else if(  dist > distance_2neighbor.getValueAt(i, j) ){ //candidate of neighbor is closer than 2nd neighbor
						    distance_2neighbor.setValueAt(i,j, dist);	
						    IDLocal[(i-yminId)*(xmaxId-xminId)+(j-xminId)].insert( ++IDLocal[(i-yminId)*(xmaxId-xminId)+(j-xminId)].begin() , *it_nn);							  
						}
						else { 
							IDLocal[(i-yminId)*(xmaxId-xminId)+(j-xminId)].push_back(*it_nn);						  
						}
					}
				}
			}
			
	}		

  
}


/**************************************/
/**************************************/

void LSbox::switch_in_and_out(){
	DimensionalBuffer<double>* temp;
	temp = inputDistance;
	inputDistance = outputDistance;
	outputDistance = temp;
}

// Comparison + Helperfunctions
/**************************************/
/**************************************/


void LSbox::set_comparison(){

	int grid_blowup = (*handler).get_grid_blowup();
	int m = (*handler).get_ngridpoints();
	double h = handler->get_h();
	LSbox* zero = handler->zeroBox;


	for (int i = outputDistance->getMinY(); i < outputDistance->getMaxY(); i++){
		for (int j = outputDistance->getMinX(); j < outputDistance->getMaxX(); j++){
			if( IDLocal[(i-yminId)*(xmaxId-xminId)+(j-xminId)].empty()) {
				outputDistance->setValueAt(i, j, inputDistance->getValueAt(i,j)); 
				continue;
			}
			else {
				if( abs(inputDistance->getValueAt(i,j)) < ( 0.7 * DELTA)){
					outputDistance->setValueAt(i, j, 0.5 * (inputDistance->getValueAt(i,j) - outputDistance->getValueAt(i,j)));
				}
				else if(inputDistance->getValueAt(i,j) > 0)
					outputDistance->setValueAt(i,j, DELTA);
				else if(inputDistance->getValueAt(i,j) < 0)
					outputDistance->setValueAt(i,j, -DELTA);
			}
		}
	}
	neighbors_old = neighbors;
}


bool LSbox::checkIntersect(LSbox* box2) {    
    if (inputDistance->getMinX() > box2->inputDistance->getMaxX() ||
    	inputDistance->getMaxX() < box2->inputDistance->getMinX() ||
    	inputDistance->getMinY() > box2->inputDistance->getMaxY() ||
    	inputDistance->getMaxY() < box2->inputDistance->getMinY())
    	return false;
    return true;
}


void LSbox::comparison(){
	if(get_status() != true ) return;
	
	DimensionalBuffer<double> distance_2neighbor(outputDistance->getMinX(), outputDistance->getMinY(),
										 	 	 outputDistance->getMaxX(), outputDistance->getMaxY());

	distance_2neighbor.clearValues(-1.0);
	outputDistance->clearValues(-1.0);
	//std::fill((*comparisonDistance).begin(),(*comparisonDistance).end(), -1.0);
	
	int loop = handler->loop;
	std::vector<LSbox*>::iterator it_nn;

	for(it_nn = neighbors_2order.begin(); it_nn != neighbors_2order.end();){		
			if( ((**it_nn).get_status() == true )) {
			if (checkIntersect(*it_nn)){
				neighbors.push_back(*it_nn);
				int x_min_new, x_max_new, y_min_new, y_max_new;
				
				if(inputDistance->getMinX() < (**it_nn).inputDistance->getMinX()) x_min_new = (**it_nn).inputDistance->getMinX();
					else x_min_new = inputDistance->getMinX();
				
				if(inputDistance->getMaxX() > (**it_nn).inputDistance->getMaxX()) x_max_new = (**it_nn).inputDistance->getMaxX();
					else x_max_new = inputDistance->getMaxX();
								
				if(inputDistance->getMinY() < (**it_nn).inputDistance->getMinY()) y_min_new = (**it_nn).inputDistance->getMinY();
					else y_min_new = inputDistance->getMinY();
					
				if(inputDistance->getMaxY() > (**it_nn).inputDistance->getMaxY()) y_max_new = (**it_nn).inputDistance->getMaxY();
					else y_max_new = inputDistance->getMaxY();
					
// 				cout << "box: intersec_xmin="<<x_min_new<< " intersec_xmax="<<x_max_new <<" intersec_ymin="<<y_min_new << " intersec_ymax="<<y_max_new<<endl;
	
				for (int i = y_min_new; i < y_max_new; i++){
					for (int j = x_min_new; j < x_max_new; j++){					
// 						after the Convolution the updated distancefunction is in the distanceBuffer2 array of each box. so we have to compare with this array. 
// 						the nearest value we save for comparison in the distanceBuffer2 array of the current grain.						
						if(abs(inputDistance->getValueAt(i,j)) < (0.7*DELTA)){
							double dist = (**it_nn).getDistance(i,j);
							if( abs(dist) < (0.7*DELTA)){								
								if( dist > outputDistance->getValueAt(i,j) ){
										if( !IDLocal[(i-yminId)*(xmaxId-xminId)+(j-xminId)].empty() ){ 
											distance_2neighbor.setValueAt(i,j,outputDistance->getValueAt(i,j));
										//Question: OutputDistance??? not Input?!
										}
										outputDistance->setValueAt(i, j, dist);
										IDLocal[(i-yminId)*(xmaxId-xminId)+(j-xminId)].insert( IDLocal[(i-yminId)*(xmaxId-xminId)+(j-xminId)].begin(), *it_nn);	
	// 								}
								}
								else if(  dist > distance_2neighbor.getValueAt(i, j) ){ //candidate of neighbor is closer than 2nd neighbor
									distance_2neighbor.setValueAt(i,j, dist);
									IDLocal[(i-yminId)*(xmaxId-xminId)+(j-xminId)].insert( ++IDLocal[(i-yminId)*(xmaxId-xminId)+(j-xminId)].begin() , *it_nn);							  
								}
								else { 
									IDLocal[(i-yminId)*(xmaxId-xminId)+(j-xminId)].push_back(*it_nn);						  
								}
							}
						}
					}
				}
			}
			
		}
		neighbors_2order.erase(it_nn);
	}

if(loop>29 && id == 48)
{
	plot_box(true,1,"Compare");
	plot_box(true,2,"Compare");
}	
// 	char buffer;
// 
// 	  // checke schnitt zum randkorn:
	checkIntersect_zero_grain();
if(loop>29 &&id == 48)
{
  
	plot_box(true,1,"Compare_zero");
	plot_box(true,2,"Compare_zero");
}	

	// 	be careful for parralisation!!!!!

	set_comparison();
if(loop>29 &&id == 48){	
	plot_box(true,1,"Compare_set");
	plot_box(true,2,"Compare_set");
}

}



double LSbox::getDistance(int i, int j){
   return inputDistance->getValueAt(i,j);
}




void LSbox::checkIntersect_zero_grain(){
	LSbox* boundary = handler->boundary;
	int grid_blowup = handler->get_grid_blowup();
	double h = handler->get_h();
	int m = handler->get_ngridpoints();
	if (!(outputDistance->getMinX() > boundary->outputDistance->getMinX() &&
		  outputDistance->getMaxX() < boundary->outputDistance->getMaxX() &&
		  outputDistance->getMinY() > boundary->outputDistance->getMinY() &&
		  outputDistance->getMaxY() < boundary->outputDistance->getMaxY())){
// 	if (checkIntersect(mid_in[0])){
		for (int i = outputDistance->getMinY(); i < outputDistance->getMaxY(); i++){
			for (int j = outputDistance->getMinX(); j < outputDistance->getMaxX(); j++){
				if ((i <= 2*grid_blowup) || (m-2*grid_blowup <= i) || (j <= 2*grid_blowup) || (m-2*grid_blowup <= j)){
					if(outputDistance->getValueAt(i,j) < boundary->outputDistance->getValueAt(i,j)){
						outputDistance->setValueAt(i,j, boundary->outputDistance->getValueAt(i,j));
						IDLocal[(i-yminId)*(xmaxId-xminId)+(j-xminId)].push_back(boundary);	
					}
				}
			}
		}
	}
}


void LSbox::add_n2o(){
	if(get_status() != true ) return;
	neighbors_2order = neighbors;
	neighbors.clear();
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

/**************************************/
// end of Comparison
/**************************************/



// Find Contour operates on inputDistance
/**************************************/
/**************************************/

void LSbox::find_contour() {
	if(get_status() != true ) return;
	switch_in_and_out();
	exist = false;
	
	// save old boundaries -> function will compute updates

    contourGrain.clear();

    MarchingSquaresAlgorithm marcher(*inputDistance);
    exist = marcher.generateContour(contourGrain);
	if(!exist) return;
    int grid_blowup = handler->get_grid_blowup();

    int xminNew = 999999, xmaxNew = -1, yminNew = 999999, ymaxNew = -1;
    for(int i=0; i<contourGrain.size(); i++)
    {
    	if(int(contourGrain[i].x + 0.5) - grid_blowup < xminNew)
    		xminNew = int(contourGrain[i].x + 0.5) - grid_blowup;
    	if(int(contourGrain[i].x + 0.5) + grid_blowup > xmaxNew)
    		xmaxNew = int(contourGrain[i].x + 0.5) + grid_blowup;
    	if(int(contourGrain[i].y + 0.5) - grid_blowup < yminNew)
    		yminNew = int(contourGrain[i].y + 0.5) - grid_blowup;
    	if(int(contourGrain[i].y + 0.5) + grid_blowup > ymaxNew)
    		ymaxNew = int(contourGrain[i].y + 0.5) + grid_blowup;
    }

    outputDistance->resize(xminNew, yminNew, xmaxNew, ymaxNew);
	outputDistance->resizeToSquare(handler->get_ngridpoints());
    
	double h = handler->get_h();
	int loop = handler->loop;
    int m = handler->get_ngridpoints();
    stringstream s;
    // compute Volume and Energy
    if ( (loop % int(ANALYSESTEP)) == 0 || loop == TIMESTEPS ) {
		vector<SPoint>::iterator volumeit=contourGrain.begin();
		energy = 0.0;
		volume = 0.0;
		double px, py;
		double gamma_hagb = 0.6;
		double theta_ref = 15.0* PI / 180.;
		double theta_mis;
		
		px= (*volumeit).x;
		py= (*volumeit).y;

		for (; volumeit!= contourGrain.end(); volumeit++){
			s << (*volumeit).x << "\t" << (*volumeit).y<<endl;
			volume += (py+(*volumeit).y)*(px-(*volumeit).x);
			px= (*volumeit).x;
			py= (*volumeit).y;
// 			cout << py << "  " << px << endl;
// 			cout << "assoziated grid point " << int(py+0.5) << "  " << int(px+0.5) << endl;
// 			if(ISOTROPIC){
// 				l= sqrt( )
// 			}
// 			else{
				if (!IDLocal[ ( (int(py+0.5)-yminId) * (xmaxId-xminId)) + (int(px+0.5) - xminId)].empty()){
					theta_mis=mis_ori( IDLocal[ ( (int(py+0.5)-yminId) * (xmaxId-xminId)) + (int(px+0.5) - xminId) ][0] ); //find the LSbox pointer to the next neighbor -> therefor find the next grid pointer
					if (theta_mis <= theta_ref)	energy += h* gamma_hagb * ( theta_mis / theta_ref) * (1.0 - log( theta_mis / theta_ref));
						else energy += h* gamma_hagb;
				}
// 			}
		}				
		stringstream dateiname;
		dateiname << "Contourline_" << id << ".gnu";
		ofstream datei;
		datei.open(dateiname.str());
		datei << s.str();
		datei.close();
		volume = abs(volume)*0.5;
		cerr<< "Volume of " << id << "= " << volume << endl;
		cerr<< "Surface Energy of " << id << "= " << abs(energy)*0.5<< endl << endl;

	}
	return;
}


/**************************************/
//  Redistancing
/**************************************/

void LSbox::redist_box() {	
	if(get_status() != true ) return;
	int grid_blowup = handler->get_grid_blowup(); 
	double h = handler->get_h();
	double slope = 1;
	double candidate, i_slope, distToZero;
	
	outputDistance->clearValues(-1.0);
	
// 	resize the outputDistance array. be careful because during this part of algorithm both arrays have not the same size!!
	int intersec_xmin, intersec_xmax, intersec_ymin, intersec_ymax;
	    
	if (inputDistance->getMinX() < outputDistance->getMinX())
		intersec_xmin = outputDistance->getMinX();
	else  intersec_xmin = inputDistance->getMinX();
	
	if (inputDistance->getMinY() < outputDistance->getMinY())
		intersec_ymin = outputDistance->getMinY();
	else  intersec_ymin = inputDistance->getMinY();
	
	if (inputDistance->getMaxX() < outputDistance->getMaxX())
		intersec_xmax= inputDistance->getMaxX();
	else  intersec_xmax = outputDistance->getMaxX();
	
	if (inputDistance->getMaxY() < outputDistance->getMaxY())
		intersec_ymax= inputDistance->getMaxY();
	else  intersec_ymax = outputDistance->getMaxY();
	
// 	cout << "box: intersec_xmin="<<intersec_xmin<< " intersec_xmax="<<intersec_xmax <<" intersec_ymin="<<intersec_ymin << " intersec_ymax="<<intersec_ymax<<endl;

	for (int i = intersec_ymin; i < outputDistance->getMaxY(); i++){
	  for (int j = intersec_xmin; j < outputDistance->getMaxX()-1; j++) {
			// x-direction forward
			if(j < intersec_xmax-1 && i < intersec_ymax ){
				if (inputDistance->getValueAt(i,j) * inputDistance->getValueAt(i,j+1) <= 0.0) {
					// interpolate
					i_slope  = ( inputDistance->getValueAt(i,j+1) - inputDistance->getValueAt(i,j) ) / h;
					distToZero = - inputDistance->getValueAt(i,j) / i_slope;
					if ( abs(outputDistance->getValueAt(i,j) ) > abs(distToZero))
						outputDistance->setValueAt(i,j,-distToZero * utils::sgn(i_slope));
				}
				candidate = outputDistance->getValueAt(i,j) + (utils::sgn(inputDistance->getValueAt(i,j+1)) * h);
				if (abs(candidate) < abs(outputDistance->getValueAt(i,j+1)))
					outputDistance->setValueAt(i,j+1, candidate);
			}
			else {
				candidate = outputDistance->getValueAt(i,j)  + (utils::sgn( outputDistance->getValueAt(i,j+1) ) * h);
				if (abs(candidate) < abs(outputDistance->getValueAt(i, j+1)))
					outputDistance->setValueAt(i,j+1, candidate);
			}
		}		
	}
	
	for (int i = intersec_ymin; i < outputDistance->getMaxY(); i++){
	  for (int j = intersec_xmax-1; j >  outputDistance->getMinX(); j--) {
	// x-direction outputDistanceward
			//check for sign change
			if(j > intersec_xmin && i < intersec_ymax){
					// calculate new distance candidate and assign if appropriate
					candidate = outputDistance->getValueAt(i,j)  + (utils::sgn( inputDistance->getValueAt(i,j-1) ) * h);
					if (abs(candidate) < abs(outputDistance->getValueAt(i,j-1)))
						outputDistance->setValueAt(i,j-1, candidate);
			}
			else {
				candidate = outputDistance->getValueAt(i,j)  + utils::sgn( outputDistance->getValueAt(i,j-1))*h;
				if (abs(candidate) < abs(outputDistance->getValueAt(i,j-1)))
					outputDistance->setValueAt(i,j-1, candidate);
			}
		}		
	}
		
	// y-direction forward
	for (int j = intersec_xmin; j < outputDistance->getMaxX(); j++) {
		for (int i = intersec_ymin; i < outputDistance->getMaxY()-1; i++) {
			if(j < intersec_xmax && i < intersec_ymax-1){
				if (inputDistance->getValueAt(i,j) * inputDistance->getValueAt(i+1,j) <= 0.0) {
					// interpolate
					i_slope  = (inputDistance->getValueAt(i+1,j) - inputDistance->getValueAt(i,j) )/ h;
					distToZero = - inputDistance->getValueAt(i,j) / i_slope;
					if ( abs(outputDistance->getValueAt(i,j) ) > abs(distToZero))
						outputDistance->setValueAt(i,j, -distToZero * utils::sgn(i_slope));
				}
				// calculate new distance candidate and assign if appropriate
				candidate = outputDistance->getValueAt(i,j)  + (utils::sgn( inputDistance->getValueAt(i+1,j) ) * h);
				if (abs(candidate) < abs(outputDistance->getValueAt(i+1,j) ))
					outputDistance->setValueAt(i+1,j, candidate);
			}
			else {
				candidate = outputDistance->getValueAt(i,j)  + (utils::sgn( outputDistance->getValueAt(i+1,j)) * h);
				if (abs(candidate) < abs(outputDistance->getValueAt(i+1,j)))
					outputDistance->setValueAt(i+1,j, candidate);
			}
		}		
	}

	for (int j = intersec_xmin; j < outputDistance->getMaxX(); j++) {
		for (int i = intersec_ymax-1; i > outputDistance->getMinY(); i--) {
			if(j < intersec_xmax && i > intersec_ymin){
				// calculate new distance candidate and assign if appropriate
				candidate = outputDistance->getValueAt(i,j)  + (utils::sgn( inputDistance->getValueAt(i-1, j) ) * h);
				if (abs(candidate) < abs(outputDistance->getValueAt(i-1, j) ))
					outputDistance->setValueAt(i-1, j, candidate);
			}
			else {
				candidate = outputDistance->getValueAt(i,j)  /+ (utils::sgn( outputDistance->getValueAt(i-1, j) ) * h);
				if (abs(candidate) < abs(outputDistance->getValueAt(i-1,j)))
					outputDistance->setValueAt(i-1, j, candidate);
			}
		}		
	}
	
	outputDistance->clampValues(-DELTA, DELTA);
	
// 	plot_box(true,1,"Redist_1");
// 	plot_box(true,2,"Redist_2");
	
	inputDistance->resize(outputDistance->getMinX(), outputDistance->getMinY(), outputDistance->getMaxX(), outputDistance->getMaxY());	
	// 	 set the references for the convolution step
}

/**************************************/
// end of redist
/**************************************/



/**************************************/
// plot the box and all its properties
/**************************************/

void LSbox::plot_box_contour(int loop, ofstream *dateiname, bool plotEnergyFunctional)
{
    vector<SPoint>::iterator contourIterator;
    if ( plotEnergyFunctional ){

	for (contourIterator= contourGrain.begin(); contourIterator != contourGrain.end(); contourIterator++){
	    *dateiname << (*contourIterator).x << "\t" << (*contourIterator).y<< "\t" << id << endl;
										//TODO change id to energy
	  
	}
    }
    else {
	for (contourIterator= contourGrain.begin(); contourIterator != contourGrain.end(); contourIterator++){
	    *dateiname << (*contourIterator).x << "\t" << (*contourIterator).y << endl;
	}
    }
    *dateiname << endl;
}


void LSbox::plot_box(bool distanceplot, int select, string simstep){
  
	cout <<" \nGrain  Info: " << endl;
	cout << " ID :" <<id << endl;
	cout << " xminIn, xmaxIn, yminIn, ymaxIn :" << inputDistance->getMinX() << " || "<< inputDistance->getMaxX() << " || " << inputDistance->getMinY() << " || " << inputDistance->getMaxY() << endl;
	cout << " xminOut, xmaxOut, yminOut, ymaxOut :" << outputDistance->getMinX() << " || "<< outputDistance->getMaxX() << " || " << outputDistance->getMinY() << " || " << outputDistance->getMaxY() << endl;
	cout << " xminId, xmaxId, yminId, ymaxId :" << xminId << " || "<< xmaxId << " || " << yminId << " || " << ymaxId << endl;
//     if (distanceplot==true) utils::print_2dim_array(distance,ymax-ymin,xmax-xmin);
// 		else cout << " no distance values in storage!" << endl;
	cout << " quaternion: " << quaternion[0] << " || "<< quaternion[1]<< " || " <<quaternion[2] << " || " << quaternion[3]<< endl;
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
		ofstream datei;
		int loop = handler->loop;
		if(select == 2) {
		filename<< "BoxDistance_"<< simstep << "out_T" << loop << "_" << id << ".gnu";
		datei.open(filename.str());
			for (int i = 0; i < handler->get_ngridpoints(); i++){
				for (int j = 0; j < handler->get_ngridpoints(); j++){
					if( i >= outputDistance->getMinY() && i < outputDistance->getMaxY() && j >=outputDistance->getMinX() && j < outputDistance->getMaxX()) {
						datei << ::std::fixed << outputDistance->getValueAt(i,j) << "\t";
					}
					else datei << ::std::fixed << -DELTA<< "\t";
				}
			datei << endl;
			}	
		}		
		if(select == 1) {
		filename<< "BoxDistance_"<< simstep << "in_T" << loop << "_" << id << ".gnu";
		datei.open(filename.str());
			for (int i = 0; i < handler->get_ngridpoints(); i++){
				for (int j = 0; j < handler->get_ngridpoints(); j++){
					if( i >= inputDistance->getMinY() && i < inputDistance->getMaxY() && j >=inputDistance->getMinX() && j < inputDistance->getMaxX()) {
						datei << ::std::fixed << inputDistance->getValueAt(i,j)<< "\t";
					}
					else datei << ::std::fixed << -DELTA<< "\t";
				}
			datei << endl;
			}	
		}		
		datei.close();
    }
}


double LSbox::mis_ori(LSbox* grain_2){
	if(get_status() != true ) {cout << "try to compute misori for are disappeared grains"; char buf; cin >> buf;}
// 	here we could work direktly with quarternions
	return (*(handler->mymath)).misorientationCubicQxQ( quaternion[0], quaternion[1], quaternion[2], quaternion[3], grain_2->quaternion[0], grain_2->quaternion[1], grain_2->quaternion[2], grain_2->quaternion[3] );
}



void LSbox::shape_distance(){
	//TODO: WORK ON THIS
	int m = outputDistance->getMaxX()-outputDistance->getMinX();
	for (int i = 0; i < outputDistance->getMaxY()-outputDistance->getMinY(); i++) {
		for (int j = 0; j < m; j++) {
			outputDistance->getRawData()[i*m +j] *= -1.0;
		}
	}
}


