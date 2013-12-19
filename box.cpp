#include "box.h"



LSbox::LSbox() {}

LSbox::LSbox(int id, int xmin, int xmax, int ymin, int ymax, double phi1,
		double PHI, double phi2, grainhdl* owner) :
		id(id), phi1(phi1),
		PHI(PHI), phi2(phi2)
{

	handler = owner;
	
	inputDistance = new DimensionalBuffer<double>(xmin, ymin, xmax, ymax);
	outputDistance = new DimensionalBuffer<double>(xmin, ymin, xmax, ymax);
	
	IDLocal.resize((xmax - xmin) * (ymax - ymin));

	local_weights = new Weightmap(owner);

}

LSbox::LSbox(int aID, voro::voronoicell_neighbor& c, double *part_pos, grainhdl* owner) : id(aID), phi1(0), PHI(0), phi2(0), nvertices(0), handler(owner) {

	int grid_blowup = owner->get_grid_blowup(); 
	double h = owner->get_h();
    // determine size of grain
    
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






LSbox::LSbox(int id, int nvertex, double* vertices, double phi1, double PHI, double phi2, grainhdl* owner) : id(id), phi1(phi1), PHI(PHI), phi2(phi2), nvertices(nvertex), handler(owner){

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
	//TODO: Fix this
	IDLocal.resize((xmaxId-xminId)*(ymaxId-yminId));	
	
	cout << "made a new box: xmin="<<xmin<< " xmax="<<xmax <<" ymin="<<ymin << " ymax="<<ymax<<endl;
	local_weights=new Weightmap(owner);
}

LSbox::~LSbox() {
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
		// 					if(((grid_blowup < i) && (i < (m- grid_blowup))) && ((grid_blowup < j) && (j < (m- grid_blowup)))) {
		// 						d=abs(d);
		// 					}
		// 					else d= abs(d);
					d= abs(d);
					if(abs(d)< abs(dmin)) dmin=d;
				}
            }
			// 			(*domain)[i][j]= dmin;
			if (abs(dmin) < DELTA)
			{
				//(*outputDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)]= dmin;
				outputDistance->setValueAt(i, j, dmin);
			}
			else
			{
				//(*outputDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)]= DELTA * utils::sgn(dmin);
				outputDistance->setValueAt(i, j, DELTA * utils::sgn(dmin));
			}
        }
	}
	
	int count = 0;
    for (j=outputDistance->getMinX();j<outputDistance->getMaxX();j++){ // ¸ber gitter iterieren
	    i=outputDistance->getMinY();
	    count = 0;
	    while( i<outputDistance->getMaxY() && count < 1) {
		    //(*outputDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)] = - abs((*outputDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)]);
	    	outputDistance->setValueAt(i, j, -abs(outputDistance->getValueAt(i,j)));
		    if ( -outputDistance->getValueAt(i,j) <=  h )
		    	count++;
		    i++;
	    } 		
	    i=outputDistance->getMaxY()-1;
	    count =0;
	    while( i>=outputDistance->getMinY() && count < 1) {
		    //(*outputDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)] = - abs((*outputDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)]);
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
		    //(*outputDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)] = - abs((*outputDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)]);
	    	outputDistance->setValueAt(i, j, -abs(outputDistance->getValueAt(i,j)));
		    if ( -outputDistance->getValueAt(i,j) <=  h )
		    	count++;
		    j++;
	    } 
	    j=outputDistance->getMaxX()-1;
	    count =0;
	    while( j>=outputDistance->getMinX() && count < 1  ) {
		    //(*outputDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)] = - abs((*outputDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)]);
	    	outputDistance->setValueAt(i, j, -abs(outputDistance->getValueAt(i,j)));
		    if ( -outputDistance->getValueAt(i,j)<=  h )
		    	count++;
		    j--;
	    } 
    }
    plot_box(true,1,"Dist_1");
	plot_box(true,2,"Dist_2");
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
    	  for (j=outputDistance->getMinX();j<outputDistance->getMaxY();j++){
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
		// 			(*domain)[i][j]= dmin;
			if (abs(dmin) < DELTA)
			{
				//(*outputDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)]= dmin;
				outputDistance->setValueAt(i, j, dmin);
			}
			else
			{
				//(*outputDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)]= DELTA * utils::sgn(dmin);
				outputDistance->setValueAt(i, j, DELTA * utils::sgn(dmin));
			}
		}
    }
	int count = 0;
    for (j=outputDistance->getMinX();j<outputDistance->getMaxX();j++){
	    i=outputDistance->getMinY();
	    count = 0;
	    while( i<outputDistance->getMaxY()  && count < 1) {
		    //(*outputDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)] = - abs((*outputDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)]);
	    	outputDistance->setValueAt(i, j, -abs(outputDistance->getValueAt(i,j)));
		    if ( -outputDistance->getValueAt(i,j) <=  h )
		    	count++;
		    i++;
	    } 		
	    i=outputDistance->getMaxY()-1;
	    count =0;
	    while( i>=outputDistance->getMinY() && count < 1) {
		    //(*outputDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)] = - abs((*outputDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)]);
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
		    //(*outputDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)] = - abs((*outputDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)]);
	    	outputDistance->setValueAt(i, j, -abs(outputDistance->getValueAt(i,j)));
		    if ( -outputDistance->getValueAt(i,j) <=  h )
		    	count++;
		    j++;
	    } 
	    j=outputDistance->getMaxX()-1;
	    count =0;
	    while( j>=outputDistance->getMinX()   && count < 1  ) {
		    //(*outputDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)] = - abs((*outputDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)]);
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
	//  set references for the convolution step

	switch_in_and_out();
	plot_box(true,1,"DistSwitch_1");
	plot_box(true,2,"DistSwitch_2");
	double* ST = handler->ST;
	int n = outputDistance->getMaxX()-outputDistance->getMinX();
	int m = outputDistance->getMaxY()-outputDistance->getMinY();
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
// 	    double rad =  DELTA* 0.7; // radius in dem ein drag wirkt
	    double weight;
// 	    int* rep = new int[3];
	    vector<LSbox*>::iterator it;
	    int intersec_xmin, intersec_xmax, intersec_ymin, intersec_ymax;
	    
	    if (xminId < outputDistance->getMinX())
		  intersec_xmin = outputDistance->getMinX();
	    else  intersec_xmin = xminId;
	    
	    if (yminId < outputDistance->getMinY())
		  intersec_ymin = outputDistance->getMinY();
	    else  intersec_ymin = yminId;
	    
	    if (xmaxId > outputDistance->getMaxX())
		  intersec_xmax= xmaxId;
	    else  intersec_xmax = outputDistance->getMaxX();
	    
	    if (ymaxId > outputDistance->getMaxY())
		  intersec_ymax= ymaxId;
	    else  intersec_ymax = outputDistance->getMaxY();
	    
	    
	    for (int i = intersec_ymin; i < intersec_ymax; i++){
		  for (int j = intersec_xmin; j < intersec_xmax; j++) {
    
		    // 		    if ( rad < abs(ref[i][j]) ) continue;
			weight = local_weights->loadWeights(IDLocal[(i-yminId)*(xmaxId-xminId) + (j-xminId)], this, handler->ST);
		    // 	      	    weight = ( 1-abs(rad - abs(ref[i][j])) ) * weight;		nur sinnvoll um einen drag zu simulieren	
			//(*outputDistance)[(i-yminOut)*(xmaxOut-xminOut)+j-xminOut] = (*inputDistance)[(i-yminIn)*(xmaxIn-xminIn)+j-xminIn] + (((*outputDistance)[(i-yminOut)*(xmaxOut-xminOut)+j-xminOut] -(*inputDistance)[(i-yminIn)*(xmaxIn-xminIn)+j-xminIn]) * weight);
			outputDistance->setValueAt(i,j,
					inputDistance->getValueAt(i,j) + (outputDistance->getValueAt(i,j) - inputDistance->getValueAt(i,j))*weight );
		  }
	   }
	   
	}
	
	IDLocal.clear();	
	get_new_IDLocalSize();
	IDLocal.resize((xmaxId-xminId)*(ymaxId-yminId));

// 	plot_box(true,1,"Convoluted_1");
// 	plot_box(true,2,"Convoluted_2");
// 	

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
	int m = outputDistance->getMaxY() - outputDistance->getMinY();
	*fftplan1 = fftw_plan_dft_r2c_2d(n,m,in,fftTemp,FFTW_ESTIMATE);
	*fftplan2 = fftw_plan_dft_c2r_2d(n,m,fftTemp,out,FFTW_ESTIMATE);
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
	int m = outputDistance->getMaxY() - outputDistance->getMinY();
// 	assert(m!=n);
	double dt = handler->get_dt();
	int n2 = floor(n/2) + 1;
	
	double nsq =  n*n; 
	double k = 2.0 * PI / n;
	double G;
	double coski;
	fftw_execute(fftplan1);
	for(int i=0;i<n2;i++) {
		coski=cos(k*i);
		for(int j=0;j<n;j++){
			// 	  G= exp((-2.0 * dt) * nsq * (2.0-cos(k*i)-cos(k*j)));			
			G = 2.0*(2.0 - coski - cos(k*j)) * nsq;
			G = 1.0/(1.0+(dt*G)) / nsq;
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
			//if( abs((*inputDistance)[(i-yminIn)*(xmaxIn-xminIn)+(j-xminIn)]) < ( 0.7 * DELTA) ) {7
			if( abs(inputDistance->getValueAt(i,j)) < ( 0.7 * DELTA) ) {
// 				 update only in a tube around the n boundary - numerical stability!s
// 				outputDistance and comparisonDistance point to the same object!
				//(*outputDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)] = 0.5 * ((*inputDistance)[(i-yminIn)*(xmaxIn-xminIn)+(j-xminIn)] - (*comparisonDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)] );
				outputDistance->setValueAt(i, j, 0.5 * (inputDistance->getValueAt(i,j) - outputDistance->getValueAt(i,j)));
			}
			//else if((*inputDistance)[(i-yminIn)*(xmaxIn-xminIn)+(j-xminIn)] > 0)
			else if(inputDistance->getValueAt(i,j) > 0)
				outputDistance->setValueAt(i,j, DELTA);
			//else if((*inputDistance)[(i-yminIn)*(xmaxIn-xminIn)+(j-xminIn)] < 0)
			else if(inputDistance->getValueAt(i,j) < 0)
				outputDistance->setValueAt(i,j, -3*DELTA);
			if ((i <= grid_blowup) || (m-grid_blowup <= i) || (j <= grid_blowup) || (m-grid_blowup <= j)) {
				outputDistance->setValueAt(i,j, -DELTA);
			}
// 			else outputDistance[(i-ymin)*(xmax-xmin)+(j-xmin)] = (*inputDistance)[(i-yminIn)*(xmaxIn-xminIn)+(j-xminIn)];
		}
	}
	
// 	perhaps better to shift this line to the convolution function:
	neighbors_old = neighbors;
// 	update the old neighborlist - for read access by the other grains in the next timestep
// 	delete [] distance_2neighbor;
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
	
// 	switch_in_and_out();
	//vector<double>* comparisonDistance = outputDistance;
	
	DimensionalBuffer<double> distance_2neighbor(outputDistance->getMinX(), outputDistance->getMinY(),
										 	 	 outputDistance->getMaxX(), outputDistance->getMaxX());

	distance_2neighbor.clearValues(-1.0);
	outputDistance->clearValues(-1.0);
	//std::fill((*comparisonDistance).begin(),(*comparisonDistance).end(), -1.0);
	
	int loop = handler->loop;
	std::vector<LSbox*>::iterator it_nn;

	for(it_nn = neighbors_2order.begin(); it_nn != neighbors_2order.end();){		
			cout << endl << "starting comparison between "<< id << "  and " << (*it_nn)->get_id()<< endl;
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
					
				cout << "box: intersec_xmin="<<x_min_new<< " intersec_xmax="<<x_max_new <<" intersec_ymin="<<y_min_new << " intersec_ymax="<<y_max_new<<endl;
	
				for (int i = y_min_new; i < y_max_new; i++){
					for (int j = x_min_new; j < x_max_new; j++){					
// 						after the Convolution the updated distancefunction is in the distanceBuffer2 array of each box. so we have to compare with this array. 
// 						the nearest value we save for comparison in the distanceBuffer2 array of the current grain.
						double dist = (**it_nn).getDistance(i,j);
						//if(abs((*inputDistance)[(i-yminIn)*(xmaxIn-xminIn)+(j-xminIn)])  < (0.7*DELTA) &&  abs(dist) < (0.7*DELTA)){
						if(abs(inputDistance->getValueAt(i,j)) < (0.7*DELTA) &&  abs(dist) < (0.7*DELTA)){
							// check if we are currently in the delta-area around the contour
 							//if( dist > (*comparisonDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)] ){
							if( dist > outputDistance->getValueAt(i,j) ){
									if( !IDLocal[(i-yminId)*(xmaxId-xminId)+(j-xminId)].empty() ){ 
										//distance_2neighbor[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)] = (*comparisonDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)];
										distance_2neighbor.setValueAt(i,j,outputDistance->getValueAt(i,j));
									}
									//(*comparisonDistance)[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)]  = dist;
									outputDistance->setValueAt(i, j, dist);
									IDLocal[(i-yminId)*(xmaxId-xminId)+(j-xminId)].insert( IDLocal[(i-yminId)*(xmaxId-xminId)+(j-xminId)].begin(), *it_nn);	
// 								}
							}
							//else if(  dist > distance_2neighbor[(i-yminOut)*(xmaxOut-xminOut)+(j-xminOut)] ){ //candidate of neighbor is closer than 2nd neighbor
							else if(  dist > distance_2neighbor.getValueAt(i, j) ){ //candidate of neighbor is closer than 2nd neighbor
								distance_2neighbor.setValueAt(i,j, dist);
								IDLocal[(i-yminId)*(xmaxId-xminId)+(j-xminId)].insert( ++IDLocal[(i-yminId)*(xmaxId-xminId)+(j-xminId)].begin() , *it_nn);							  
							}
							else { 
								// probably there are more than 3 grains nearer than DELTA to this gridpoint
								IDLocal[(i-yminId)*(xmaxId-xminId)+(j-xminId)].push_back(*it_nn);						  
							}
						}
					}
				}
			}
			
		}
		neighbors_2order.erase(it_nn);
	}
// 	plot_box(true,1,"Compare_1");
// 	plot_box(true,2,"Compare_2");
	char buffer;
// 
// 	  // checke schnitt zum randkorn:
	checkIntersect_zero_grain();
// 	plot_box(true,1,"Compare_1_zero");
// 	plot_box(true,2,"Compare_2_zero");

	// 	be careful for parralisation!!!!!

	set_comparison();
	plot_box(true,1,"Compare_1_set");
	plot_box(true,2,"Compare_2_set");

	
// 	cout << "comparison complete for grain: "<<id << endl;
// 	cin>> buffer;
	// 	write the compared values to the distanceBuffer1 array

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
		for (int i = inputDistance->getMinY(); i < inputDistance->getMaxY(); i++){
			for (int j = inputDistance->getMinX(); j < inputDistance->getMaxX(); j++){
				if ((i <= 2* grid_blowup) || (m-2*grid_blowup <= i) || (j <= 2*grid_blowup) || (m-2*grid_blowup <= j)){
					if(outputDistance->getValueAt(i,j) < boundary->outputDistance->getValueAt(i,j)){
						outputDistance->setValueAt(i,j,boundary->outputDistance->getValueAt(i,j));
					}
				}
			}
		}
	}
}


void LSbox::add_n2o(){
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
	switch_in_and_out();
	exist = false;
	
	// save old boundaries -> function will compute updates

    contourGrain.clear();
    int grid_blowup = handler->get_grid_blowup(); 
    double h = handler->get_h();
    int loop = handler->loop;
    int m = handler->get_ngridpoints();
	
	
    stringstream s;
    int first_i, first_j;
    int current_i, current_j;
    int next_i, next_j;
    char direction = 1; // x+
    // directions 0 = y-  //  2 = y+  //  3 x-  (y- = up // y + down; (0,0)left upper corner)
	double pointx; 
    double pointy;
    
	
	int dist = inputDistance->getMaxY() - inputDistance->getMinY();
	int i = inputDistance->getMinY() + int(dist/2);
    // look for distToZero in row y
    for (int j = inputDistance->getMinX(); j < inputDistance->getMaxX()-1; j++) {
        //if ((*inputDistance)[(i-yminIn)*(xmaxIn-xminIn)+(j-xminIn)]  * (*inputDistance)[(i-yminIn)*(xmaxIn-xminIn)+(j-xminIn+1)]  <= 0) {
        if (inputDistance->getValueAt(i,j)  * inputDistance->getValueAt(i,j+1)  <= 0) {
            first_i = i; 	first_j = j; 
            current_i = i; 	current_j = j;
	    next_i =i; 		next_j = j+1;
	    exist= true;
		cout << "zero found" << endl;
		double slope =  inputDistance->getValueAt(current_i, current_j) - inputDistance->getValueAt(next_i,next_j);		
		SPoint point;
		point.x= current_j + (inputDistance->getValueAt(current_i, current_j)/slope);
		point.y = current_i;
		 cout << current_i <<"   " << current_j << endl;
		 cout << point.y <<"   " << point.x << endl;
		 char buf;
		 cin>> buf;
		contourGrain.emplace_back(point);
	    break;
        }
	}
	if (!exist) {
		cout << "search in y-direction" << endl;
		int dist = inputDistance->getMaxX() - inputDistance->getMinX();
		int j = inputDistance->getMinX() + int(dist/2);
		for (int i = inputDistance->getMinY(); i < inputDistance->getMaxY()-1; i++) {
			//if ((*inputDistance)[(i-yminIn)*(xmaxIn-xminIn)+(j-xminIn)]  * (*inputDistance)[(i-yminIn+1)*(xmaxIn-xminIn)+(j-xminIn)]  <= 0) {
			if (inputDistance->getValueAt(i,j)  * inputDistance->getValueAt(i+1,j)  <= 0) {
				first_i = i; 	first_j = j; 
				current_i = i; 	current_j = j;
				next_i =i+1; 	next_j = j;
				exist= true;
				direction = 2;
				SPoint point;
				cout << "boundary found"<< endl;
				
				double slope =  inputDistance->getValueAt(next_i,next_j) - inputDistance->getValueAt(current_i, current_j);
				
				cout << slope << endl;
				slope = -1.0 *slope;
				point.x= current_j;
				point.y = current_i+ (inputDistance->getValueAt(current_i, current_j)/slope);
				
				contourGrain.emplace_back(point);
					cout << current_i <<"   " << current_j << endl;
					cout << point.y <<"   " << point.x << endl;
					char buf;
					cin>> buf;
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
    SPoint point;

    // 	 begin search
	int xmaxOut = 0, xminOut = m, ymaxOut = 0, yminOut = m;
    while (newZero) {		
	  current_i = next_i;
	  current_j = next_j;
	  
	  // check for size change
	  if (current_j-grid_blowup <= xminOut && current_j-grid_blowup >= 0) 		
		xminOut = current_j - grid_blowup;
	  else if (current_j+grid_blowup >= xmaxOut && current_j + grid_blowup <= m) 	
		xmaxOut = current_j + grid_blowup;
	  if (current_i <= yminOut+grid_blowup && current_i-grid_blowup >= 0) 		
		yminOut = current_i - grid_blowup;
	  else if (current_i+grid_blowup >= ymaxOut && current_i + grid_blowup <= m) 	
		ymaxOut = current_i + grid_blowup ;
	  
	  // change search directions
	  //(1 = left turn; -1 right turn)  

	  //sgn = utils::sgn((*inputDistance)[(current_i-yminIn)*(xmaxIn-xminIn)+(current_j-xminIn)]);
	  sgn = utils::sgn(inputDistance->getValueAt(current_i,current_j));
	  if (sgn == 0) { 
		if (direction == 0) 	{next_i = current_i-1;next_j = current_j;}
		else if (direction == 2)  {next_i = current_i+1;next_j = current_j;}
		else if (direction == 1)  {next_j = current_j+1;next_i = current_i;}
		else if (direction == 3)  {next_j = current_j-1;next_i = current_i;}
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
// 		  cout <<"   "<< current_i<<"   "<< current_j <<"   "<< next_i <<"   "<< next_j<< endl;
		  
//            cerr << current_i  <<"  "<< current_j <<"  "<< next_i <<"  "<< next_j <<endl ;
		  //if ( (*inputDistance)[(current_i-yminIn)*(xmaxIn-xminIn)+(current_j-xminIn)] * (*inputDistance)[(next_i-yminIn)*(xmaxIn-xminIn)+(next_j-xminIn)] <= 0) {
		  if ( inputDistance->getValueAt(current_i, current_j) * inputDistance->getValueAt(next_i, next_j) <= 0) {
			  foundnext = true;
			  //if( ((*inputDistance)[(current_i-yminIn)*(xmaxIn-xminIn)+(current_j-xminIn)] * (*inputDistance)[(next_i-yminIn)*(xmaxIn-xminIn)+(next_j-xminIn)] )!= 0.0 )
			  if( inputDistance->getValueAt(current_i, current_j) * inputDistance->getValueAt(next_i, next_j) != 0.0 )
			  {
					//double slope =  (*inputDistance)[(next_i-yminIn)*(xmaxIn-xminIn)+(current_j-xminIn)] - (*inputDistance)[(current_i-yminIn)*(xmaxIn-xminIn)+(current_j-xminIn)];
					double slope =  inputDistance->getValueAt(current_i, current_j) - inputDistance->getValueAt(next_i,next_j) ;
// 					cout << slope << endl;
					
					if (direction == 1) {	 
					//point.x= current_j + (*inputDistance)[(current_i-yminIn)*(xmaxIn-xminIn)+(current_j-xminIn)]/slope;
					point.x= current_j + (inputDistance->getValueAt(current_i, current_j)/slope);
					point.y = current_i;
					}
					else if (direction == 3){	 
					//point.x = current_j - (*inputDistance)[(current_i-yminIn)*(xmaxIn-xminIn)+(current_j-xminIn)]/slope;
						point.x = current_j - inputDistance->getValueAt(current_i, current_j)/slope;
					point.y = current_i;
					}
					else if (direction == 0) {	 				
					//point.y =  current_i - (*inputDistance)[(current_i-yminIn)*(xmaxIn-xminIn)+(current_j-xminIn)]/slope;
						point.y =  current_i - inputDistance->getValueAt(current_i, current_j)/slope;
					point.x = current_j;
					}
					else if (direction == 2){	
					//point.y =  current_i + (*inputDistance)[(current_i-yminIn)*(xmaxIn-xminIn)+(current_j-xminIn)]/slope;
					point.y =  current_i + inputDistance->getValueAt(current_i, current_j)/slope;
					point.x = current_j;
					}		 
			  }
			  else {
				cerr << "levelset on gridpoint  " << current_i << "\t" << current_j<<"\t" << inputDistance->getValueAt(current_i, current_j)<<"\t" << inputDistance->getValueAt(next_i, current_j)<< endl;
				
				//if ((*inputDistance)[(current_i-yminIn)*(xmaxIn-xminIn)+(current_j-xminIn)] == 0.0) { point.y =  current_i;  point.x = current_j; }
				if (inputDistance->getValueAt(current_i, current_j) == 0.0)
					{ point.y =  current_i;  point.x = current_j; }
				//else if ((*inputDistance)[(next_i-yminIn)*(xmaxIn-xminIn)+(current_j-xminIn)]  == 0.0)
				else if (inputDistance->getValueAt(next_i, current_j)  == 0.0)
					{ point.y =  next_i;  point.x = next_j; }
			  }
			  cout << "new zero  " << point.y <<"   " << point.x << endl;
			  contourGrain.emplace_back(point);	

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
	  // check if completed contourline
	  if (current_j == first_j && current_i == first_i) {
		  newZero = false;	//springen aus der while-schleife, wir haben eine geschlosse kurve gefunden
	  }
    }
    cout << " contour found!! " << endl;
	
    if (xminOut < 0) {cout <<"undefined box size for xmin: "<< xminOut << endl; abort();} //xmin = 0;
	if (xmaxOut > m) {cout <<"undefined box size for xmax: "<< xmaxOut << endl; abort();} //xmax = m;
	if (yminOut < 0) {cout <<"undefined box size for ymin: "<< yminOut << endl; abort();} //ymin = 0;
	if (ymaxOut > m) {cout <<"undefined box size for ymax: "<< xmaxOut << endl; abort();} // ymax = m;
	
	
    outputDistance->resize(xminOut, yminOut, xmaxOut, ymaxOut);
	outputDistance->resizeToSquare(handler->get_ngridpoints());
    
    // compute Volume and Energy
    if ( (loop % int(ANALYSESTEP)) == 0 || loop == TIMESTEPS ) {
		vector<SPoint>::iterator volumeit=contourGrain.begin();
		energy = 0;
		volume = 0;
		double px, py;
		double gamma_hagb = 0.6;
		double theta_ref = 15.0* PI / 180.;
		double theta_mis;
		
		px= (*volumeit).x;
		py= (*volumeit).y;
		cout << py << "  " << px << endl;
		cout << (*(++volumeit)).y  << "  " << (*(++volumeit)).x << endl;
		volumeit++;
		for (; volumeit!= contourGrain.end(); volumeit++){
			s << (*volumeit).x << "\t" << (*volumeit).y<<endl;      
			volume += (py+(*volumeit).y)*(px-(*volumeit).x);
			px= (*volumeit).x;
			py= (*volumeit).y;
// 			cout << px << "  " << py << endl;
// 			theta_mis=mis_ori( IDLocal[ ( (int(py+0.5)-yminId) * (xmaxId-xminId)) + (int(px+0.5) - xminId) ][0] ); //find the LSbox pointer to the next neighbor -> therefor find the next grid pointer
// 			if (theta_mis <= theta_ref)	energy += h* gamma_hagb * ( theta_mis / theta_ref) * (1.0 - log( theta_mis / theta_ref));
// 				else energy += h* gamma_hagb;
		}		
		cout << (*(--volumeit)).y  << "  " << (*(--volumeit)).x << endl;
		cout << py << "  " << px << endl;
		
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


/**************************************/
//  Redistancing
/**************************************/

void LSbox::redist_box() {	
	int grid_blowup = handler->get_grid_blowup(); 
	double h = handler->get_h();
// 	plot_box(false);

	double slope = 1;
	double candidate, i_slope, distToZero;
	
	//resizeToSquareOut();
	//TODO: Check this

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
			if(j <= intersec_xmax && i <= intersec_ymax){
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
			if(j >= intersec_xmin && i <= intersec_ymax){
// 				if (inputDistance->getValueAt(i,j) * inputDistance->getValueAt(i,j-1) <= 0.0) {
// 					// interpolate
// 					i_slope  = ( inputDistance->getValueAt(i,j-1)  - inputDistance->getValueAt(i,j) ) / h;
// 					distToZero = - inputDistance->getValueAt(i,j) / i_slope;
// 					if ( abs(inputDistance->getValueAt(i,j) ) > abs(distToZero))
// 						inputDistance->setValueAt(i,j, -distToZero * utils::sgn(i_slope));
// 					}
					// calculate new distance candidate and assign if appropriate
					candidate = outputDistance->getValueAt(i,j)  + (utils::sgn( inputDistance->getValueAt(i,j-1) ) * h);
					if (abs(candidate) < abs(outputDistance->getValueAt(i,j-1)))
						outputDistance->setValueAt(i,j-1, candidate);
			}
			else {
				candidate = outputDistance->getValueAt(i,j)  + utils::sgn( inputDistance->getValueAt(i,j-1))*h;
				if (abs(candidate) < abs(outputDistance->getValueAt(i,j-1)))
					outputDistance->setValueAt(i,j-1, candidate);
			}
		}		
	}

		
	// y-direction forward
	for (int j = intersec_xmin; j < outputDistance->getMaxX(); j++) {
		for (int i = intersec_ymin; i < outputDistance->getMaxY()-1; i++) {
			if(j <= intersec_xmax && i <= intersec_ymax){
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
				candidate = outputDistance->getValueAt(i,j)  -h;//+ (utils::sgn( inputDistance->getValueAt(i+1,j)) * h);
				if (abs(candidate) < abs(outputDistance->getValueAt(i+1,j)))
					outputDistance->setValueAt(i+1,j, candidate);
			}
		}		
	}

	for (int j = intersec_xmin; j < outputDistance->getMaxX(); j++) {
		for (int i = intersec_ymax-1; i > outputDistance->getMinY(); i--) {
			if(j <= intersec_xmax && i >= intersec_ymin){
				if (inputDistance->getValueAt(i,j) * inputDistance->getValueAt(i-1,j) <= 0.0) {
					// interpolate
					i_slope  = (inputDistance->getValueAt(i-1,j)  - inputDistance->getValueAt(i,j))/ h;
					distToZero = - inputDistance->getValueAt(i,j) / i_slope;
					if ( abs(outputDistance->getValueAt(i,j) ) > abs(distToZero))
						outputDistance->setValueAt(i,j, -distToZero * utils::sgn(i_slope));
				}
				// calculate new distance candidate and assign if appropriate
				candidate = outputDistance->getValueAt(i,j)  + (utils::sgn( inputDistance->getValueAt(i-1, j) ) * h);
				if (abs(candidate) < abs(outputDistance->getValueAt(i-1, j) ))
					outputDistance->setValueAt(i-1, j, candidate);
			}
			else {
				candidate = outputDistance->getValueAt(i,j)  -h ; //+ (utils::sgn( outputDistance->getValueAt(i-1, j) ) * h);
				if (abs(candidate) < abs(outputDistance->getValueAt(i-1,j)))
					outputDistance->setValueAt(i-1, j, candidate);
			}
		}		
	}

	plot_box(true,1,"Redist_1");
	plot_box(true,2,"Redist_2");
	//TODO: Analyze this
		char buf;
	cin>>buf;
	
	inputDistance->resize(outputDistance->getMinX(), outputDistance->getMinY(), outputDistance->getMaxX(), outputDistance->getMaxY());

	
	// 	 set the references for the convolution step

}

/**************************************/
// end of redist
/**************************************/



/**************************************/
// plot the box and all its properties
/**************************************/

void LSbox::plot_box_contour(int loop)
{

  stringstream filename;
  filename<< "TempBox_"<< id <<"_T"<< loop << ".gnu";
  
  ofstream datei;
  datei.open(filename.str());
  vector<SPoint>::iterator contourIterator;

    for (contourIterator= contourGrain.begin(); contourIterator != contourGrain.end(); contourIterator++){
	datei << (*contourIterator).x << "\t" << (*contourIterator).y<< endl;
    }
    datei << endl;
  
datei.close();

  
}

void LSbox::plot_box(bool distanceplot, int select, string simstep){
	cout <<" \nGrain  Info: " << endl;
	cout << " ID :" <<id << endl;
    cout << " xminIn, xmaxIn, yminIn, ymaxIn :" << inputDistance->getMinX() << " || "<< inputDistance->getMaxX() << " || " << inputDistance->getMinY() << " || " << inputDistance->getMaxY() << endl;
	 cout << " xminOut, xmaxOut, yminOut, ymaxOut :" << outputDistance->getMinX() << " || "<< outputDistance->getMaxX() << " || " << outputDistance->getMinY() << " || " << outputDistance->getMaxY() << endl;
	  cout << " xminId, xmaxId, yminId, ymaxId :" << xminId << " || "<< xmaxId << " || " << yminId << " || " << ymaxId << endl;
//     if (distanceplot==true) utils::print_2dim_array(distance,ymax-ymin,xmax-xmin);
// 		else cout << " no distance values in storage!" << endl;

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
		filename<< "BoxDistance_"<< simstep << "_" << id << ".gnu";
		ofstream datei;
		datei.open(filename.str());
		
		if(select == 2) {
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
// 	here we could work direktly with quarternions
	return misorientationCubic(phi1,PHI,phi2,grain_2->get_phi1(), grain_2->get_PHI(), grain_2->get_phi2());
}



void LSbox::shape_distance(){
	//TODO: WORK ON THIS
	int m = outputDistance->getMaxX()-outputDistance->getMinX();
	for (int i = 0; i < outputDistance->getMaxY()-outputDistance->getMinY(); i++) {
		for (int j = 0; j < m; j++) {
			outputDistance->getRawData()[i*m +j] *= -4.0;
		}
	}
}


