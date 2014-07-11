#include "box.h"
#include "Settings.h"
#include "dimensionalBufferReal.h"
#include "pooledDimensionalBufferReal.h"
#include "contourSector.h"
#include "grahamScan.h"

LSbox::LSbox() :exist(false), quaternion(NULL), inputDistance(NULL), outputDistance(NULL), local_weights(NULL){  }

LSbox::~LSbox() {
	if(quaternion!=NULL) delete [] quaternion;
	delete inputDistance;
	delete outputDistance;
	if(local_weights!=NULL) delete local_weights;
}


LSbox::LSbox(int id, double phi1, double PHI, double phi2, grainhdl* owner) :
		id(id), handler(owner)
{	
	exist = true;
	quaternion = new double[4];
	double euler[3] = {phi1,PHI,phi2};
	(*(handler->mymath)).euler2quaternion( euler, quaternion );

	inputDistance = new DimensionalBufferReal(0, 0, 0, 0, 0, 0);
	outputDistance = new DimensionalBufferReal(0, 0, 0, 0, 0, 0);

	local_weights = new Weightmap(owner);
	boundaryGrain = false;
}

LSbox::LSbox(int aID, voro::voronoicell_neighbor& c, double *part_pos, grainhdl* owner) : id(aID), nvertices(0), handler(owner) {
	
	int grid_blowup = owner->get_grid_blowup(); 
	exist = true;
	double h = owner->get_h();
    // determine size of grain
	quaternion = new double[4];
	
		if(Settings::UseTexture){
		double newOri[3];
		(*(handler->mymath)).newOrientationFromReference( handler->bunge, handler->deviation, newOri );
		(*(handler->mymath)).euler2quaternion(newOri, quaternion);
	}
	else (*(handler->mymath)).randomOriShoemakeQuat( quaternion );
	
    int xmax = 0; 
	int xmin = handler->get_ngridpoints(); 
	int ymax = 0; 
	int ymin = xmin;
	int zmax = 0;
	int zmin = xmin;
	
    vektor x1(2), x2(2), x3(2);
    vector<double> vv;
    exist = true;
	c.vertices(part_pos[3*(id-1)],part_pos[3*(id-1)+1],part_pos[3*(id-1)+2],vv);

//	TODO how to genrate the 3D Hull or DISTANCEFUNCTION????
    GrahamScan scanner(c, id, part_pos);
//    TODO contourGrain??
    scanner.generateCovnexHull(contourGrain);

    double x, y, z;
    for (unsigned int k=0; k < contourGrain.size(); k++){
//   TODO 	perhaps compare against all points defining th c object
    		y=contourGrain[k].y;
    		x=contourGrain[k].x;
    		z=contourGrain[k].z;
    		if (y/h < ymin) ymin = y/h;
    		if (y/h > ymax) ymax = y/h+1;
    		if (x/h < xmin) xmin = x/h;
    		if (x/h > xmax) xmax = x/h+1;
    		if (z/h < zmin) zmin = z/h;
    		if (z/h > zmax) zmax = z/h+1;
        }
    xmax += 2*grid_blowup;
    ymax += 2*grid_blowup;
    zmax += 2*grid_blowup;

    inputDistance = new DimensionalBufferReal(xmin, ymin, xmax, ymax, zmin, zmax);
	outputDistance = new DimensionalBufferReal(xmin, ymin, xmax, ymax, zmin, zmax);
	
 	inputDistance->resizeToCube(handler->get_ngridpoints());
 	outputDistance->resizeToCube(handler->get_ngridpoints());
	inputDistance->clearValues(0.0);
	outputDistance->clearValues(0.0);
	
	get_new_IDLocalSize();
	IDLocal.resize(xminId, yminId, xmaxId, ymaxId, zminId, zmaxId);
	
	local_weights=new Weightmap(owner);
	boundaryGrain = false;
// 	cout << "made a new box: xmin="<<xmin<< " xmax="<<xmax <<" ymin="<<ymin << " ymax="<<ymax<<endl;
}


LSbox::LSbox(int id, int nvertices, double* vertices, double q1, double q2, double q3, double q4, grainhdl* owner) : id(id), nvertices(nvertices), handler(owner){
	exist = true;
	quaternion = new double[4];
	quaternion[0]=q1 ;
	quaternion[1]=q2 ;
	quaternion[2]=q3 ;
	quaternion[3]=q4 ;
	contourGrain.resize(nvertices);

	int grid_blowup = owner->get_grid_blowup();
	double h = owner->get_h();
    // determine size of grain
    int xmax = 0;
	int xmin = handler->get_ngridpoints();
	int ymax = 0;
	int ymin = xmin;
	int zmax = 0;
	int zmin = xmin;

	double y,x,z;
	exist = true;

	for (unsigned int k=0; k < nvertices; k++){
		z=vertices[(2*k)+2];
		y=vertices[(2*k)+1];
		x=vertices[2*k];
		contourGrain[k].x =vertices[2*k];
		contourGrain[k].y =vertices[2*k +1];
		contourGrain[k].z =vertices[2*k +2];
		if (y/h < ymin) ymin = y/h;
		if (y/h > ymax) ymax = y/h+1;
		if (x/h < xmin) xmin = x/h;
		if (x/h > xmax) xmax = x/h+1;
		if (z/h < zmin) zmin = z/h;
		if (z/h > zmax) zmax = z/h+1;
    }
	xmax += 2*grid_blowup;
	ymax += 2*grid_blowup;
    zmax += 2*grid_blowup;

    inputDistance = new DimensionalBufferReal(xmin, ymin, xmax, ymax, zmin, zmax);
	outputDistance = new DimensionalBufferReal(xmin, ymin, xmax, ymax, zmin, zmax);

	
 	inputDistance->resizeToCube(handler->get_ngridpoints());
 	outputDistance->resizeToCube(handler->get_ngridpoints());
	inputDistance->clearValues(0.0);
	outputDistance->clearValues(0.0);

	get_new_IDLocalSize();
	IDLocal.resize(xminId, yminId, xmaxId, ymaxId, zminId, zmaxId);

// 	cout << "made a new box: xmin="<<xmin<< " xmax="<<xmax <<" ymin="<<ymin << " ymax="<<ymax<<endl;
	local_weights=new Weightmap(owner);
	boundaryGrain = false;
}


int LSbox::getID() {
    return id;
}

bool LSbox::isNeighbour(LSbox* candidate){
	vector<characteristics>::iterator it;
	for(it = grainCharacteristics.begin(); it!=grainCharacteristics.end(); it++){
		if((*it).directNeighbour==candidate) return true;
	}
	return false;
}


double LSbox::mis_ori(LSbox* grain_2){
	if(get_status() != true ) {
		cout << "try to compute misori for are disappeared grains :   " ;
		cout << "grains are: " << id << "  " << grain_2->get_id() << endl;
// 		char buf; cin >> buf;
	}
// 	here we could work direktly with quarternions
	return (*(handler->mymath)).misorientationCubicQxQ( quaternion[0], quaternion[1], quaternion[2], quaternion[3], grain_2->quaternion[0], grain_2->quaternion[1], grain_2->quaternion[2], grain_2->quaternion[3] );
}



void LSbox::get_new_IDLocalSize(){
	xmaxId = outputDistance->getMaxX();
	xminId = outputDistance->getMinX();
	ymaxId = outputDistance->getMaxY();
	yminId = outputDistance->getMinY();
	zmaxId = outputDistance->getMaxZ();
	zminId = outputDistance->getMinZ();
}

void LSbox::executeFFTW(fftw_plan fftplan)
{
  fftw_execute(fftplan);
}

void LSbox::executeFFTW(fftwf_plan fftplan)
{
  fftwf_execute(fftplan);
}


void LSbox::destroyFFTWs(fftw_plan fwdPlan, fftw_plan bwdPlan){
	  fftw_destroy_plan(fwdPlan);
	  fftw_destroy_plan(bwdPlan);
}

void LSbox::destroyFFTWs(fftwf_plan fwdPlan, fftwf_plan bwdPlan){
	  fftwf_destroy_plan(fwdPlan);
	  fftwf_destroy_plan(bwdPlan);
}


void LSbox::makeFFTPlans(double *in, double* out,fftw_complex *fftTemp, fftw_plan *fftplan1, fftw_plan *fftplan2)
{ /* creates plans for FFT and IFFT */
	int n = outputDistance->getMaxX() - outputDistance->getMinX();
	int m = outputDistance->getMaxY() - outputDistance->getMinY();
	int l = outputDistance->getMaxZ() - outputDistance->getMinZ();

	*fftplan1 = fftw_plan_dft_r2c_3d(n,m,l,in,fftTemp,FFTW_ESTIMATE);
	*fftplan2 = fftw_plan_dft_c2r_3d(n,m,l,fftTemp,out,FFTW_ESTIMATE);
	/*
	The flags argument is usually either FFTW_MEASURE or FFTW_ESTIMATE. FFTW_MEASURE
	instructs FFTW to run and measure the execution time of several FFTs in order to find the
	best way to compute the transform of size n. This process takes some time (usually a few
	seconds), depending on your machine and on the size of the transform. FFTW_ESTIMATE,
	on the contrary, does not run any computation and just builds a reasonable plan that is
	probably sub-optimal. In short, if your program performs many transforms of the same size
	and initialization time is not important, use FFTW_MEASURE; otherwise use the estimate. */
}

void LSbox::makeFFTPlans(float *in, float* out,fftwf_complex *fftTemp, fftwf_plan *fftplan1, fftwf_plan *fftplan2)
{ /* creates plans for FFT and IFFT */
	int n = outputDistance->getMaxX() - outputDistance->getMinX();
	int m = outputDistance->getMaxY() - outputDistance->getMinY();
	int l = outputDistance->getMaxZ() - outputDistance->getMinZ();

	*fftplan1 = fftwf_plan_dft_r2c_3d(n,m,l,in,fftTemp,FFTW_ESTIMATE);
	*fftplan2 = fftwf_plan_dft_c2r_3d(n,m,l,fftTemp,out,FFTW_ESTIMATE);
}


//TODO 3D Convolution Generator
void LSbox::conv_generator(fftwp_complex *fftTemp, fftwp_plan fftplan1, fftwp_plan fftplan2)
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
	double dt = handler->get_dt();
	int n2 = floor(n/2) + 1;
	int nn = (*handler).get_realDomainSize();
	double nsq =  nn*nn;
	double k = 2.0 * PI / n;
	double G;
	double coski;
	int j2;
	int i2;

	executeFFTW(fftplan1);
	//	Forward DFT

	switch (Settings::ConvolutionMode){
		case E_LAPLACE : {
//			for(int i=0;i<n2;i++) {
//				coski=cos(k*i);
//				for(int j=0;j<n;j++){
//					G = 2.0*(2.0 - coski - cos(k*j)) * nsq;
//					G = 1.0/ (1.0+(dt*G)) / (n*n);
//					fftTemp[i+n2*j][0] = fftTemp[i+n2*j][0]*G;
//					fftTemp[i+n2*j][1] = fftTemp[i+n2*j][1]*G;
//				}
//			}
			break;
		}
		case E_LAPLACE_RITCHARDSON : {
//			Ritchardson Extrapolation
//			for(int i=0;i<n2;i++) {
//				coski=cos(k*i);
//				for(int j=0;j<n;j++){
//					G = 2.0*(2.0 - coski - cos(k*j)) * nsq;
//					G = (4.0/pow(1+1.5*(dt)/40*G,40) - 1.0 / pow(1+3.0*(dt)/40*G,40)) / 3.0 / (double)(n*n);
//					fftTemp[i+n2*j][0] = fftTemp[i+n2*j][0]*G;
//					fftTemp[i+n2*j][1] = fftTemp[i+n2*j][1]*G;
//				}
//			}
			break;
		}
		case E_GAUSSIAN : {
			nsq=n*n*n;
//			Convolution with Normaldistribution
			for(int k=0;k<n2;k++) {
				k2 = mymin(k,n-k);
				for(int i=0;i<n2;i++) {
					i2 = mymin(i,n-i);
					for(int j=0;j<n;j++){
						j2 = mymin(j,n-j);
						G = exp(-(static_cast<double>(i2*i2+j2*j2+k2*k2))*8.0*dt*PI*PI*PI) / nsq;
						//G = exp(-(static_cast<double>(i2*i2+j2*j2))*4.0*dt*PI*PI) / nsq;
						fftTemp[i+n2*j][0] = fftTemp[i+n2*j][0]*G;
						fftTemp[i+n2*j][1] = fftTemp[i+n2*j][1]*G;
					}
				}
			}
			break;
		}
	}

	executeFFTW(fftplan2);
	//	Inverse DFT
}




//TODO 3D
void LSbox::distancefunction(){

	int grid_blowup = handler->get_grid_blowup();
	double h = handler->get_h();
	int i=0,j=0,k=0;
	SPoint to_test;
	int contour_size = contourGrain.size();
    double* constant = new double[contour_size];
    double* multiple = new double[contour_size];

    //Precalculate values to speed up the PIP iterations
    for (int i = 1; i < contour_size; i++)
	{
		if (contourGrain[i].x == contourGrain[i].y)
		{
			constant[i] = contourGrain[i].x;
			multiple[i] = 0;
		}
		else
		{
			constant[i] = contourGrain[i].x
					- (contourGrain[i].y * contourGrain[j].x) / (contourGrain[j].y - contourGrain[i].y)
					+ (contourGrain[i].y * contourGrain[i].x) / (contourGrain[j].y - contourGrain[i].y);
			multiple[i] = (contourGrain[j].x - contourGrain[i].x) / (contourGrain[j].y - contourGrain[i].y);
		}
		j = i;
	}

	for (i = outputDistance->getMinY(); i < outputDistance->getMaxY(); i++)
	{
		for (j = outputDistance->getMinX(); j < outputDistance->getMaxX(); j++)
		{
			to_test.x = (j - grid_blowup) * h;
			to_test.y = (i - grid_blowup) * h;

			bool isInside = false;

			for (int k = 1, l =0; k < contour_size; k++)
			{
				if ((contourGrain[k].y < to_test.y && contourGrain[l].y > to_test.y) ||
					(contourGrain[l].y < to_test.y && contourGrain[k].y > to_test.y))
				{
					bool new_val = (to_test.y * multiple[k] + constant[k] < to_test.x);
					isInside = (isInside != new_val);	//isInside = isInside XOR new_val
				}
				l = k;
			}

			double minDist = 1000000.0;
			for(int k=1, l=0; k<contour_size; k++)
			{
				SPoint	u = contourGrain[k]-contourGrain[l];
				double lambda = (to_test - contourGrain[l]).dot(u);
				lambda /= u.dot(u);

				double dist;
				if(lambda < 0)
				{
					dist = (to_test - contourGrain[l]).len();
				}
				else if (lambda > 1)
				{
					dist = (contourGrain[k] - to_test).len();
				}
				else
				{
					dist = (to_test - (contourGrain[l] + u*lambda) ).len();
				}
				minDist = min(minDist,dist);
				l = k;
			}
			if(minDist > handler->delta)
				minDist = handler->delta;
			outputDistance->setValueAt(i,j, isInside ? minDist : -minDist);
		}
	}

	delete [] constant;
	delete [] multiple;
//	plot_box(true,2,"Dist",true);
}

/*************************************/
void LSbox::switchInNOut(){


	DimensionalBufferReal* temp;

	temp = inputDistance;
	inputDistance = outputDistance;
	outputDistance = temp;
}


void LSbox::add_n2o_2(){
	if(get_status() != true ) return;
	vector<characteristics> ::iterator it, it_ngC;
	vector<LSbox*> ::iterator it_com, it_nC;
	bool just_in;
	neighbors_2order.clear();
	for(it = grainCharacteristics.begin(); it != grainCharacteristics.end(); it++){
		if((*it).directNeighbour->get_status()==true)
			neighbourCandidates.push_back((*it).directNeighbour);
		for( it_ngC = (*it).directNeighbour->grainCharacteristics.begin(); it_ngC != (*it).directNeighbour->grainCharacteristics.end(); it_ngC++){
			if((*it_ngC).directNeighbour->get_status()==true)
				if(checkIntersect((*it_ngC).directNeighbour)){
					neighbourCandidates.push_back((*it_ngC).directNeighbour);
				}
		}
	}
	for(it_nC = neighbourCandidates.begin(); it_nC != neighbourCandidates.end(); it_nC++){
		just_in=false;
		if((*it_nC)==this) continue;
		if((*it_nC)==handler->boundary) continue;
		for(it_com= neighbors_2order.begin(); it_com != neighbors_2order.end(); it_com++){
			if((*it_com)==(*it_nC)) {
				just_in = true;
				break;
			}
		}
		if((!just_in)) neighbors_2order.push_back((*it_nC));
	}
	neighbourCandidates.clear();
}



// Comparison + Helperfunctions
/**************************************/
/**************************************/

double LSbox::getDistance(int i, int j, int k){
   return inputDistance->getValueAt(i,j,k);
}


void LSbox::set_comparison(){
	for (int k = outputDistance->getMinZ(); k < outputDistance->getMaxZ(); k++){
		for (int i = outputDistance->getMinY(); i < outputDistance->getMaxY(); i++){
			for (int j = outputDistance->getMinX(); j < outputDistance->getMaxX(); j++){

				if(abs(inputDistance->getValueAt(i,j,k)) < 0.7*handler->delta ) {
					outputDistance->setValueAt(i, j, k, 0.5 * (inputDistance->getValueAt(i,j,k) - outputDistance->getValueAt(i,j,k)));
				}
				else outputDistance->setValueAt(i, j, k, inputDistance->getValueAt(i,j,k));
			}
		}
	}
}


bool LSbox::checkIntersect(LSbox* box2) {
    if (inputDistance->getMinX() > box2->inputDistance->getMaxX() ||
    	inputDistance->getMaxX() < box2->inputDistance->getMinX() ||
    	inputDistance->getMinY() > box2->inputDistance->getMaxY() ||
    	inputDistance->getMaxY() < box2->inputDistance->getMinY() ||
    	inputDistance->getMinZ() > box2->inputDistance->getMaxZ() ||
    	inputDistance->getMaxZ() < box2->inputDistance->getMinZ())
    	return false;
    return true;
}


void LSbox::comparison(ExpandingVector<char>& mem_pool){
	if(get_status() != true ) return;
	int m = handler->get_ngridpoints();
	int grid_blowup = handler->get_grid_blowup();

	mem_pool.expand( (outputDistance->getMaxX() - outputDistance->getMinX()) *
					 (outputDistance->getMaxY() - outputDistance->getMinY()) *
					 (outputDistance->getMaxZ() - outputDistance->getMinZ()) * sizeof(double) );
	PooledDimensionalBuffer distance_2neighbor(&mem_pool[0], mem_pool.size(),
			outputDistance->getMinX(), outputDistance->getMinY(),
			outputDistance->getMaxX(), outputDistance->getMaxY(),
			outputDistance->getMinZ(), outputDistance->getMaxZ());


	distance_2neighbor.clearValues(-1.0);
	outputDistance->clearValues(-1.0);

	int loop = handler->loop;
	std::vector<LSbox*>::iterator it_nn;

	for(it_nn = neighbors_2order.begin(); it_nn != neighbors_2order.end(); it_nn++){
		int x_min_new, x_max_new, y_min_new, y_max_new, z_min_new, z_max_new;

		if(inputDistance->getMinX() < (**it_nn).inputDistance->getMinX()) x_min_new = (**it_nn).inputDistance->getMinX();
			else x_min_new = inputDistance->getMinX();

		if(inputDistance->getMaxX() > (**it_nn).inputDistance->getMaxX()) x_max_new = (**it_nn).inputDistance->getMaxX();
			else x_max_new = inputDistance->getMaxX();

		if(inputDistance->getMinY() < (**it_nn).inputDistance->getMinY()) y_min_new = (**it_nn).inputDistance->getMinY();
			else y_min_new = inputDistance->getMinY();

		if(inputDistance->getMaxY() > (**it_nn).inputDistance->getMaxY()) y_max_new = (**it_nn).inputDistance->getMaxY();
			else y_max_new = inputDistance->getMaxY();

		if(inputDistance->getMinZ() < (**it_nn).inputDistance->getMinZ()) z_min_new = (**it_nn).inputDistance->getMinZ();
					else y_min_new = inputDistance->getMinY();

		if(inputDistance->getMaxZ() > (**it_nn).inputDistance->getMaxZ()) z_max_new = (**it_nn).inputDistance->getMaxZ();
			else z_max_new = inputDistance->getMaxZ();

		for (int k = z_min_new; k < z_max_new; k++){
			for (int i = y_min_new; i < y_max_new; i++){
				for (int j = x_min_new; j < x_max_new; j++){
					if(abs(inputDistance->getValueAt(i,j,k)) < handler->delta ) {
						double dist = (**it_nn).getDistance(i,j,k);
						if(abs(dist) < handler->delta ) {
							if( dist > outputDistance->getValueAt(i,j,k) ){
								if( IDLocal.getValueAt(i,j,k).total_elements == 0 ){
									distance_2neighbor.setValueAt(i,j,k,outputDistance->getValueAt(i,j,k));
								}
								outputDistance->setValueAt(i, j, k, dist);
								IDLocal.getValueAt(i,j,k).insertAtPosition(E_FIRST_POSITION, *it_nn);
							}
							else if(  dist > distance_2neighbor.getValueAt(i, j ,k) ){ //candidate of neighbor is closer than 2nd neighbor
								distance_2neighbor.setValueAt(i,j,k, dist);
								IDLocal.getValueAt(i,j,k).insertAtPosition(E_SECOND_POSITION, *it_nn);
							}
							else {
								IDLocal.getValueAt(i,j,k).insertAtPosition(E_LAST_POSITION, *it_nn);
							}
						}
					}
				}
			}
		}
	}
	if (BoundaryIntersection()){
		boundaryGrain = true;
		boundaryCondition();
	}
	else boundaryGrain = false;
	set_comparison();
//	plot_box(true,2,"Com",true);
//	if(id==201) plot_box(true,2,"Combig",false);

}


bool LSbox::BoundaryIntersection(){
  int xMinBoundary = handler->get_grid_blowup()+ handler->getBoundaryGrainTube();
  int yMinBoundary = xMinBoundary;
  int zMinBoundary = xMinBoundary;

  int xMaxBoundary = handler->get_ngridpoints() - handler->get_grid_blowup() - handler->getBoundaryGrainTube();
  int yMaxBoundary = xMaxBoundary;
  int zMaxBoundary = xMaxBoundary;

  if(outputDistance->getMinX() > xMinBoundary &&  outputDistance->getMaxX() < xMaxBoundary &&
     outputDistance->getMinY() > yMinBoundary &&  outputDistance->getMaxY() < yMaxBoundary &&
     outputDistance->getMinZ() > zMinBoundary &&  outputDistance->getMaxZ() < zMaxBoundary )
	return false;
  else return true;
}


//TODO comparison to boundary
void LSbox::boundaryCondition(){
	LSbox* boundary = handler->boundary;
	int grid_blowup = handler->get_grid_blowup();
	double h = handler->get_h();
	int m = handler->get_ngridpoints();
	int distXMin, distXMax, distX;
	int distYMin, distYMax, distY;
	double dist;
	int xMin, xMax, yMin, yMax;

	for (int i=inputDistance->getMinY(); i< inputDistance->getMaxY() ; i++ ){
		for (int j=inputDistance->getMinX(); j<inputDistance->getMaxX() ; j++ ){
			distXMin	= 	-(j-grid_blowup);
			distYMin 	= 	-(i-grid_blowup);
			distXMax	=	(j-(m-grid_blowup));
			distYMax 	= 	(i-(m-grid_blowup));

			if(abs(distXMin) < abs(distXMax)) distX = distXMin;
			  else distX = distXMax;
			if(abs(distYMin) < abs(distYMax)) distY = distYMin;
			  else distY = distYMax;

			if(distX > 0 && distY > 0)
				dist = sqrt((double)distX*distX + distY*distY);

			else if(distX < 0 && distY > 0)
				dist = distY;
			else if (distX > 0 && distY < 0)
				dist = distX;
			else if(distX < 0 && distY < 0)
				dist=max(distX,distY);
			else if (distX == 0){
				if (distY == 0) dist = 0;
				else if(distY < 0) dist = 0;
				else if(distY > 0) dist = distY;
			}
			else if (distY == 0){
				if(distX < 0) dist = 0;
				else if(distX > 0) dist =distX;
			}

			  if( dist*h > outputDistance->getValueAt(i,j) ){
				  outputDistance->setValueAt(i, j, dist*h);
				  IDLocal.getValueAt(i,j).insertAtPosition(E_FIRST_POSITION, boundary);
			  }
			  else {
				  IDLocal.getValueAt(i,j).insertAtPosition(E_LAST_POSITION, boundary);
			  }
//			}
		}
	}
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
	int intersec_xmin, intersec_xmax, intersec_ymin, intersec_ymax, intersec_zmin, intersec_zmax;

	if (inputDistance->getMinX() < outputDistance->getMinX())
		intersec_xmin = outputDistance->getMinX();
	else  intersec_xmin = inputDistance->getMinX();

	if (inputDistance->getMinY() < outputDistance->getMinY())
		intersec_ymin = outputDistance->getMinY();
	else  intersec_ymin = inputDistance->getMinY();

	if (inputDistance->getMinZ() < outputDistance->getMinZ())
		intersec_zmin = outputDistance->getMinZ();
	else  intersec_zmin = inputDistance->getMinZ();

	if (inputDistance->getMaxX() < outputDistance->getMaxX())
		intersec_xmax= inputDistance->getMaxX();
	else  intersec_xmax = outputDistance->getMaxX();

	if (inputDistance->getMaxY() < outputDistance->getMaxY())
		intersec_ymax= inputDistance->getMaxY();
	else  intersec_ymax = outputDistance->getMaxY();

	if (inputDistance->getMaxZ() < outputDistance->getMaxZ())
		intersec_zmax= inputDistance->getMaxZ();
	else  intersec_zmax = outputDistance->getMaxZ();



	for (int k = intersec_zmin; k < outputDistance->getMaxZ()-1; k++) {
		for (int i = intersec_ymin; i < outputDistance->getMaxY(); i++){
		  for (int j = intersec_xmin; j < outputDistance->getMaxX()-1; j++) {
				// x-direction forward
				if(j < intersec_xmax-1 && i < intersec_ymax ){
					if (inputDistance->getValueAt(i,j,k) * inputDistance->getValueAt(i,j+1,k) <= 0.0) {
						// interpolate
						i_slope  = ( inputDistance->getValueAt(i,j+1) - inputDistance->getValueAt(i,j,k) ) / h;
						distToZero = - inputDistance->getValueAt(i,j,k) / i_slope;
						if ( abs(outputDistance->getValueAt(i,j,k) ) > abs(distToZero))
							outputDistance->setValueAt(i,j,k,-distToZero * utils::sgn(i_slope));
					}
					candidate = outputDistance->getValueAt(i,j,k) + (utils::sgn(inputDistance->getValueAt(i,j+1),k) * h);
					if (abs(candidate) < abs(outputDistance->getValueAt(i,j+1,k)))
						outputDistance->setValueAt(i,j+1,k, candidate);
				}
				else {
					candidate = outputDistance->getValueAt(i,j,k)  + (utils::sgn( outputDistance->getValueAt(i,j+1,k) ) * h);
					if (abs(candidate) < abs(outputDistance->getValueAt(i, j+1,k)))
						outputDistance->setValueAt(i,j+1,k, candidate);
				}
			}
		}


		for (int i = intersec_ymin; i < outputDistance->getMaxY(); i++){
		  for (int j = intersec_xmax-1; j >  outputDistance->getMinX(); j--) {
		// x-direction outputDistanceward
				//check for sign change
				if(j > intersec_xmin && i < intersec_ymax){
						// calculate new distance candidate and assign if appropriate
						candidate = outputDistance->getValueAt(i,j,k)  + (utils::sgn( inputDistance->getValueAt(i,j-1,k) ) * h);
						if (abs(candidate) < abs(outputDistance->getValueAt(i,j-1,k)))
							outputDistance->setValueAt(i,j-1,k, candidate);
				}
				else {
					candidate = outputDistance->getValueAt(i,j,k)  + utils::sgn( outputDistance->getValueAt(i,j-1,k))*h;
					if (abs(candidate) < abs(outputDistance->getValueAt(i,j-1,k)))
						outputDistance->setValueAt(i,j-1,k, candidate);
				}
			}
		}

		// y-direction forward
		for (int j = intersec_xmin; j < outputDistance->getMaxX(); j++) {
			for (int i = intersec_ymin; i < outputDistance->getMaxY()-1; i++) {
				if(j < intersec_xmax && i < intersec_ymax-1){
					if (inputDistance->getValueAt(i,j,k) * inputDistance->getValueAt(i+1,j,k) <= 0.0) {
						// interpolate
						i_slope  = (inputDistance->getValueAt(i+1,j,k) - inputDistance->getValueAt(i,j,k) )/ h;
						distToZero = - inputDistance->getValueAt(i,j,k) / i_slope;
						if ( abs(outputDistance->getValueAt(i,j,k) ) > abs(distToZero))
							outputDistance->setValueAt(i,j,k, -distToZero * utils::sgn(i_slope));
					}
					// calculate new distance candidate and assign if appropriate
					candidate = outputDistance->getValueAt(i,j,k)  + (utils::sgn( inputDistance->getValueAt(i+1,j,k) ) * h);
					if (abs(candidate) < abs(outputDistance->getValueAt(i+1,j,k) ))
						outputDistance->setValueAt(i+1,j,k, candidate);
				}
				else {
					candidate = outputDistance->getValueAt(i,j,k)  + (utils::sgn( outputDistance->getValueAt(i+1,j,k)) * h);
					if (abs(candidate) < abs(outputDistance->getValueAt(i+1,j,k)))
						outputDistance->setValueAt(i+1,j,k, candidate);
				}
			}
		}

		for (int j = intersec_xmin; j < outputDistance->getMaxX(); j++) {
			for (int i = intersec_ymax-1; i > outputDistance->getMinY(); i--) {
				if(j < intersec_xmax && i > intersec_ymin){
					// calculate new distance candidate and assign if appropriate
					candidate = outputDistance->getValueAt(i,j,k)  + (utils::sgn( inputDistance->getValueAt(i-1, j,k) ) * h);
					if (abs(candidate) < abs(outputDistance->getValueAt(i-1, j,k) ))
						outputDistance->setValueAt(i-1, j,k, candidate);
				}
				else {
					candidate = outputDistance->getValueAt(i,j,k)  + (utils::sgn( outputDistance->getValueAt(i-1, j,k) ) * h);
					if (abs(candidate) < abs(outputDistance->getValueAt(i-1,j,k)))
						outputDistance->setValueAt(i-1, j,k, candidate);
				}
			}
		}
	}
	//TODO redist into the third direction
	// for all layers the redist is done - compare into the depth to do
	outputDistance->clampValues(-handler->delta, handler->delta);

	inputDistance->resize(outputDistance->getMinX(), outputDistance->getMinY(), outputDistance->getMaxX(), outputDistance->getMaxY(), outputDistance->getMinZ(), outputDistance->getMaxZ());
	// 	 set the references for the convolution step
//	 plot_box(true,2,"Redist",true);
}

/**************************************/
// end of redist
/**************************************/


void LSbox::convolution(ExpandingVector<char>& mem_pool)
{
	double h = handler->get_h();
	if(get_status() != true ) return;
	//  set references for the convolution step

	double* ST = handler->ST;
	int n = outputDistance->getMaxX()-outputDistance->getMinX();
	int m = outputDistance->getMaxY()-outputDistance->getMinY();
	int l = outputDistance->getMaxZ()-outputDistance->getMinZ();
	int dt 	= handler->get_dt();


	//fftw_complex *fftTemp;
	int desired_size = n*m*(floor(l/2)+1)*sizeof(fftwp_complex);
	mem_pool.expand(desired_size);

	fftwp_plan fwdPlan, bwdPlan;
	
	fftwp_complex *fftTemp = (fftwp_complex*) &mem_pool[0];

		
#pragma omp critical
{
	makeFFTPlans(inputDistance->getRawData(),outputDistance->getRawData(), fftTemp, &fwdPlan, &bwdPlan);
}
	conv_generator(fftTemp,fwdPlan,bwdPlan);
//	plot_box(true,2,"Convoluted_ISO_",true);
#pragma omp critical
{
	destroyFFTWs(fwdPlan, bwdPlan);
}
	/*********************************************************************************/
	// Velocity Corrector Step: 
	/*********************************************************************************/
	// hier soll energycorrection gerechnet werden.
	// in der domainCl steht die urspr�nglich distanzfunktion, in dem arry die gefaltete
//TEST CODE 
//	if(handler->loop > 0)
//		constructBoundarySectors(/*handler->loop % Settings::AnalysisTimestep == 0*/ false);
//TEST CODE
//	if(!Settings::IsIsotropicNetwork && handler->loop!=0){
//	    vector<LSbox*>::iterator it;
//	    int intersec_xmin, intersec_xmax, intersec_ymin, intersec_ymax;
//		double weight, gamma;
//		double val;
//		double dist2OrderNeigh;
//		int nActiveGrains;
//		vector<LSbox*> IDs;
//		vector<LSbox*> IDsActive;
//
//	    if (xminId < outputDistance->getMinX())
//		  intersec_xmin = outputDistance->getMinX();
//	    else  intersec_xmin = xminId;
//
//	    if (yminId < outputDistance->getMinY())
//		  intersec_ymin = outputDistance->getMinY();
//	    else  intersec_ymin = yminId;
//
//	    if (xmaxId > outputDistance->getMaxX())
//		  intersec_xmax= outputDistance->getMaxX();
//	    else  intersec_xmax = xmaxId;
//
//	    if (ymaxId > outputDistance->getMaxY())
//		  intersec_ymax= outputDistance->getMaxY();
//	    else  intersec_ymax = ymaxId;
//
//
//	    for (int i = intersec_ymin; i < intersec_ymax; i++){
//			for (int j = intersec_xmin; j < intersec_xmax; j++) {
//				val = inputDistance->getValueAt(i,j);
//				if(val > -handler->delta){
//					if(IDLocal.getValueAt(i,j).total_elements >= 2){
//						if( i==32 && j==32)
//							double c=0;
//						if(IDLocal.getValueAt(i,j).getElementAt(1)->get_status() == true &&
//						   IDLocal.getValueAt(i,j).getElementAt(1)->inputDistance->isPointInside(i,j))
//						{
//							dist2OrderNeigh = IDLocal.getValueAt(i,j).getElementAt(1)->inputDistance->getValueAt(i,j);
//							weight = local_weights->loadWeights(IDLocal.getValueAt(i,j).local_chunks, IDLocal.getValueAt(i,j).total_elements, this, handler->ST);
//							gamma = getGBEnergyTimesGBMobility(i,j);
//							if( dist2OrderNeigh > -handler->delta) {
//								weight = -(dist2OrderNeigh/ handler->delta * (gamma - weight) )+ weight;
//								weight = weight *Settings::TriplePointDrag;
//							}
//							else weight = gamma;
//
//
//							// the weight is a function of the distance to the 2 order neighbor
//							outputDistance->setValueAt(i,j, val + (outputDistance->getValueAt(i,j) - val) * weight );b
//						}
//					}
//					else if(IDLocal.getValueAt(i,j).total_elements == 1){
//
//						gamma = getGBEnergyTimesGBMobility(i,j);
//						outputDistance->setValueAt(i,j, val + (outputDistance->getValueAt(i,j) - val) * gamma );
//					}
//
////					else if (IDLocal.getValueAt(i,j).total_elements > 2){
////						nActiveGrains = IDLocal.getValueAt(i,j).total_elements;
////						IDs.clear();
////						for (int ii = 0;ii < nActiveGrains; ii++){
////							if(isNeighbour(IDLocal.getValueAt(i,j).getElementAt(ii))) {
////								IDs.push_back(IDLocal.getValueAt(i,j).getElementAt(ii));
////							}
////						}
////						if (IDs.size()==2) weight=local_weights->loadWeights(IDs, this,handler->ST);
////						else if (IDs.size()==3){
////							weight =0;
////							IDsActive.clear();
////							IDsActive.push_back(IDs[0]);IDsActive.push_back(IDs[1]);
////							weight += local_weights->isTriplePoint(IDsActive);
////							IDsActive.clear();
////							IDsActive.push_back(IDs[1]);IDsActive.push_back(IDs[2]);
////							weight += local_weights->isTriplePoint(IDsActive);
////							IDsActive.clear();
////							IDsActive.push_back(IDs[0]);IDsActive.push_back(IDs[2]);
////							weight += local_weights->isTriplePoint(IDsActive);
////							weight /= 3;
////						}
////						else weight = 0.5*handler-> hagb;
////
////						outputDistance->setValueAt(i,j, val + (outputDistance->getValueAt(i,j) - val) * weight);
////					}
//				}
//				// this should avoid spikes, depending on thrird order neighbour interaction, occuring by periodicity of the the convoluted function
//				else  outputDistance->setValueAt(i,j, - handler->delta);
//
////				double bulkenergy = 50;
//				outputDistance->setValueAt(i,j,outputDistance->getValueAt(i,j)/*+(dt/2*bulkenergy)*/);
//			}
//		}
//	}


	IDLocal.clear();
		
// 	IDLocal.clear(); 
	get_new_IDLocalSize();
	IDLocal.resize(xminId, yminId, xmaxId, ymaxId, zminId, zmaxId);
// 	if(id == 15 && handler->loop >90)plot_box(true,2,"Convoluted_2_");
// 	plot_box(true,1,"Convoluted_",true);
//	plot_box(true,2,"Convoluted_",true);


}



//
//LSbox::LSbox(int id, int nedges, double* edges, double phi1, double PHI, double phi2, grainhdl* owner) : id(id), nvertices(nedges), handler(owner){
//	exist = true;
//
//	quaternion = new double[4];
//	double euler[3] = {phi1,PHI,phi2};
//	(*(handler->mymath)).euler2quaternion( euler, quaternion );
//
//	int grid_blowup = owner->get_grid_blowup();
//	double h = owner->get_h();
//    // determine size of grain
//    int xmax = 0;
//	int xmin = handler->get_ngridpoints();
//	int ymax = 0;
//	int ymin = xmin;
//
//	vektor x1(2), x2(2);
//	exist = true;
//
//	for (unsigned int k=0; k < nedges; k++){
//	  x1[0]=edges[(4*k)+1]; x1[1]=edges[4*k];
//	  x2[0]=edges[(4*k)+3]; x2[1]=edges[(4*k)+2];
//
//		//	for convention:
//		//	x[i][j]:
//		//	i = Zeilenindex(y-direction)
//		// 	j = Spaltenindex(x-direction)
//
//		// check for "Zeilen" Minima/Maxima
//		if (x1[0]/h < ymin) ymin = x1[0]/h;
//		if (x2[0]/h < ymin) ymin = x2[0]/h;
//
//		if (x1[0]/h > ymax) ymax = (x1[0]/h);
//		if (x2[0]/h > ymax) ymax = (x2[0]/h);
//
//		// check for "Spalten" Minima/Maxima
//		if (x1[1]/h < xmin) xmin = x1[1]/h;
//		if (x2[1]/h < xmin) xmin = x2[1]/h;
//
//		if (x1[1]/h > xmax) xmax = (x1[1]/h);
//		if (x2[1]/h > xmax) xmax = (x2[1]/h);
//    }
//	xmax += 2*grid_blowup;
//	ymax += 2*grid_blowup;
//
//
//	inputDistance = new DimensionalBufferReal(xmin, ymin, xmax, ymax);
//	outputDistance = new DimensionalBufferReal(xmin, ymin, xmax, ymax);
//
// 	inputDistance->resizeToSquare(handler->get_ngridpoints());
// 	outputDistance->resizeToSquare(handler->get_ngridpoints());
//	inputDistance->clearValues(0.0);
//	outputDistance->clearValues(0.0);
//
//	get_new_IDLocalSize();
//	IDLocal.resize(xminId, yminId, xmaxId, ymaxId);
//
//// 	cout << "made a new box: xmin="<<xmin<< " xmax="<<xmax <<" ymin="<<ymin << " ymax="<<ymax<<endl;
//	local_weights=new Weightmap(owner);
//	boundaryGrain = false;
//}



//void LSbox::distancefunctionToEdges(int nedges, double* edges){
//// 	plot_box(false);
//	int grid_blowup = handler->get_grid_blowup();
//	double h = handler->get_h();
//	int i,j,k;
//	double d, dmin,lambda;
//	vektor u(2), a(2), p(2), x1(2), x2(2);
//
//	for (i=outputDistance->getMinY();i<outputDistance->getMaxY();i++){ // ¸ber gitter iterieren
//	  for (j=outputDistance->getMinX();j<outputDistance->getMaxX();j++){
//            dmin = 1000.;
//            p[0] = (i-grid_blowup)*h; p[1] = (j-grid_blowup)*h;
//
//            for(int k=0; k < nedges; k++) {
//				x1[0]=edges[(4*k)+1]; x1[1]=edges[4*k];
//				x2[0]=edges[(4*k)+3]; x2[1]=edges[(4*k)+2];
//				if (x1 != x2){
//					a = x1;
//					u = x2-x1;
//					lambda=((p-a)*u)/(u*u);
//
//					if(lambda <= 0.) 				d = (p-x1).laenge();
//					if((0. < lambda) && (lambda < 1.)) 		d = (p-(a+(u*lambda))).laenge();
//					if(lambda >= 1.) 				d = (p-x2).laenge();
//					d= abs(d);
//					if(abs(d)< abs(dmin)) dmin=d;
//				}
//            }
//			if (abs(dmin) < handler->delta){
//				outputDistance->setValueAt(i, j, dmin);
//			}
//			else{
//				outputDistance->setValueAt(i, j, handler->delta * utils::sgn(dmin));
//			}
//        }
//	}
//
//	int count = 0;
//    for (j=outputDistance->getMinX();j<outputDistance->getMaxX();j++){ // ¸ber gitter iterieren
//	    i=outputDistance->getMinY();
//	    count = 0;
//	    while( i<outputDistance->getMaxY() && count < 1) {
//	    	outputDistance->setValueAt(i, j, -abs(outputDistance->getValueAt(i,j)));
//		    if ( -outputDistance->getValueAt(i,j) <=  h )
//		    	count++;
//		    i++;
//	    }
//	    i=outputDistance->getMaxY()-1;
//	    count =0;
//	    while( i>=outputDistance->getMinY() && count < 1) {
//	    	outputDistance->setValueAt(i, j, -abs(outputDistance->getValueAt(i,j)));
//		    if ( -outputDistance->getValueAt(i,j) <= h )
//		    	count++;
//		    i--;
//	    }
//    }
//
//    for (i=outputDistance->getMinY();i<outputDistance->getMaxY();i++){
//	    j=outputDistance->getMinX();
//	    count = 0;
//	    while( j<outputDistance->getMaxX() && count < 1 ) {
//	    	outputDistance->setValueAt(i, j, -abs(outputDistance->getValueAt(i,j)));
//		    if ( -outputDistance->getValueAt(i,j) <=  h )
//		    	count++;
//		    j++;
//	    }
//	    j=outputDistance->getMaxX()-1;
//	    count =0;
//	    while( j>=outputDistance->getMinX() && count < 1  ) {
//	    	outputDistance->setValueAt(i, j, -abs(outputDistance->getValueAt(i,j)));
//		    if ( -outputDistance->getValueAt(i,j)<=  h )
//		    	count++;
//		    j--;
//	    }
//    }
//}






//double LSbox::getGBEnergyTimesGBMobility(int i,int j, int k){
//	LSbox* neighbour = IDLocal.getValueAt(i,j,k).getElementAt(0);
//	vector<characteristics>::iterator it;
//	for (it = grainCharacteristics.begin(); it != grainCharacteristics.end(); it++){
//		if (neighbour == (*it).directNeighbour){
//			return (*it).energyDensity*(*it).mobility ;
//		}
//	}
//	neighbour = IDLocal.getValueAt(i,j).getElementAt(1);
//	for (it = grainCharacteristics.begin(); it != grainCharacteristics.end(); it++){
//		if (neighbour == (*it).directNeighbour){
//			return (*it).energyDensity*(*it).mobility ;
//		}
//	}
//	return 1.0;
//
//}
//
//double LSbox::getGBEnergyTimesGBMobility(LSbox* neighbour){
//	vector<characteristics>::iterator it;
//	for (it = grainCharacteristics.begin(); it != grainCharacteristics.end(); it++){
//		if (neighbour == (*it).directNeighbour){
//			return (*it).energyDensity*(*it).mobility ;
//		}
//	}
//	return 1.;
//}
//
//double LSbox::getGBEnergy(LSbox* neighbour){
//	vector<characteristics>::iterator it;
//	for (it = grainCharacteristics.begin(); it != grainCharacteristics.end(); it++){
//		if (neighbour == (*it).directNeighbour){
//			return (*it).energyDensity;
//		}
//	}
//}


/**************************************/
/**************************************/

/*
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

  
}*/

// >>>>>>> review_single

/**************************************/
/**************************************/



/*
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
// 	plot_box(false,1,"no");
}*/


/**************************************/
// end of Comparison
/**************************************/



// Find Contour operates on inputDistance
/**************************************/
/**************************************/

//void LSbox::find_contour() {
//	exist = false;
//
//	contourGrain.clear();
////	vector<GrainJunction> Junctions;
//    MarchingSquaresAlgorithm marcher(*inputDistance, IDLocal, this);
//    junctions.clear();
//    exist = marcher.generateContour(contourGrain, junctions);
//
//	if(!exist) return;
//    int grid_blowup = handler->get_grid_blowup();
//	int m = handler->get_ngridpoints();
//
//    int xminNew = m, xmaxNew = 0, yminNew = m, ymaxNew = 0;
//    for(int i=0; i<contourGrain.size(); i++)
//    {
//    	if(int(contourGrain[i].x + 0.5) - grid_blowup < xminNew)
//    		xminNew = int(contourGrain[i].x + 0.5) - grid_blowup;
//    	if(int(contourGrain[i].x + 0.5) + grid_blowup > xmaxNew)
//    		xmaxNew = int(contourGrain[i].x + 0.5) + grid_blowup;
//    	if(int(contourGrain[i].y + 0.5) - grid_blowup < yminNew)
//    		yminNew = int(contourGrain[i].y + 0.5) - grid_blowup;
//    	if(int(contourGrain[i].y + 0.5) + grid_blowup > ymaxNew)
//    		ymaxNew = int(contourGrain[i].y + 0.5) + grid_blowup;
//    }
//
//
//	double h = handler->get_h();
//	int loop = handler->loop;
//
//	if(xminNew < 0 || yminNew < 0 || ymaxNew > m|| xmaxNew > m) {
//		cout <<endl << "Timestep: " <<handler->loop << endl << endl;
//
//		cout << "WARNING - undefined Boxsize in Box: "<< id <<" in Timestep: "<<loop << "!!" <<endl;
//		cout << "Number of gridpoints: " << m << endl;
//		cout << yminNew << " || " << xminNew << " || " << ymaxNew  << " || " << xmaxNew << endl;
//		for(int i=0; i<contourGrain.size(); i++){
//			cout << contourGrain[i].y << "   " << contourGrain[i].x << endl;
//		}
//	}
//
//    // compute Volume and Energy
//	if ( (loop % int(Settings::AnalysisTimestep)) == 0 || loop == Settings::NumberOfTimesteps ) {
//		computeVolumeAndEnergy();
//
//	}
//	else updateFirstOrderNeigbors();
//
//	if(grainCharacteristics.size() < 2 ) {
////		cout << endl << "Timestep: " <<handler->loop << endl;
////		cout << "GRAIN: " << id << " has a positive Volume but less than 2 neighbors" << endl;
////		plot_box(true, 2, "error_grain", true);
////		plot_box(true, 1, "error_grain", true);
////		plot_box_contour(handler->loop, true);
//		exist =false;
//	}
//	outputDistance->resize(xminNew, yminNew, xmaxNew, ymaxNew);
// 	outputDistance->resizeToSquare(handler->get_ngridpoints());
//
//	return;
//}

//revise algorithm to update first order neighbors
//void LSbox::updateFirstOrderNeigbors(){
//	grainCharacteristics.clear();
//	vector<characteristics>::iterator it;
//	double h = handler->get_h();
//	double line_length;
//	double thetaMis=0, mu;
//	double theta_ref = 15.0 * PI / 180.0;
//	double gamma_hagb = handler->hagb;
//	int i;
//	for(i=0; i<contourGrain.size() - 1; i++){
//	  double px =contourGrain[i].x;
//	  double py =contourGrain[i].y;
//	  int pxGrid = int(px);
//	  int pyGrid = int(py);
//	  if(inputDistance->getValueAt(pyGrid,pxGrid) > 0) {
//			pxGrid = int(px+1);
//			if(inputDistance->getValueAt(pyGrid,pxGrid) > 0){
//				pyGrid = int(py +1);
//				if(inputDistance->getValueAt(pyGrid,pxGrid) > 0){
//					pxGrid = int(px);
//				}
//			}
//		}
//
//	  for (it = grainCharacteristics.begin(); it != grainCharacteristics.end(); it++){
//		 if (IDLocal.getValueAt(py,px).getElementAt(0) == (*it).directNeighbour)
//			  break;
//	  }
//	  if (it == grainCharacteristics.end()){
//		  double energyLocal;
//		  if(Settings::IsIsotropicNetwork){
//			  energy = 1.0;
//		  }
//		  else{
//			  thetaMis = mis_ori( IDLocal.getValueAt(pyGrid, pxGrid).getElementAt(0));
//			  if (thetaMis <= theta_ref)
//				  energyLocal = gamma_hagb * ( thetaMis / theta_ref) * (1.0 - log( thetaMis / theta_ref));
//			  else
//				  energyLocal = gamma_hagb;
//		  }
//		  if(Settings::UseMobilityFactor && (!Settings::IsIsotropicNetwork)) mu = GBmobilityModel(thetaMis);
//		  else mu = 1;
//		  grainCharacteristics.push_back(characteristics( IDLocal.getValueAt(pyGrid, pxGrid).getElementAt(0), 0, energyLocal,thetaMis, mu));
//		  it = grainCharacteristics.end();
//		  it--;
//	  }
//	  line_length = sqrt(((contourGrain[i].x-contourGrain[i+1].x)*(contourGrain[i].x-contourGrain[i+1].x)) + ((contourGrain[i].y-contourGrain[i+1].y)*(contourGrain[i].y-contourGrain[i+1].y)));
//	  //save length in GrainCharaczeristics
//	  it->length += (line_length*h);
//	}
//}

//void LSbox::computeVolumeAndEnergy()
//{
//	double dA = volume;
//	volume = 0;
//	perimeter = 0;
//	energy = 0;
//	double line_length, mu;
//	double h = handler->get_h();
//	double thetaMis=0;
//	double theta_ref = 15.0 * PI / 180.0;
//	double gamma_hagb = handler->hagb;
//	double dt = handler->get_dt();
//	int i;
//	grainCharacteristics.clear();
//	vector<characteristics>::iterator it;
//
//	for(i=0; i<contourGrain.size() - 1; i++){
//
//		// Gaussian Trapez Formula:
//		volume += (contourGrain[i].y+contourGrain[i+1].y)*(contourGrain[i].x-contourGrain[i+1].x);
//		double px =contourGrain[i].x;
//		double py =contourGrain[i].y;
//		int pxGrid = int(px);
//		int pyGrid = int(py);
//		// evaluate the ID at a outer point -> there must be one, we are at the boundray
//		if(inputDistance->getValueAt(pyGrid,pxGrid) > 0){
//			pxGrid = int(px+1);
//			if(inputDistance->getValueAt(pyGrid,pxGrid) > 0){
//				pyGrid = int(py +1);
//				if(inputDistance->getValueAt(pyGrid,pxGrid) > 0){
//					pxGrid = int(px);
//				}
//			}
//		}
//		if(Settings::IsIsotropicNetwork){
//		  contourGrain[i].energy = 1.0;
//		}
//		else if(Settings::ResearchMode ==1 && Settings::MicrostructureGenMode != E_GENERATE_TESTCASE){
//			theta_ref = 42 * PI /180;
//			thetaMis = mis_ori( IDLocal.getValueAt(pyGrid, pxGrid).getElementAt(0));
//			if (thetaMis <= theta_ref)	contourGrain[i].energy = 0.3;
//			else contourGrain[i].energy = gamma_hagb;
//		}
//		//! Handles the colouring for the E_GENERATE_TESTCASE
//		else if(Settings::ResearchMode == 1 && Settings::MicrostructureGenMode == E_GENERATE_TESTCASE){
//
//			//!cout << "the current grain: " << get_id() << " the grain in vicinity: " << IDLocal.getValueAt(pyGrid, pxGrid).getElementAt(0)->get_id() << endl;
//			contourGrain[i].energy = handler->weightsMatrix[get_id()][IDLocal.getValueAt(pyGrid, pxGrid).getElementAt(0)->get_id()];
//
//		}
//		else{
//			thetaMis = mis_ori( IDLocal.getValueAt(pyGrid, pxGrid).getElementAt(0));
//		  if (thetaMis <= theta_ref)
//		    contourGrain[i].energy = gamma_hagb * ( thetaMis / theta_ref) * (1.0 - log( thetaMis / theta_ref));
//		  else
//		    contourGrain[i].energy = gamma_hagb;
//		}
//
//		// 	Check if the direct neigbourGrain is already in the vector
//		for (it = grainCharacteristics.begin(); it != grainCharacteristics.end(); it++){
//			if (IDLocal.getValueAt(pyGrid, pxGrid).getElementAt(0) == (*it).directNeighbour)
//				break;
//		}
//		if (it == grainCharacteristics.end()){
//			if(Settings::UseMobilityFactor && (!Settings::IsIsotropicNetwork)) mu = GBmobilityModel(thetaMis);
//			else mu = 1;
//			grainCharacteristics.push_back(characteristics(IDLocal.getValueAt(pyGrid, pxGrid).getElementAt(0), 0, contourGrain[i].energy,thetaMis, mu));
//			it = grainCharacteristics.end();
//			it--;
//		}
//
//		line_length = sqrt(((contourGrain[i].x-contourGrain[i+1].x)*(contourGrain[i].x-contourGrain[i+1].x)) + ((contourGrain[i].y-contourGrain[i+1].y)*(contourGrain[i].y-contourGrain[i+1].y)));
//
//		//save length in GrainCharaczeristics
//		line_length *=h;
//		it->length += line_length;
//
//		// save Grain properties:
//		perimeter += line_length;
//		energy += (contourGrain[i].energy * line_length);
//
//	}
//
//	volume = abs(volume) *h*h;
//	contourGrain[contourGrain.size()-1].energy = contourGrain[0].energy;
//
//	//!
//	//! Evaluating the area variation in the current time step and
//	//! saving this variation together with the current number
//	//! of neighbours in a vector for further analyses.
//	//! The area variation is normalized by a factor coming from the
//	//! Neumann-Mullins equation.
//	//!
//
//	dA = volume - dA;
//
//	dA /= Settings::AnalysisTimestep*handler->get_dt();
//	dA *= (3/PI);
//	VolEvo.push_back(VolEvolution(dA,grainCharacteristics.size()));
//
//// 	for (it = grainCharacteristics.begin(); it != grainCharacteristics.end(); it++){
//// 		 printf("%d\t %lf\t %lf\t%lf\n", (*it).directNeighbour->get_id(),(*it).length,(*it).energyDensity,(*it).mis_ori);
//// 	}
//}





/**************************************/
// plot the box and all its properties
/**************************************/
//void LSbox::resizeGrid(double shrinkFactor){
//
//  int realDomainSizen = handler->get_realDomainSize()*(1-shrinkFactor)+1;
////   int ngridpointsn = realDomainSizen+2*handler->get_grid_blowup();
//  double h = handler->get_h();
//  double hn = 1.0/(realDomainSizen);
//
//
//
//  int minXnew = outputDistance->getMinX()*(h/hn);
//  int maxXnew = outputDistance->getMaxX()*(h/hn)+1;
//  int minYnew = outputDistance->getMinY()*(h/hn);
//  int maxYnew = outputDistance->getMaxY()*(h/hn)+1;
//  double xl,xr,yo,yu;
//
//  double pointx, pointy;
//
//
//// resize to complete superposition
//      if ( minXnew * hn > outputDistance->getMinX()*h)
//	minXnew--;
//      if ( minYnew * hn > outputDistance->getMinY()*h)
//	minYnew--;
//      if ( maxXnew * hn < outputDistance->getMaxX()*h)
//	maxXnew++;
//      if ( maxYnew * hn < outputDistance->getMaxY()*h)
//	maxYnew++;
//
////   plot_box(true, 2, "before_resize", true);
//
//  inputDistance->resize(minXnew, minYnew, maxXnew, maxYnew);
//  for (int i = minYnew; i < maxYnew; i++){
//    for (int j = minXnew; j < maxXnew; j++){
//	pointx = j*(hn/h);
//	pointy = i*(hn/h);
//
//	xl= int(pointx);
//	xr= int(pointx+1);
//	yo= int(pointy+1);
//	yu= int(pointy);
//
//
//	if (xr > outputDistance->getMaxX()-2||yo > outputDistance->getMaxY()-2||yu < outputDistance->getMinY()||xl < outputDistance->getMinX()){
//	  inputDistance->setValueAt(i,j, -handler->delta);
//	  continue;
//	}
//	double ro,ru,newDistVal;
//	ro = 1/(xr-xl)*((xr-pointx)*outputDistance->getValueAt(yo, xl)+(pointx-xl)*outputDistance->getValueAt(yo, xr));
//	ru = 1/(xr-xl)*((xr-pointx)*outputDistance->getValueAt(yu, xl)+(pointx-xl)*outputDistance->getValueAt(yu, xr));
//	newDistVal = 1/(yo-yu)*((yo-pointy)*ru+(pointy-yu)*ro);
//	if (newDistVal != newDistVal) {
//	  char waitbuffer;
//	  cerr << " nan " << endl;
//	  cin >> waitbuffer;
//	}
//	inputDistance->setValueAt(i,j, newDistVal);
//
//    }
//  }
//  outputDistance->resize(minXnew, minYnew, maxXnew, maxYnew);
//
//
////   plot_box(true, 1, "after_resize", true);
//    //plot_box for all boxes and compare with prior !
//
//
//}


//void LSbox::plot_box_contour(int timestep, bool plot_energy, ofstream* dest_file)
//{
//	if(exist == false)
//		return;
//    ofstream* output_file = dest_file;
//    if(dest_file == NULL)
//    {
//    	output_file = new ofstream();
//        stringstream filename;
//        filename<<"Contourline_"<< id;
//        filename<<"_Timestep_"<<timestep;
//        filename<<".gnu";
//        output_file->open(filename.str());
//    }
//	ofstream& file = *output_file;
//    if ( plot_energy)
//    {
//		for(const auto& iterator : contourGrain)
//		{
////			file << (iterator.x-handler->get_grid_blowup()) *handler->get_h()  << "\t" << (iterator.y-handler->get_grid_blowup()) *handler->get_h()  << endl;
//			file << iterator.x << "\t" << iterator.y<< "\t" << iterator.energy << endl;
//		}
//    }
//    else {
//		for(const auto& iterator : contourGrain)
//		{
//			file << (iterator.x-handler->get_grid_blowup()) *handler->get_h()  << "\t" << (iterator.y-handler->get_grid_blowup()) *handler->get_h()  << endl;
//		}
//    }
//    file<<endl;
//    if(dest_file == NULL)
//    {
//    	file.close();
//    	delete output_file;
//    }
//}


void LSbox::plot_box(bool distanceplot, int select, string simstep, bool local){
  
	cout <<" \nGrain  Info: " << endl;
	cout << " ID :" <<id << endl;
	cout << " xminIn, xmaxIn, yminIn, ymaxIn :" << inputDistance->getMinX() << " || "<< inputDistance->getMaxX() << " || " << inputDistance->getMinY() << " || " << inputDistance->getMaxY() << endl;
	cout << " xminOut, xmaxOut, yminOut, ymaxOut :" << outputDistance->getMinX() << " || "<< outputDistance->getMaxX() << " || " << outputDistance->getMinY() << " || " << outputDistance->getMaxY() << endl;
	cout << " xminId, xmaxId, yminId, ymaxId :" << xminId << " || "<< xmaxId << " || " << yminId << " || " << ymaxId << endl;
//     if (distanceplot==true) utils::print_2dim_array(distance,ymax-ymin,xmax-xmin);
// 		else cout << " no distance values in storage!" << endl;
	cout << " quaternion: " << quaternion[0] << " || "<< quaternion[1]<< " || " <<quaternion[2] << " || " << quaternion[3]<< endl;
	if(grainCharacteristics.empty()!=true){
		cout << " List of Neighbors : ";
		vector<characteristics> :: iterator it;
		for (it = grainCharacteristics.begin(); it != grainCharacteristics.end(); it++){
			cout << (*it).directNeighbour->getID() << " || ";
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
		if(select == 2 && !local) {
		filename<< "BoxDistance_"<< simstep << "out_T" << loop << "_" << id << ".gnu";
		datei.open(filename.str());
			for (int i = 0; i < handler->get_ngridpoints(); i++){
				for (int j = 0; j < handler->get_ngridpoints(); j++){
					if( i >= outputDistance->getMinY() && i < outputDistance->getMaxY() && j >=outputDistance->getMinX() && j < outputDistance->getMaxX()) {
						datei << ::std::fixed << outputDistance->getValueAt(i,j) << "\t";
					}
					else datei << ::std::fixed << -handler->delta<< "\t";
				}
			datei << endl;
			}	
		}		
		
		if(select == 2 && local) {
		filename<< "BoxDistance_"<< simstep << "out_T" << loop << "_" << id << ".gnu";
		datei.open(filename.str());
			for (int i = outputDistance->getMinY(); i < outputDistance->getMaxY(); i++){
				for (int j = outputDistance->getMinX(); j < outputDistance->getMaxX(); j++){
					if( i >= outputDistance->getMinY() && i < outputDistance->getMaxY() && j >=outputDistance->getMinX() && j < outputDistance->getMaxX()) {
						datei << ::std::fixed << outputDistance->getValueAt(i,j) << "\t";
					}
					else datei << ::std::fixed << -handler->delta<< "\t";
				}
			datei << endl;
			}	
		}
		if(select == 1 && local) {
		filename<< "BoxDistance_"<< simstep << "in_T" << loop << "_" << id << ".gnu";
		datei.open(filename.str());
			for (int i = inputDistance->getMinY(); i < inputDistance->getMaxY(); i++){
				for (int j = inputDistance->getMinX(); j < inputDistance->getMaxX(); j++){
					if( i >= inputDistance->getMinY() && i < inputDistance->getMaxY() && j >=inputDistance->getMinX() && j < inputDistance->getMaxX()) {
						datei << ::std::fixed << inputDistance->getValueAt(i,j)<< "\t";
					}
					else datei << ::std::fixed << -handler->delta<< "\t";
				}
			datei << endl;
			}	
		}		
		if(select == 1 && !local) {
		filename<< "BoxDistance_"<< simstep << "in_T" << loop << "_" << id << ".gnu";
		datei.open(filename.str());
			for (int i = 0; i < handler->get_ngridpoints(); i++){
				for (int j = 0; j < handler->get_ngridpoints(); j++){
					if( i >= inputDistance->getMinY() && i < inputDistance->getMaxY() && j >=inputDistance->getMinX() && j < inputDistance->getMaxX()) {
						datei << ::std::fixed << inputDistance->getValueAt(i,j)<< "\t";
					}
					else datei << ::std::fixed << -handler->delta<< "\t";
				}
			datei << endl;
			}	
		}		
		datei.close();
    }
    if(!contourGrain.empty()){
		for(int i=0; i<contourGrain.size() - 1; i++){
			double px =contourGrain[i].x;
			double py =contourGrain[i].y;
			cout << py << "   "<< px << endl;
		}
    }
}



//double LSbox::GBmobilityModel(double thetaMis){
//	double theta_ref = 15.0 * PI / 180.0;
//	double theta_ref_2 = 35.0 * PI / 180.0;
//	double tubeWidth= 10.0 * PI / 180.0;
//	double lagbM = 0.1;
//	double hagbM = 0.5;
//	double mu;
//	if(thetaMis < theta_ref) mu = 0.1;
//	else mu = hagbM;
//// 	else if(thetaMis > theta_ref && thetaMis < theta_ref_2) mu = hagbM;
//// 	else if(thetaMis > (theta_ref_2+tubeWidth)) mu = hagbM;
//// 	else  mu = hagbM + (0.5 *sin((2*PI*thetaMis/ (2*tubeWidth ))));
//// 	cout << thetaMis << "  " << theta_ref << "  "<< theta_ref_2 << "  "<< mu <<endl;
//// 	char buf;
//// 	cin >> buf;
//	return mu;
//}





//void LSbox::saveGrain(ofstream* destfile ){
//	 ofstream& file = *destfile;
//	 file << id << "\t" << contourGrain.size() << "\t" << quaternion[0] << "\t" << quaternion[1]<< "\t" << quaternion[2]<< "\t" << quaternion[3] << endl;
//
//	 for(const auto& iterator : contourGrain){
//		file << iterator.x << "\t" << iterator.y<<  endl;
//	 }
//}

//int getIdxSector(vector<ContourSector>& sectors, SPoint& point)
//{
//	for(int i=0; i<sectors.size(); i++)
//		if(sectors[i].isPointWithinSectoinRadiuses(point))
//			return i;
//	return -1;
//}
//void LSbox::constructBoundarySectors(bool test_plot)
//{
//	vector<ContourSector> sectors;
//	//Step 1 = merged circles...
//	ContourSector::INNER_CIRCLE_RADIUS = 3;
//	for(int i=0; i<junctions.size(); i++)
//	{
//		bool merged = false;
//		for(int j=0; j<sectors.size() && merged == false; j++)
//		{
//			if(sectors[j].mergeWith(&junctions[i]))
//				merged = true;
//		}
//		if(merged == false)
//			sectors.push_back(ContourSector(&junctions[i]));
//	}
//	//Step 2 - identify points on the contourline
//	int realContourSize = contourGrain.size() - 1;
//	int currentSector = getIdxSector(sectors, contourGrain[0]);
//	for(int i=0; i<realContourSize; i++)
//	{
//		int P = getIdxSector(sectors, contourGrain[i]);
//		if(P == -1)
//		{
//			if(currentSector == -1)
//				continue;
//			else	//currentSector != -1
//			{
//				sectors[currentSector].setLeftContourPoint((i-1+realContourSize)%realContourSize);
//				currentSector = -1;
//				continue;
//			}
//		}
//		else	//P!=-1
//		{
//			if(currentSector == -1)
//			{
//				sectors[P].setRightContourPoint((i+realContourSize)%realContourSize);
//				currentSector = P;
//				continue;
//			}
//			else if (currentSector == P)
//			{
//				continue;
//			}
//			else	//currentSector != -1 && currentSector != P
//			{
//				sectors[currentSector].setLeftContourPoint((i-1+realContourSize)%realContourSize);
//				sectors[P].setRightContourPoint((i+realContourSize)%realContourSize);
//				currentSector = P;
//				continue;
//			}
//		}
//	}
//
//	//! test_plot = false
//	test_plot = false;
//	if(test_plot)
//	{
//		plot_box_contour(handler->loop, true);
//    	ofstream output_file;
//        stringstream filename;
//        filename<<"Sectors_"<< id;
//        filename<<"_Timestep_"<<handler->loop;
//        filename<<".gnu";
//        output_file.open(filename.str());
//        for(int i=0; i<junctions.size(); i++)
//        {
//        	output_file <<junctions[i].coordinates.x <<" "
//        	        	<<junctions[i].coordinates.y <<endl<<endl;
//        }
//        output_file<<endl;
//        for(int i=0; i<sectors.size(); i++)
//        {
//        	sectors[i].debugPrintSector(contourGrain, output_file);
//        }
//
//	}
//}
