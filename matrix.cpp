#include "matrix.h"

using namespace std;
using namespace voro;


matrix::matrix(){}

matrix::matrix(int m, int n) : m(m), n(n) {
	int id=0;
    x=new vektor*[m];
    for (int i=0;i<m;i++) x[i]=new vektor(n);
}
matrix::matrix(int m, int n, int id) : m(m), n(n), id(id) {
	x=new vektor*[m];
    for (int i=0;i<m;i++) x[i]=new vektor(n);
}

matrix::matrix(int m, int n, int id, double startval) : m(m), n(n), id(id) {
	x=new vektor*[m];
    for (int i=0;i<m;i++) x[i]=new vektor(n, startval);
}



matrix::~matrix() {
    for (int i=0;i<m;i++) delete x[i];
    delete [] x;
// 	cout << "destroyed matrix with  id: "<< id << endl;
}

matrix::matrix(const matrix& v) : m(v.m), n(v.n), id(v.id) {
    x=new vektor*[m];
    for (int i=0;i<m;i++) {
        x[i]=new vektor(n);
        *(x[i])=v[i];
    }
}

vektor& matrix::operator[](int i) {
    if (0<= i < m) return *x[i];
    else {
        outOfBoundsException e(i, i);
        cout << "invalid index " << e.what() << endl;
    }
}

const vektor& matrix::operator[](int i) const {
    if (0<= i < m) return *x[i];
    
    else {
        outOfBoundsException e(i, i);
        cout << "invalid index " << e.what() << endl;
    }
    
}

// "=" Operator
double& matrix::operator=(const matrix& A){
    if (this != &A) {
        assert(m == A.m && n == A.n);
        for (int i = 0; i < m; i++) *x[i] = *(A.x[i]);
    }
}


// "+" Operator
matrix matrix::operator+(const matrix& A){
    matrix C(m,n);
    assert(m == A.m && n == A.n);
    for (int i = 0; i < m; i++) *(C.x[i]) = (*x[i]) + *(A.x[i]);
    return C;
}



// "-" Operator
matrix matrix::operator-(matrix& A){
    matrix C(m,n);
	C.grains = this->grains;
	C.id = id;
    assert(m == A.m && n == A.n);
    for (int i = 0; i < m; i++) *(C.x[i]) = (*x[i]) - *(A.x[i]);
    return C;
}

// "*" Matrixmultiplikation
matrix matrix::operator*(const matrix& A){
    matrix C(m,n);
    assert(m == A.n && n == A.m);
    for (int i = 0; i < m; i++) { // alle EintrÃ¤ge der Ergebnismatrix 0 setzen
        for (int j = 0; j < n; j++) C[i][j] = 0;
    }
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            for (int k = 0; k < n; k++) C[i][j] +=  (*x[i])[k] * A[k][j];
        }
    }
    return C;
}

// Matrix * Vektor
vektor matrix::operator*(const vektor& v){
    vektor erg(m);
    assert(n == v.n);
    for (int i = 0; i < n; i++) erg[i] = 0; // alle EintrÃ¤ge der Ergebnismatrix 0 setzen
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++) erg[i] +=  (*x[i])[j] * v[j];
    }
    return erg;
}

// Matrix/ Matrix
matrix matrix::operator/(matrix& A){
    assert(m == A.m && n == A.n);
    matrix C(m,n);
    for (int i = 0; i < m; i++) { // alle EintrÃ¤ge der Ergebnismatrix 0 setzen
        for (int j = 0; j < n; j++) C[i][j] = 0;
    }
    vektor e(n);
    for (int i = 0; i < n; i++) e[i] = 0;
    vektor x(n);
    vektor zerg(n);
    for (int i = 1; i <= n; i++) {
        e[i] = 1;
        x = e/A;
        zerg = *this * x;
        for (int j = 0; j < n; j++) {
            C[i-1][j]=zerg[j];
        }
        e[i] = 0;
    }
    return C;
}

ostream& operator<<(ostream &os, const matrix& A) {
    int i;
	for(i=0; i< (A.m); i++){
		os << *A.x[i] << endl;
	}
    return os;
}

void matrix::save_matrix(const char* filename){
    ofstream datei;
    datei.open(filename);
    datei << *this << endl;
    datei.close();
}

double matrix::entry(const int i, const int j){
    return((*x[i])[j]);
}

void matrix::mult_with_scalar(const double d){
  for (int i = 0; i < m; i++) 
	for (int j = 0; j < n; j++) 
	  (*this)[i][j] *= d;
}
    

matrix matrix::distancefunction(voronoicell_neighbor& c, int *gridIDs, double *part_pos, int grid_blowup, double h){
    int i,j,k;
	double d, dmin,lambda;
	int m=get_m();
	int n=get_n();
	
	vektor u(2), a(2), p(2), x1(2), x2(2);
    vector<double> vv;
	c.vertices(part_pos[3*id],part_pos[3*id+1],part_pos[3*id+2],vv);
	
	double domain_vertices[] = {0.,0.,1.,0.,1.,1.,0.,1.,0.,0.}; // array of vertices to loop over
	
	for (i=0;i< m;i++){ // über gitter iterieren
        for (j=0;j< n;j++){
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
//                     (ID_mat[0][i*m +j])->getID()    
                        if(id==gridIDs[i*m +j] && ((grid_blowup < i) && (i < (m- grid_blowup))) && ((grid_blowup < j) && (j < (m- grid_blowup)))) d=abs(d);
                        else d= -abs(d);
                        
                        if(abs(d)< abs(dmin)) dmin=d;
                    }
                }
            }
            //  if (abs(dmin) < DELTA || dmin > 0.0) (*x[i])[j]= dmin;
            //else (*x[i])[j] = -DELTA;
            (*x[i])[j]= dmin;
        }
	}
	return(*this);
}



/*************************************************************************/
// BOX Functions
/*************************************************************************/

bool matrix::addBox(LSbox* aBox){
    vector<LSbox*>::iterator it;
    // check if new box interesects with any of the current boxes
    bool intersect = false;
    if (!grains.empty())
        for (it = grains.begin(); it != grains.end(); it++) {
            intersect = aBox->checkIntersect(*it);
            if (intersect) break;
        }
    if (!intersect) {
        grains.push_back(aBox);
        aBox->setDomain(this);
        return true;
    }
    return false;
}


vector<LSbox*> matrix::getBoxList() {
    return grains;
}



bool matrix::grainCheck(double h, int grid_blowup, vector<LSbox*>& buffer){
	char buffer1;
	bool exist;
	vector<LSbox*>::iterator it,it2;
    for(it = grains.begin(); it != grains.end();){
		if((**it).get_status()==false) {
			cout << "try to delete grain " << (**it).get_id() << " in domain "<< id << endl;
			// the grain has disappeared the timestep before
			//now we can clean up the memory
			int buffer =(**it).get_id();
			delete (*it); grains.erase(it);
			cout << "successful delete" << endl;
			cout << buffer  << endl;
			(*it)->setZeros(h, grid_blowup);
		}
		else {(**it).setZeros(h,grid_blowup); it++;}
		
		// test if the grain diappears in the current timestep
		// if false, we must update the neighbors in the next comparison step
		// this needs us to keep the object!!!
    }
	
	
    // try to add boxes from buffer
    if (!buffer.empty()) {
        for (it = buffer.begin(); it != buffer.end(); it++) {
            cout << "trying to add box " << (*it)->getID() << " to Domain " << id;
            bool insert = true;
            for (it2 = grains.begin(); it2 != grains.end(); it2++) {
                if ((*it)->checkIntersect(*it2)) {
                    insert = false;
                    break;
                }
            }
            if (insert) {
                grains.push_back(*it);
				(*it)->setDomain(this);
				(*it)->copy_distances_to_domain();
				//(*it)->free_memory_distance();
                buffer.erase(it); it--;
                cout << ": success" << endl;
            } else {
                cout << ": failed" << endl;
            }
        }
    }

    
    // check for intersects
 
    if (!grains.empty()) 
		for (it = grains.begin(); it != grains.end()-1; ++it) {
			for (it2 = it+1; it2 != grains.end(); ++it2) {
				// on intersect ad box to buffer and erase from grain list
				if ((*it)->checkIntersect(*it2)) {
					cout << "found intersecting box " << (*it)->getID() << " in Domain " << id << endl;
					(*it)->copy_distances();
					buffer.push_back(*it);
					grains.erase(it); it--;
					break;
				}
			}
		}
    else return false; // falls grains.empty() == true
    return true;

}


void matrix::euler(double dt, double h){
	vector<LSbox*>::iterator it;
    for(it = grains.begin(); it != grains.end(); it++)
    {
       (*it)->euler_forward(dt, h);
    }

}



/*********************************************************************************/
// FOURIER TRANSFORMATION + Helperfunctions
//
/*********************************************************************************/

// array must be allocated before!
void matrix::matrix_to_array(double *u){
	int n = get_n();
	int m = get_m();
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			u[i*m+j] = (*x[i])[j];
}

void matrix::array_to_matrix(double *u){
	m = get_m();
	n = get_n();
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			(*x[i])[j] = u[i*m+j];
}

void matrix::makeFFTPlans(double *u, fftw_complex *fftTemp, fftw_plan *fftplan1, fftw_plan *fftplan2)
{ /* creates plans for FFT and IFFT */
	//double *u
	int n = get_n();
	int m = get_m();
	*fftplan1 = fftw_plan_dft_r2c_2d(m,n,u,fftTemp,FFTW_ESTIMATE);
	*fftplan2 = fftw_plan_dft_c2r_2d(m,n,fftTemp,u,FFTW_ESTIMATE);
}

void matrix::conv_generator(double *u, fftw_complex *fftTemp, fftw_plan fftplan1, fftw_plan fftplan2, double dt)
{
	/* Function returns in u the updated value of u as described below..
	u -> (G_{dt})*u
	Assumptions:
	fftplan1 converts u to it's FT (in fftTemp), and
	fftplan2 converts the FT (after pointwise multiplication with G)
	back to a real-valued level set function at u.
	Memory is already allocated in fftTemp
	(necessary to create the plans) */
	
	int n = get_n();
	int m = get_m();
	int i, j;
	int n2 = floor(n/2) + 1;
	double nsq = (double) n * (double) n;
	double k = 2.0 * PI / (double)n;
	double G;
	
	fftw_execute(fftplan1);
	
	for(i=0;i<n2;i++)
		for(j=0;j<n;j++){
			// 	  G= exp((-2.0 * dt) * nsq * (2.0-cos(k*i)-cos(k*j)));
			
			G = 2.0*(2.0 - cos(k*i) - cos(k*j)) * nsq;
			G = 1.0/(1.0+dt*G) / nsq;
			//        USE this line for Richardson-type extrapolation
			//       G = (4.0/pow(1+1.5*(dt)/40*G,40) - 1.0 / pow(1+3.0*(dt)/40*G,40)) / 3.0 / (double)(n*n);
			
			/* normalize G by n*n to pre-normalize convolution results */
			fftTemp[i+n2*j][0] = fftTemp[i+n2*j][0]*G;
			fftTemp[i+n2*j][1] = fftTemp[i+n2*j][1]*G;
		}
	fftw_execute(fftplan2);
}

void matrix::convolution(const double dt, LSbox ***ID){
	int n = get_n();
	int m = get_m();
	double *u, *v;
	
	fftw_complex *fftTemp;
	fftw_plan fwdPlan, bwdPlan;
	matrix diff(m,n);
	
	fftTemp = (fftw_complex*) fftw_malloc(n*(floor(n/2)+1)*sizeof(fftw_complex));
	u = (double*)	fftw_malloc(m*n*sizeof(double));
	v = (double*)	fftw_malloc(m*n*sizeof(double));
	matrix_to_array(u);
	char buffer;
	makeFFTPlans(u, fftTemp,&fwdPlan,&bwdPlan);
	conv_generator(u,fftTemp,fwdPlan,bwdPlan,dt);
	
	/*********************************************************************************/
	// Velocity Corrector Step: 
	/*********************************************************************************/
	// hier soll energycorrection gerechnet werden.
	// in der matrix steht die ursprünglich distanzfunktion, in dem arry die gefaltete
	// energy_correction();	
	// funktion muss umgeschrieben werden
	
	array_to_matrix(u);
	
	fftw_destroy_plan(fwdPlan);
	fftw_destroy_plan(bwdPlan);
	
}
/*********************************************************************************/
/*********************************************************************************/


/*********************************************************************************/
// berechenung der differenz von ausgangslage und bewegter distanzfunktion
// punktweise repräsentiert die distanz, die krümmung, denn kraft = masse * beschleunigung
// die masse ist normiert also 1, die breschleunigung ist kappa. die arbeit ist also (delta d * kappa)
/*********************************************************************************/

matrix matrix::energy_correction(const LSbox ***&ID){
/*	assert(A.n == B.n);
	assert(A.m == B.m);
	
	// boxweise rechnen:
	// boxen sollen dazu neighbor informationen enthalten
	
	temp = A;
	temp =temp-B;
	
	//for (int i = 0; i < m; i++)
	//	for (int j = 0; j < n; j++) {
// 	wie komme ich an die korrekten ID???
// 		temp[i,j] = temp[i,j] * ST[A.id + PARTICLES* B.id];
			
	//	}
	*/
	matrix temp(m,n); 
	return (temp);
}






/*********************************************************************************/
// discrete_convolution is a function to compute a discrete convolution by a grid_blowup times grid_blowup
// kernel in space coordinates. to handle the distancefunction correctly at the boundary we expand
// the domain by grid_blow gridpoints at each boundary.
/*********************************************************************************/

//TO DO: umschreiben auf boxconzept


// Status: Funktion unvollständig!!

bool matrix::discrete_convolution(const double dt, const double h, const int grid_blowup, double (*kernel)(double,int,int,int)){
    int m= get_m();
    int n= get_n();
    matrix erg(m,n,id);
    bool exist = false;
    double conv_rad = grid_blowup/2;
    double tube = double(DELTA)+(h*grid_blowup); // sinnlos oder? vergrößert ja den schlauch?? was war hier meine idee?
	const double outside_domain = -2.0;
    //   double tube = double(DELTA)-conv_rad;
    //   erg = *this;
    for (int i=0; i< m; i++)
        for (int j=0; j< n; j++)
			if((*x[i])[j] > - DELTA  && ((grid_blowup < i) && (i < (m- grid_blowup))) && ((grid_blowup < j) && (j < (n- grid_blowup)))) {		
// 				berechne krümmung im punkt i,j (kappa)
// 				vn= mu + gamma*kappa
// 				x[i][j]= vn* dt + x[i][j]
            }
            
    *this= erg;
	return(exist);
}



/*********************************************************************************/
//maximum of two matrices
/*********************************************************************************/

void matrix::maximum(const matrix &A, const matrix &B){
	assert(A.n == B.n);
	assert(A.m == B.m);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++) {
			if (A[i][j] > B[i][j]) (*x[i])[j] = A[i][j];
			else (*x[i])[j] = B[i][j];
		}
}

int matrix::minimumInPoint(std::list<matrix> distances, int m, int n, int neglect){
	
	std::list<matrix>::iterator it;
	int minID = -1;
	double minVal = 100000; // just a value that can not be reached
	
	for(it = distances.begin(); it != distances.end(); it++){
		if ((*it).id != neglect) {
			if (abs((*it)[m][n]) < minVal) {
				minVal = abs((*it)[m][n]);
				minID = (*it).id;
			}
		}
	}
	
	return minID;
}


/*********************************************************************************/
// Comparison Step
// Comparison step 0.5 *(A_k(x) - max A_i(x) | i!=k)*/
/*********************************************************************************/

void matrix::comparison(std::list<matrix> distances, int grid_blowup){
	std::list<matrix>::iterator it = distances.begin();
	std::list<matrix>::iterator it_current_grain;
	
	int m = get_m();
	int n = get_n();
	matrix Max(m,n);

	if (id == (*it).id) Max = *(++it);
	else Max = *it;
	
	for (it = distances.begin(); it != distances.end(); it++ ){
		if (id != (*it).id) Max.maximum(Max,*it);
		else it_current_grain = it;
	}
	
	Max = ((*it_current_grain)-Max);
	Max.mult_with_scalar(0.5);
	*this = Max;
	
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if ((i <= grid_blowup) || (m-grid_blowup <= i) || (j <= grid_blowup) || (n-grid_blowup <= j)) {
				(*this)[i][j] = INTERIMVAL;
			}
		}
	}
}

void matrix::set_border_to_INTERIMVAL(int grid_blowup)
{
  for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if ((i <= grid_blowup) || (m-grid_blowup <= i) || (j <= grid_blowup) || (n-grid_blowup <= j)) {
				(*this)[i][j] = INTERIMVAL;
			}
		}
	}
}
        
/*********************************************************************************/
// Redistancing für alle Boxen -> greift auf LSbox::redistancing zu:
/*********************************************************************************/

void matrix::clear_domain(double value){
	for (int i = 0; i < m; i++) { // alle EintrÃ¤ge der Ergebnismatrix 0 setzen
        for (int j = 0; j < n; j++) (*this)[i][j] = value;
	}
}


void matrix::redistancing_for_all_boxes(double h, int grid_blowup){
	vector<LSbox*>::iterator it;
	for(it = grains.begin(); it != grains.end(); it++)
	{
		// find zeros and new box size
		(*it)->redistancing(h, grid_blowup);
		cout << "box complete" << endl;
	}
			
}
        
        
/*********************************************************************************/
// Redistancing Standard:

// four runs through a 2-dim grid
// forward x, backward x, forward y, backward y
// remember that iteration i is running in y-direction and vice versa!!!!
/*********************************************************************************/
        
void matrix::redistancing(double h, int grid_blowup){
	int n = get_n();
	int m = get_m();
	matrix *temp = new matrix(m,n,id);
	
	double limiter = -INTERIMVAL;
	double slope = 1;
	// x-direction forward
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n-1; j++) {
			if (j==0) (*temp)[i][j] = -limiter;
			(*temp)[i][j+1] = limiter * utils::sgn((*this)[i][j+1]); // set temp to limiter initially
			
			// check for sign change
			if ((*this)[i][j] * (*this)[i][j+1] < 0.0) {
				// interpolate
				double i_slope  = ((*this)[i][j+1] - (*this)[i][j]) / h;
				double zero = -(*this)[i][j] / i_slope;
				if ( abs((*temp)[i][j]) > abs(-zero)) (*temp)[i][j] = -zero * utils::sgn(i_slope);
			}
			// calculate new distance candidate and assign if appropriate
			double candidate = (*temp)[i][j] + (utils::sgn((*this)[i][j+1]) * h);
			if (abs(candidate) < abs((*temp)[i][j+1])) (*temp)[i][j+1] = candidate;
		}
	}
	
	// x-direction backward
	for (int i = 0; i < m; i++) {
		for (int j = n-1; j > 0; j--) {
			// calculate new distance candidate and assign if appropriate
			double candidate = (*temp)[i][j] + (utils::sgn((*this)[i][j-1]) * h); // replace with the "a"-slope stuff...
			if (abs(candidate) < abs((*temp)[i][j-1])) (*temp)[i][j-1] = candidate;
		}
	}
	
	// y-direction forward
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < m-1; i++) {	
			// check for sign change
			if ((*this)[i][j] * (*this)[i+1][j] < 0.0) {                
				// interpolate
				double i_slope  = ((*this)[i+1][j] - (*this)[i][j]) / h;
				double zero = -(*this)[i][j] / i_slope;
				if ( abs((*temp)[i][j]) > abs(-zero)) (*temp)[i][j] = -zero * utils::sgn(i_slope);
			}
			// calculate new distance candidate and assign if appropriate
			double candidate = (*temp)[i][j] + (utils::sgn((*this)[i+1][j]) * h);
			if (abs(candidate) < abs((*temp)[i+1][j])) (*temp)[i+1][j] = candidate;
		}
	}
	
	// y-direction backward
	for (int j = 0; j < n; j++) {
		for (int i = m-1; i > 0; i--) {	
			// calculate new distance candidate and assign if appropriate
			double candidate = (*temp)[i][j] + (utils::sgn((*this)[i-1][j]) * h); // replace with the "a"-slope stuff...
			if (abs(candidate) < abs((*temp)[i-1][j])) (*temp)[i-1][j] = candidate;
		}
	}
	
}



/*********************************************************************************/
// Redistancing another old Version for the standard case:

// this version needs only to runs over the grid
// therefor we always look in two directions


// BESTE VARIANTE!!!
/*********************************************************************************/
        
void matrix::redistancing_2(double h, int grid_blowup){
	int n = get_n();
	int m = get_m();
	matrix *temp = new matrix(m,n,id);
	
	double limiter = -INTERIMVAL;
	double slope = 1;
	
	// temporary matrix
	(*temp)[0][0]=-limiter;
	for (int i = 0; i < m-1; i++) {
		for (int j = 0; j < n-1; j++){
			// rectangle comparison from upper left corner
			// needs forward declaration in the first row:
			if (i==0 || j==m-2) {
				(*temp)[i][j+1]=-limiter;
			}
			// sign change in x direction
			if ((*this)[i][j] * (*this)[i][j+1] <= 0.0) {
				// interpolate
				double slope  = ((*this)[i][j+1] - (*this)[i][j])  / h;
				double zero_x = -(*this)[i][j] / slope;
				if ( abs((*temp)[i][j]) > abs(-zero_x)) (*temp)[i][j] = -zero_x * utils::sgn(slope);
			}
			
			// sign change in y direction
			if ((*this)[i][j] * (*this)[i+1][j] <= 0.0) {
				// interpolate
				double slope  = ((*this)[i+1][j] - (*this)[i][j])  / h;
				double zero_y = -(*this)[i][j] / slope;
				if (abs((*temp)[i][j]) > abs(-zero_y)) (*temp)[i][j] = -zero_y * utils::sgn(slope);
			}
			
			// calculate new distance candidate and assign if appropriate
			double candidate_x = (*temp)[i][j] + (utils::sgn((*this)[i][j+1]) * h);
			if (abs(candidate_x) < abs((*temp)[i][j+1])) (*temp)[i][j+1] = candidate_x;
			
			// y direction
			// initial "forward"-value in y-direction, depending on sign of respective compared value
			
			
			double candidate_y = (*temp)[i][j] + (utils::sgn((*this)[i+1][j]) * h);
			(*temp)[i+1][j] = limiter * utils::sgn((*this)[i+1][j]);
			if (abs(candidate_y) < abs((*temp)[i+1][j])) (*temp)[i+1][j] = candidate_y;
		}
	}
	
	(*temp)[m-1][n-1]=(*temp)[m-1][n-2]-h;
	// assign temporary matrix to this matrix
	
	for (int i = m-1; i > 0; i--) {
		for (int j = n-1; j > 0; j--){
			double candidate_x = (*temp)[i][j] + (utils::sgn((*this)[i][j-1]) * h);
			if (abs(candidate_x) < abs((*temp)[i][j-1])) (*temp)[i][j-1] = candidate_x;
			
			double candidate_y = (*temp)[i][j] + (utils::sgn((*this)[i-1][j]) * h);
			if (abs(candidate_y) < abs((*temp)[i-1][j])) (*temp)[i-1][j] = candidate_y;
		}
	}
	(*temp)[0][0]=(*temp)[0][1]-h;
	*this=*temp;
	delete temp;
	
	
}     
        
        
        
        
        
        
/*********************************************************************************/
// Redistancing_Advanced -> This method uses a Slopefield to 
// 		computed the distances with a slopefactor
/*********************************************************************************/

void matrix::redistancing_advanced(double h, int grid_blowup, std::list<matrix> distances, double** borderSlopes, double** slopeField) {
	int n = get_n();
	int m = get_m();
	matrix *temp = new matrix(m,n,id);
	
	double limiter = DELTA;
	double slope = 1;
				
	// x-direction forward
	for (int i = 0; i < m; i++) {
		slope = 1;
		for (int j = 0; j < n-1; j++) {
			if (j==0) (*temp)[i][j] = -limiter;
			(*temp)[i][j+1] = limiter * utils::sgn((*this)[i][j+1]); // set temp to limiter initially
			
			// check for sign change
			if ((*this)[i][j] * (*this)[i][j+1] < 0.0) {
				// find grain with minimal distance to [i][j]
				int rightID = (*this).id;
				int leftID = minimumInPoint(distances, i, j, rightID);
				slope = borderSlopes[leftID][rightID];				
				if (slope == 0) slope = 1;				
				// interpolate
				double i_slope  = ((*this)[i][j+1] - (*this)[i][j]) / h;
				double zero = -(*this)[i][j] / i_slope;
				if ( abs((*temp)[i][j]) > abs(-zero)) (*temp)[i][j] = -zero * utils::sgn(i_slope);
			}
			// calculate new distance candidate and assign if appropriate
			double candidate = (*temp)[i][j] + (utils::sgn((*this)[i][j+1]) * h * slope); 
			if (abs(candidate) < abs((*temp)[i][j+1])) (*temp)[i][j+1] = candidate;
		}
	}
	
	// y-direction forward
	for (int j = 0; j < n; j++) {
		slope = 1;
		for (int i = 0; i < m-1; i++) {			
			// check for sign change
			if ((*this)[i][j] * (*this)[i+1][j] < 0.0) {
				// find grain with minimal distance to [i][j]
				int bottomID = (*this).id;
				int topID = minimumInPoint(distances, i, j, bottomID);
				slope = borderSlopes[topID][bottomID];				
				if (slope == 0) slope = 1;				
				// interpolate
				double i_slope  = ((*this)[i+1][j] - (*this)[i][j]) / h;
				double zero = -(*this)[i][j] / i_slope;
				if ( abs((*temp)[i][j]) > abs(-zero)) (*temp)[i][j] = -zero * utils::sgn(i_slope);
			}
			// calculate new distance candidate and assign if appropriate
			double candidate = (*temp)[i][j] + (utils::sgn((*this)[i+1][j]) * h * slope);
			if (abs(candidate) < abs((*temp)[i+1][j])) (*temp)[i+1][j] = candidate;
		}
	}
	
	// x-direction backward
	for (int i = 0; i < m; i++) {
		slope = 1;
		for (int j = n-1; j > 0; j--) {			
			// check for sign change
			if ((*this)[i][j] * (*this)[i][j-1] < 0.0) {
				// find grain with minimal distance to [i][j]
				int leftID = (*this).id;
				int rightID = minimumInPoint(distances, i, j, leftID);
				slope = borderSlopes[leftID][rightID];
				
				if (slope == 0) slope = 1;
			}			
			// calculate new distance candidate and assign if appropriate
			double candidate = (*temp)[i][j] + (utils::sgn((*this)[i][j-1]) * h * slope); // replace with the "a"-slope stuff...
			if (abs(candidate) < abs((*temp)[i][j-1])) (*temp)[i][j-1] = candidate;
		}
	}
	
	
	// y-direction backward
	for (int j = 0; j < n; j++) {
		slope = 1;
		for (int i = m-1; i > 0; i--) {			
			// check for sign change
			if ((*this)[i][j] * (*this)[i-1][j] < 0.0) {
				// find grain with minimal distance to [i][j]
				int topID = (*this).id;
				int bottomID = minimumInPoint(distances, i, j, topID);
				slope = borderSlopes[topID][bottomID];
				
				if (slope == 0) slope = 1;
			}			
			// calculate new distance candidate and assign if appropriate
			double candidate = (*temp)[i][j] + (utils::sgn((*this)[i-1][j]) * h * slope); // replace with the "a"-slope stuff...
			if (abs(candidate) < abs((*temp)[i-1][j])) (*temp)[i-1][j] = candidate;
		}
	}
	
	*this = *temp;
	delete temp;
}
	
	
	
int matrix::get_m() const { return m; };
int matrix::get_n() const { return n; };
int matrix::get_id() const { return id; };
int matrix::get_nr_of_grains() { return grains.size(); };
	
