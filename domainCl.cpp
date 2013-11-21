#include "domainCl.h"

using namespace std;
using namespace voro;


domainCl::domainCl()
{
  val	=	(double*) fftw_malloc ( m*n*sizeof (double)); 
    x	= new double*[m];
}

domainCl::domainCl(int m, int n) : m(m), n(n) {
      
    id	=	0;
    val	=	(double*) fftw_malloc ( m*n*sizeof (double)); 
    x	= new double*[m];
    for (int i=0;i<m;i++) 
      x[i]=&val[i*n];
  
}

//keine ID übergabe erforderlich?

domainCl::domainCl(int m, int n, int id) : m(m), n(n), id(id) {
    val = (double*) fftw_malloc (m*n*sizeof(double)); 
    x	=  new double*[m];
    for (int i=0;i<m;i++) 
      x[i]	=	&val[i*n];
}

domainCl::domainCl(int m, int n, int id, double startval) : m(m), n(n), id(id) {
    val = (double*) fftw_malloc (m*n*sizeof(double)); 
    x	= new double*[m];
    std::fill_n(val, m*n, startval);
    for (int i=0;i<m;i++) 
      x[i]=&val[i*n];
}

domainCl::domainCl(int m, int n, int id, double startval, grainhdl* owner) : m(m), n(n), id(id), owner(owner) {
    val = (double*) fftw_malloc ( m*n*sizeof(double)); 
    x	= new double*[m];
    std::fill_n(val, m*n, startval);
    for (int i=0;i<m;i++) 
      x[i]	=	&val[i*n];
}

domainCl::~domainCl() {
     delete [] x;
     fftw_free(val);
}

domainCl::domainCl(const domainCl& v) : m(v.m), n(v.n), id(v.id) {
//  cerr << "m " ;
    val = (double*) fftw_malloc ( m*n*sizeof(double)); 
    x	=  new double*[m];
    for (int i=0;i<m;i++) 
      x[i]=&val[i*n];

    for (int i=0;i<m;i++) {
      for (int j = 0; j < n; j++){
		val[i*n+j] = v[i][j];
      }
	
    }
}

double* domainCl::operator[](int i) {
    if (0	<= i &&i< m) {
      return x[i];
    }
    

  
}

const double* domainCl::operator[](int i) const {
    if (0	<= i &&i< m) {
      return x[i];
    }    
     
}
// double& domainCl::operator()(int x,int y)  {
//     return val[x*m+n]; 
// }
// const double& domainCl::operator()(int x,int y) const {
//     return val[x*m+n];
//         
// }

// "=" Operator
double& domainCl::operator=(const domainCl& A){
//   cerr << "here" << endl;
    if (this != &A) {
    assert(m == A.m && n == A.n);
	for (int i=0;i<m;i++) {
	  for (int j=0;j<n;j++)
	    val[i*n+j]=A[i][j];
	}
    }
}



ostream& operator<<(ostream &os, const domainCl& A) {
    int i;
	for(i=0; i< (A.m); i++){
	  for (int j=0; j< A.n; j++){
		os << A.val[i*A.n+j]<< "\t";
	  }
	  os << endl;
	}
    return os;
}

void domainCl::save_domainCl(const char* filename){
    ofstream datei;
    datei.open(filename);
    datei << *this << endl;
    datei.close();
}

double domainCl::entry(const int i, const int j){
    return(x[i][j]);
}

void domainCl::mult_with_scalar(const double d){
  for (int i = 0; i < m*n; i++) 
	  val[i] *= d;
}
 
 domainCl domainCl::operator-(domainCl& A){
    domainCl C(m,n,id);
    C.grains = this->grains;
    assert(m == A.m && n == A.n);
    for (int i=0; i<m; i++) 
      for (int j=0; j<n; j++)
		C[i][j] = (*this)[i][j] - A[i][j];
    return C;
}


domainCl domainCl::distancefunction(voronoicell_neighbor& c, int *gridIDs, double *part_pos, int grid_blowup, double h){
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
            p[0]=(i-grid_blowup)*h; 
			p[1]=(j-grid_blowup)*h;
            for(int ii=0;ii<c.p;ii++) {
				for(int jj=0;jj<c.nu[ii];jj++) {
                 
                    k=c.ed[ii][jj];
                    
                    x1[0]=vv[3*ii];x1[1]=vv[3*ii+1];
                    x2[0]=vv[3*k]; x2[1]=vv[3*k+1];
                    
                    if (x1!=x2){
                        a=x1;
                        u = x2-x1;
                        lambda=((p-a)*u)/(u*u);
                        
                        if(lambda <= 0)
			  d= (p-x1).laenge();
                        if((0. < lambda) && (lambda < 1.)) 	
			  d= (p-(a+(u*lambda))).laenge();
                        if(lambda >= 1.) 					    
			  d= (p-x2).laenge();
//                     (ID_mat[0][i*m +j])->getID()    
                        if(id==gridIDs[i*m +j] && ((grid_blowup < i) && (i < (m- grid_blowup))) && ((grid_blowup < j) && (j < (m- grid_blowup)))) 
			  d=abs(d);
                        else d= -abs(d);
                        
                        if(abs(d)< abs(dmin)) dmin=d;
                    }
                }
            }
            //  if (abs(dmin) < DELTA || dmin > 0.0) (*x[i])[j]= dmin;
            //else (*x[i])[j] = -DELTA;
            x[i][j]= dmin;
        }
	}
	return(*this);
}



/*************************************************************************/
// BOX Functions
/*************************************************************************/

bool domainCl::addBox(LSbox* aBox){
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


vector<LSbox*> domainCl::getBoxList() {
    return grains;
}



bool domainCl::grainCheck(double h, int grid_blowup, vector<LSbox*>& buffer, int loop){
	char buffer1;
	bool exist;
	vector<LSbox*>::iterator it,it2;
    for(it = grains.begin(); it != grains.end();){
		if((**it).get_status()==false) {
			cout << "try to delete grain " << (**it).get_id() << " in domain "<< id << endl;
			// the grain has disappeared the timestep before
			//now we can clean up the memory
			delete (*it); grains.erase(it);
			cout << "successful delete" << endl;
		}
		else {(**it).setZeros(h,grid_blowup, loop); it++;}
		
		// test if the grain diappears in the current timestep
		// if false, we must update the neighbors in the next comparison step
		// this needs us to keep the object!!!
    }
	
	
    // try to add boxes from buffer
    if (!buffer.empty()) {
        for (it = buffer.begin(); it != buffer.end(); ) {
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

                buffer.erase(it);
                cout << ": success" << endl;
            } else {
				it++;
                cout << ": failed" << endl;
            }
        }
    }

    bool del = false;
    // check for intersects
    if (!grains.empty()) 		
		for (it = grains.begin(); it != grains.end()-1;) {
			for (it2 = it+1; it2 != grains.end(); ++it2) {
				// on intersect add box to buffer and erase from grain list				
				if ((*it)->checkIntersect(*it2)) {
					cout << "found intersecting box " << (*it)->getID() << " in Domain " << id << endl;
					(*it)->copy_distances();
					buffer.push_back(*it);
					grains.erase(it); 
					del = true;
					break;
				}
// 				if(it2 != grains.end()) it++;
			}
			if (!del) it++;
			del=false;
		}
    else return false; // falls grains.empty() == true
    return true;

}





/*********************************************************************************/
// FOURIER TRANSFORMATION + Helperfunctions
//
/*********************************************************************************/

void domainCl::makeFFTPlans(double *u, fftw_complex *fftTemp, fftw_plan *fftplan1, fftw_plan *fftplan2)
{ /* creates plans for FFT and IFFT */
	int nn= owner->get_ngridpoints();
	*fftplan1 = fftw_plan_dft_r2c_2d(nn,nn,u,fftTemp,FFTW_ESTIMATE);
	*fftplan2 = fftw_plan_dft_c2r_2d(nn,nn,fftTemp,u,FFTW_ESTIMATE);
}

void domainCl::conv_generator(double *u, fftw_complex *fftTemp, fftw_plan fftplan1, fftw_plan fftplan2, double dt)
{
	/* Function returns in u the updated value of u as described below..
	u -> (G_{dt})*u
	Assumptions:
	fftplan1 converts u to it's FT (in fftTemp), and
	fftplan2 converts the FT (after pointwise multiplication with G)
	back to a real-valued level set function at u.
	Memory is already allocated in fftTemp
	(necessary to create the plans) */
	
	int n= owner->get_ngridpoints();
	int i, j;
	int n2 = floor(n/2) + 1;
	double nsq = n *  n;
	double k = 2.0 * PI / n;
	double G;
	double coski;
	fftw_execute(fftplan1);
	
	for(i=0;i<n2;i++) {
	  coski=cos(k*i);
		for(j=0;j<n;j++){
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





void domainCl::convolution(const double dt, double *ST, LSbox ***ID, domainCl &ref, LSbox* zeroBox, int grid_blowup, weightmap* my_weights){
	int n= owner->get_ngridpoints();

	double *u, *v;
	double vn, vnn;
	
	fftw_complex *fftTemp;
	fftw_plan fwdPlan, bwdPlan;
	
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
	// in der domainCl steht die ursprünglich distanzfunktion, in dem arry die gefaltete
	if(!ISOTROPIC){
		double rad =  DELTA;
		double weight;
// 		int* rep = new int[3];
		vector<LSbox*>::iterator it;
		for(it = grains.begin(); it != grains.end(); it++){
// 			bool out=false;
			for (int i = (**it).ymin; i < (**it).ymax; i++){
				for (int j = (**it).xmin; j < (**it).xmax; j++) {
					if( ID[0][i*m +j] != zeroBox ){
						if (ID[0][i*m +j]!=ID[1][i*m +j] && ID[1][i*m +j]!=ID[2][i*m +j]){
							if(rad < abs(ref[i][j])) weight = 1.0;
							else {
								weight = (*my_weights).load_weights(m, ST, ID,i,j,(**it).get_id());
								weight = ( 1-abs(rad - abs(ref[i][j])) )* weight;
							}
// 							weight = (*my_weights).load_weights(m, ST, ID,i,j,(**it).get_id());
// 							cout << weight << endl;
							(*this)[i][j] = ref[i][j] + (((*this)[i][j] -ref[i][j]) * weight);
						}
					}
				}
			}
		}
	}
}
/*********************************************************************************/
/*********************************************************************************/



/*********************************************************************************/
//maximum of two matrices
/*********************************************************************************/

void domainCl::maximum(const domainCl &A, const domainCl &B){
	assert(A.n == B.n);
	assert(A.m == B.m);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++) {
			if (	A[i][j] > B[i][j]) x[i][j] = A[i][j];
			else 	x[i][j] = B[i][j];
		}
}

int domainCl::minimumInPoint(std::list<domainCl> distances, int m, int n, int neglect){
	
	std::list<domainCl>::iterator it;
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

void domainCl::comparison(std::list<domainCl> distances, int grid_blowup){
	
	std::list<domainCl>::iterator it = distances.begin();
	std::list<domainCl>::iterator it_current_distance;
	int n= owner->get_ngridpoints();
	
	domainCl Max(m,n), cur_Max(m,n,id);

	if (id == (*it).id) Max = *(++it);
	else Max = *it;
	for (it = distances.begin(); it != distances.end(); it++ ){
		  if (id != (*it).id) Max.maximum(Max,*it);
		 else it_current_distance = it;
	}
	
	cur_Max = ((*it_current_distance)-Max);
	cur_Max.mult_with_scalar(0.5);
	*this = cur_Max;
	
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if ((i <= grid_blowup) || (m-grid_blowup <= i) || (j <= grid_blowup) || (n-grid_blowup <= j)) {
				(*this)[i][j] = INTERIMVAL;
			}
		}
	}
}

void domainCl::set_border_to_INTERIMVAL(int grid_blowup)
{
// 	double h = owner->get_h();
//   for (int i = 0; i < m; i++) {
// 		for (int j = 0; j < n; j++) {
// 			if ((i <= grid_blowup) || (m-grid_blowup <= i) || (j <= grid_blowup) || (n-grid_blowup <= j)) {
// 				(*this)[i][j] = -DELTA;
// 			}
// 			if ((i == grid_blowup+1) || (m-grid_blowup-1== i) || (j == grid_blowup+1) || (n-grid_blowup-1 == j))
// 				(*this)[i][j] = h;
// 		}
// 	}
}
        
/*********************************************************************************/
// Redistancing für alle Boxen -> greift auf LSbox::redistancing zu:
/*********************************************************************************/

void domainCl::clear_domain(double value){
	std::fill_n(val, m*n, value);
}


void domainCl::redistancing_for_all_boxes(double h, int grid_blowup){
	vector<LSbox*>::iterator it;
	for(it = grains.begin(); it != grains.end(); it++)
	{
		// find zeros and new box size
		(*it)->redist_box(h, grid_blowup);
	}
			
}
        
        
/*********************************************************************************/
// Redistancing Standard:

// four runs through a 2-dim grid
// forward x, backward x, forward y, backward y
// remember that iteration i is running in y-direction and vice versa!!!!
/*********************************************************************************/
        
void domainCl::redistancing(double h, int grid_blowup){
	int n = get_n();
	int m = get_m();
	
	domainCl *temp = new domainCl(m,n,id,-1.0);
	double slope = 1;
	double candidate, i_slope,zero;
	// x-direction forward
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n-1; j++) {
			//check for sign change
			if ((*this)[i][j] * (*this)[i][j+1] <= 0.0) {
				// interpolate
				i_slope  = ((*this)[i][j+1] - (*this)[i][j]) / h;
				zero = -(*this)[i][j] / i_slope;
				if ( abs((*temp)[i][j]) > abs(zero)) (*temp)[i][j] = -zero * utils::sgn(i_slope);

			}
			// calculate new distance candidate and assign if appropriate
			candidate = (*temp)[i][j] + (utils::sgn((*this)[i][j+1]) * h);
			if (abs(candidate) < abs((*temp)[i][j+1])) (*temp)[i][j+1] = candidate;
		}
	}
	
	// x-direction backward
	for (int i = 0; i < m; i++) {
		for (int j = n-1; j > 0; j--) {
			if ((*this)[i][j] * (*this)[i][j-1] <= 0.0) continue;
			// calculate new distance candidate and assign if appropriate
			candidate = (*temp)[i][j] + (utils::sgn((*this)[i][j-1]) * h); // replace with the "a"-slope stuff...
			if (abs(candidate) < abs((*temp)[i][j-1])) (*temp)[i][j-1] = candidate;
		}
	}
	
	// y-direction forward
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < m-1; i++) {	
			// check for sign change
			if ((*this)[i][j] * (*this)[i+1][j] <= 0.0) {   
					// interpolate
					i_slope  = ((*this)[i+1][j] - (*this)[i][j]) / h;
					zero = -(*this)[i][j] / i_slope;
					if ( abs((*temp)[i][j]) > abs(zero)) (*temp)[i][j] = -zero * utils::sgn(i_slope);
			}
			// calculate new distance candidate and assign if appropriate
			candidate = (*temp)[i][j] + (utils::sgn((*this)[i+1][j]) * h);
			if (abs(candidate) < abs((*temp)[i+1][j])) (*temp)[i+1][j] = candidate;
		}
	}
	
	// y-direction backward
	for (int j = 0; j < n; j++) {
		for (int i = m-1; i > 0; i--) {	
			// calculate new distance candidate and assign if appropriate
			candidate = (*temp)[i][j] + (utils::sgn((*this)[i-1][j]) * h); // replace with the "a"-slope stuff...
			if (abs(candidate) < abs((*temp)[i-1][j])) (*temp)[i-1][j] = candidate;
		}
	}
	
	*this=*temp;
	delete temp;
	
}



/*********************************************************************************/
// Redistancing another old Version for the standard case:

// this version needs only to runs over the grid
// therefor we always look in two directions


// BESTE VARIANTE!!!
/*********************************************************************************/
        
void domainCl::redistancing_2(double h, int grid_blowup){
	int n = get_n();
	int m = get_m();
	domainCl *temp = new domainCl(m,n,id,-1.0);
	
	double limiter = -INTERIMVAL;
	double slope = 1;
	
	// temporary domainCl
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
	// assign temporary domainCl to this domainCl
	
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
        
	
int domainCl::get_m() const { return m; };
int domainCl::get_n() const { return n; };
int domainCl::get_id() const { return id; };
int domainCl::get_nr_of_grains() { return grains.size(); };
	
