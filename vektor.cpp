#include "vektor.h"
 
using namespace std;


vektor::vektor(int n) : n(n) {
        x=new double[n];
        for (int i=0;i<n;i++) x[i]=0;
    }

vektor::vektor(int n, double startval) : n(n) {
        x=new double[n];
        for (int i=0;i<n;i++) x[i]=startval;
    }


vektor::vektor(const vektor& v) : n(v.n) {
        x=new double[n];
        for (int i=0;i<n;i++) x[i]=v.x[i];
}


vektor::~vektor() { 
	delete [] x; 
}
/*****************************************************************************/
/*****************************************************************************/
// overloaded operators
/*****************************************************************************/


double& vektor::operator[](int i) {
    if (0<= i < n) return x[i];
    else {
        outOfBoundsException e(i, i);
        cout << "invalid index " << e.what() << endl;
    }
}

const double& vektor::operator[](int i) const {
    if (0<= i < n) return x[i];
    else {
        outOfBoundsException e(i, i);
        cerr << "invalid index " << e.what() << endl;
    }
}

    // "=" Operator Zuweisung
double& vektor::operator=(const vektor& w){
        if (this!= &w) {
            assert(n == w.n);
            for (int i = 0; i < n; i++) x[i] = w.x[i];
        }
}

    // "+" Operator
vektor vektor::operator+(const vektor& w){
    vektor erg(n);
    assert(n == w.n);
    for (int i = 0; i < n; i++) erg[i] = x[i] + w.x[i];
    return erg;
}

    // "-" Operator
vektor vektor::operator-(const vektor& w){
    vektor erg(n);
    assert(n == w.n);
    for (int i = 0; i < n; i++) erg[i] = x[i] - w.x[i];
    return erg;
}





bool vektor::operator==(const vektor& w){
	 assert(n == w.n);
       for (int i = 0; i < n; i++) 	if( abs(x[i]-w.x[i]) > EPS) return(false);
	 return(true);
}

bool vektor::operator!=(const vektor& w){
	 assert(n == w.n);
       for (int i = 0; i < n; i++) 	{
		 if( abs(x[i]-w.x[i]) > EPS ) return(true);
	   }
	 return(false);
}


/*****************************************************************************/
/*****************************************************************************/
// friend functions
/*****************************************************************************/


ostream& operator<<(ostream &os, const vektor& w) {
    int i;
	for(i=0; i< (w.n-1); i++){
		os << ::std::fixed << w.x[i] << "\t";
	}
	os << ::std::fixed << w.x[w.n-1];
      return os;
}

//berechnet Abstand einer Geraden, gegeben durch den Noralvektor n und einen bel. Punkt a, zum Punkt p
double hessedistanz(vektor& n, const vektor& a, const vektor& p){
	double d;
	assert(a.n==p.n);
	assert(n.n==p.n);
	n.Normiere();
	d=(n*p)-(n*a);
	return(d);
}

/*****************************************************************************/
/*****************************************************************************/
// class functions
/*****************************************************************************/


    // Skalarprodukt
double vektor::operator*(const vektor& w){ // kein "&" hinter double
    double erg = 0.0;
    for (int i = 0; i < n; i++) erg += x[i] * w.x[i];
    return erg;
}


    // Kreuzprodukt für dim=2,3
vektor vektor::XProdukt(const vektor& v, const vektor& w){ 
	assert(v.n==w.n && ( n == 2 || n == 3 ));
	if(n==2){
		x[1]= v.x[0]-w.x[0];
		x[0]= -(v.x[1]-w.x[1]);
		Normiere();
		}
	if(n==3){
		x[0] = v.x[1]*w.x[2] - v.x[2]*w.x[1];
		x[1] = v.x[2]*w.x[0] - v.x[0]*w.x[2];
		x[2] = v.x[0]*w.x[1] - v.x[1]*w.x[0];
		Normiere();
		}
}


double vektor::distance_to_line(vektor& x1, vektor& x2){
  vektor lfp(2);
  vektor n(2), a(2);
  n= x2-x1;
  a=x1;
  double lambda;
  lambda= ((*this*n)- (n*a))/ (n*n);
  lfp = x1 + x2*lambda;
  return((*this-lfp).laenge());  
}



int vektor::get_n() const { return n; };
