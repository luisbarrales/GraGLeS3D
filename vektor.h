#ifndef VEKTOR_h
#define VEKTOR_h

#include "ggLS.h"

using namespace std;

 class domainCl;
// class particle_information;

class vektor{
	int	n;
	double	*x;

  public:
	/*****************************************************************************/
	//friends
	/*****************************************************************************/
    friend class domainCl;
	
	/*****************************************************************************/
	// friend functions
	/*****************************************************************************/
	friend ostream& operator<<(ostream& os, const vektor& ws);
	friend double hessedistanz(vektor& n, const vektor& a, const vektor& p);
	
	
	/*****************************************************************************/
	// constructors etc.
	/*****************************************************************************/
	vektor(int n);
	vektor(int n, double startval);
	vektor(const vektor& v);
	~vektor();
	 
	/*****************************************************************************/
	// overloaded operators
	/*****************************************************************************/
    double& 		operator[]	(int i);
    const double& 	operator[]	(int i) const;
    double& 		operator=	(const vektor& w);
    vektor  		operator+	(const vektor& w);
    vektor  		operator-	(const vektor& w);
    double  		operator*	(const vektor& w); // skalarprodukt
	
	template <class T>
	vektor operator*(const T d){
	  vektor erg(n);
	  for(int i=0; i<n; i++) 
	    erg[i]=d*x[i];
	  return(erg);
	} // skalar multiplikation vektor* skalar = vektor
	
	vektor  operator*(const domainCl& A);
	vektor  operator/(domainCl& A);
	bool 	operator==(const vektor& w);
	bool 	operator!=(const vektor& w);
	
	/*****************************************************************************/
	// class functions
	/*****************************************************************************/
	inline void skaliere(const double laenge){
	  for (int i=0; i < n; i++)
		x[i]/=laenge;
	}
	inline double laenge(){
	  double sum = 0; 
	  for(int i=0; i<n; i++) sum+=(x[i]*x[i]);
	  return(sqrt(sum));
	}
	
	
	inline void Normiere() {
	  skaliere(laenge());
	}

	vektor XProdukt(const vektor& v, const vektor& w);
	int get_n() const;
	double distance_to_line(vektor& x1, vektor& x2);
	
};

#endif
