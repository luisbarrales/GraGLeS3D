// Single Voronoi cell example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#include "voro++.hh"
#include <iostream>
#include <vector>

using namespace voro;

// This function returns a random floating point number between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {
	double x,y,z,rsq,r;
	int k;
	voronoicell v;

	// Initialize the Voronoi cell to be a cube of side length 2, centered
	// on the origin
	v.init(-1,1,-1,1,0,0);

	// Cut the cell by 250 random planes which are all a distance 1 away
	// from the origin, to make an approximation to a sphere
	for(int i=0;i<10;i++) {
		x=2*rnd()-1;
		y=2*rnd()-1;
		z=0;//2*rnd()-1;
		rsq=x*x+y*y+z*z;
		if(rsq>0.01&&rsq<1) {
			r=1/sqrt(rsq);x*=r;y*=r;z*=r;
// 			cout << x << " ||  " << y << " ||  " << z << endl;
			v.plane(x,y,z,1);
		}
	}
	
	vector<int> fv;
	v.face_vertices(fv);
// 	for(int i=0; i<int( 3*v.p);i++) cout << fv[i] << " || ";
// 	cout << endl << endl;
	vector<double> vv;
	v.vertices(vv);
// 	for(int i=0; i<int( 3*v.p);i++) cout << vv[i] << " || ";
// 	cout << endl << endl;
// 	
	v.centroid(x,y,z);
	for(int i=0;i<v.p;i++) {
	  for(int j=0;j<v.nu[i];j++) {
		k=v.ed[i][j]; 
		cout << k << " : "<< vv[3*k] << " || " << vv[3*k +1] << " || "  << vv[3*k +2] << endl;
		cout << k << " : "<< fv[3*k] << " || " << fv[3*k +1] << " || "  << fv[3*k +2] << endl << endl;
	  }
	  cout << endl << endl;
	}


	for (int i=0; i < 3*v.p; i++) cout << 0.5*v.pts[i] << " || ";
	cout << endl << "centroid:" << endl;
	cout << x << " ||  " << y << " ||  " << z << endl;
	// Output the Voronoi cell to a file, in the gnuplot format
	v.draw_gnuplot(0,0,0,"single_cell.gnu");
	
}
