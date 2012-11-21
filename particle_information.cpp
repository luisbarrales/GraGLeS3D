#include "particle_information.h"

using namespace std;

particle_information::particle_information(){
 particles=0;
 center=NULL;
}

particle_information::particle_information(int particles) : particles(particles) {
  center = new double[3*particles];
}


particle_information::~particle_information() { 
  delete [] center;
}

vektor particle_information::get_coordinates(int id){
 vektor coord(3);
 coord[0]=center[3*id];
 coord[1]=center[3*id+1];
 coord[2]=center[3*id+2];
 cout << id << " | "<< coord << endl;
 return coord;
}

void particle_information::save_info(int id,double rx,double ry,double rz){
  center[3*id]=rx; center[3*id+1]=ry; center[3*id+2]=rz;
}


