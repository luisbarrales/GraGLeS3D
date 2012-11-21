#ifndef particle_information_h
#define particle_information_h

#include "levelsetproject.h"


class particle_information {
  
  int particles;
  double *center;

public:
  friend class vektor;
  
  particle_information();
  particle_information(int particles);
  ~particle_information();
  
  vektor get_coordinates(int id);
  void save_info(int id,double rx,double ry,double rz);
};


#endif
