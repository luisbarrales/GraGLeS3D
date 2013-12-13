#ifndef utilities_h
#define utilities_h

#include "ggLS.h"


namespace utils {
  
template <class T >
void print_vector(const vector<T> &w) ;

template <class T >
void print_vector(const vector<T> &w) {
	stringstream os;
		for(int i=0; i< (w.size()-1); i++){
		os << ::std::fixed << w[i] << "\t";
	}
	os << ::std::fixed << w[w.size()-1];
	std::cout << os.str() << endl;
}

template <class T >
void print_2dim_array( T* arr, int m , int n);

template <class T >
void print_2dim_array( T* arr, int m , int n){
	stringstream os;
	for (int i = 0; i < m; i++){
	  for (int j = 0; j < n; j++) {
		os << ::std::fixed << arr[i*n+j] << "\t";
		}
		os << endl;
	}
	std::cout << os.str() << endl;
}


template <class T>
void save_2dim_array( T *arr, int m , int n, const char* filename);


template <class T>
void save_2dim_array( T *arr, int m , int n, const char* filename){
    ofstream datei;
    datei.open(filename);
// 	stringstream os;
	  for (int i = 0; i < m; i++){
		for (int j = 0; j < n; j++){
		  datei << ::std::fixed << arr[i*n+j] << "\t";
		}
		datei << endl;
	}
// 	datei << os << endl;
    datei.close();
}

void plotContour(const char *fileName, const char *plotfiles, int len);


inline double rnd() {return double(rand())/RAND_MAX;}

inline int sgn(double x)
{
    if (x==0.0)
        return 0.0;
    else
        return (x>.0) ? 1. : -1.;
}

void plotGnu(const char *fileName, const char *plotfiles, int len);
    void plotGnuPNG(const char *fileName, const char *plotfiles, int len);
    void PNGtoGIF(const char *fileName);

};
#endif
