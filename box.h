#ifndef BOX_H
#define BOX_H

#include "ggLS.h"

// #include "contour.h"
using namespace std;

class LSbox;
class grainhdl;
class Weightmap;

template<class T>
class DimensionalBuffer;

class MarchingSquaresAlgorithm;

struct SPoint;

struct characteristics{
    LSbox* directNeighbour;
    double length;
    double energyDensity;
    double mis_ori;
	double mobility;
    characteristics(LSbox* directNeighbour, double length, double energyDensity, double mis_ori) : directNeighbour(directNeighbour), length(length), energyDensity(energyDensity), mis_ori(mis_ori), mobility(1)
	{}
	characteristics(LSbox* directNeighbour, double length, double energyDensity, double mis_ori, double mu) : directNeighbour(directNeighbour), length(length), energyDensity(energyDensity), mis_ori(mis_ori), mobility(mu)
	{}
	characteristics( const characteristics& other ) :
		directNeighbour(other.directNeighbour), length(other.length), energyDensity(other.energyDensity), mis_ori(other.mis_ori), mobility(other.mobility)
	{}
};

/*!
 * \class LSBox
 * \brief Class encapsulating a Level Set Box.
 *
 * LSbox class contains the coordinates of the box in the actual grid. <br>
 * For each point that the LSbox covers, it stores: <br>
 * - Distances to the actual grain boundry. <br>
 * - List of pointers to other LSBoxes that influence the point.
 *
 * The class also stores the coordinates of the LSBox, a weightmap,Euler angles that represent
 * the orientation, the volume of the grain, <br>
 * the energy of the grain and a pointer to the \b grainhdl object.
 *
 */
class LSbox {
	unsigned int id;
	int xmaxId ,xminId , ymaxId, yminId; 
	
	bool exist;
	vector<vector<LSbox*>> IDLocal;
	Weightmap* local_weights;

	int nvertices;
	double* quaternion;
	double energy, volume, perimeter;
	grainhdl* handler;
	vector<SPoint> contourGrain;
	vector<characteristics> grainCharacteristics;
	
	DimensionalBuffer<double>* inputDistance;
	DimensionalBuffer<double>* outputDistance;

public:
	friend class grainhdl;
    LSbox();
    ~LSbox();
	vector<LSbox*> neighbourCandidates;
	vector<LSbox*> neighbors_2order;
	LSbox(int id, double phi1, double PHI, double phi2, grainhdl* owner);
    LSbox(int aID, voro::voronoicell_neighbor& c, double *part_pos, grainhdl* owner);
	LSbox(int id, int nvertex, double* vertices, double phi1, double PHI, double phi2, grainhdl* owner);
	void distancefunction(int nvertex, double* vertices);
    void distancefunction(voro::voronoicell_neighbor& c, double *part_pos);
    void redist_box();
	void find_contour();
	int  getID();
	void resizeToSquareIn();
    void resizeToSquareOut();
    void determineIDs();
    void comparison();
	double getDistance(int i, int j);
    void set_comparison();
    void add_n2o();
	void add_n2o_2();
    void computeVolumeAndEnergy();
	double getGBEnergyTimesGBMobility(int i,int j);
	double getGBEnergyTimesGBMobility(LSbox* neighbour);
	double getGBEnergy(LSbox* neighbour);
	
    bool checkIntersect(LSbox* box2);   	
	void free_memory_distance();
      
	void convolution();
	void get_new_IDLocalSize();

	void plot_box_contour(int timestep = -1, bool plot_energy=false, ofstream* dest_file = NULL);
	void saveGrain(ofstream* dest_file);

	void plot_box(bool distanceplot, int select, string simstep, bool local);
	double mis_ori(LSbox* grain_2);
	void checkIntersect_zero_grain();
	void inversDistance();
	
	void makeFFTPlans(double *in, double* out,fftw_complex *fftTemp, fftw_plan *fftplan1, fftw_plan *fftplan2);
	void conv_generator(fftw_complex *fftTemp, fftw_plan fftplan1, fftw_plan fftplan2);
	void switchInNOut();
	void boundaryCondition();
	void updateFirstOrderNeigbors();
	double GBmobilityModel(double thetaMis);
	bool isNeighbour(LSbox* candidate);
	bool BoundaryIntersection();

		
	inline bool get_status(){ return exist;}
	inline int get_id() 	{ return id; }
	inline double get_vol() {return volume ;}
	inline double getEnergy(){ return energy;}
};
#endif
