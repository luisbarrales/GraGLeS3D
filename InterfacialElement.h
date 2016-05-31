/*
 * InterfacialElement.h
 *
 *  Created on: Nov 24, 2015
 *      Author: cm654063
 */

#ifndef INTERFACIALELEMENTS_H
#define INTERFACIALELEMENTS_H

#include <vector>
#include "triangle.h"
struct NeighborList;

class LSbox;
class grainhdl;
class Settings;
class MarchingCubesAlgorithm;
class myQuaternion;
class GrainHull;

using namespace std;

struct GBInfo {
	double energy;
	double mobility;
	GBInfo() {
	}
	GBInfo(double _mobility, double _energy) :
		energy(_energy), mobility(_mobility) {
	}
	GBInfo(const GBInfo& other) :
		mobility(other.mobility), energy(other.energy) {
	}
};

class InterfacialElement {
protected:
	vector<Triangle> m_Triangles;
	vector<Vector3d> m_barycenterTriangles;
	vector<Vector3d> m_UnitNormalTriangles;
	int m_Key_NeighborList;
	double m_mobility;
	double m_energy;
	GrainHull* m_owner;
public:
	InterfacialElement();
	InterfacialElement(int key, GrainHull* owner);
	virtual ~InterfacialElement();
	double computeMobilityMisori(double misori);
	double computeReadShockleyEnergy(double misori);

	virtual void computeEnergy()= 0;
	virtual void computeMobility()= 0;
	virtual Vector3d get_Position()=0;

	void addBaryCenter(Vector3d current) {
		m_barycenterTriangles.push_back(current);
	}
	void addTriangle(Triangle current) {
		m_Triangles.push_back(current);
	}
	inline int get_m_Key_NeighborList() {
		return m_Key_NeighborList;
	}
	inline GBInfo get_GBInfo() {
		return GBInfo(m_energy, m_mobility);
	}
	inline double get_energy() {
		return m_energy;
	}
	inline double get_mobility() {
		return m_mobility;
	}
};
class QuadrupleJunction;

class HighOrderJunction: public InterfacialElement {
	Vector3d m_position;
	vector<int> m_neighborIDs;
	//TODO:
public:
	friend class GrainHull;
	HighOrderJunction(int key, GrainHull *owner);
	HighOrderJunction(QuadrupleJunction* A, QuadrupleJunction* B, GrainHull *owner);
	~HighOrderJunction();
	void computeEnergy();
	void computeMobility();
	void computePosition();
	void mergeWith(QuadrupleJunction* B);
	void mergeWith(HighOrderJunction* B);
	inline vector<int> get_NeighborIDs() {
		return m_neighborIDs;
	}
	;
	inline Vector3d get_Position() {
		return m_position;
	}
	;
};

class QuadrupleJunction: public InterfacialElement {
	Vector3d m_position;
	int m_neighborID[3];
public:
	friend class GrainHull;
	QuadrupleJunction(int key, GrainHull *owner);
	~QuadrupleJunction();
	void computeEnergy();
	void computeMobility();
	void computePosition();
	inline int get_FirstNeighbor() {
		return m_neighborID[0];
	}
	;
	inline int get_SecondNeighbor() {
		return m_neighborID[1];
	}
	;
	inline int get_ThirdNeighbor() {
		return m_neighborID[2];
	}
	;
	inline Vector3d get_Position() {
		return m_position;
	}
	;
};

class TripleLine: public InterfacialElement {
	vector<InterfacialElement*> m_vertices;
	int m_neighborID[2];
public:
	friend class GrainHull;
	TripleLine(int key, GrainHull *owner);
	TripleLine(int neighbor1, int neighbor2, GrainHull* owner);
	~TripleLine();
	void computeEnergy();
	void computeMobility();
	void findAdjacentJunctions(vector<QuadrupleJunction*> ,
			vector<HighOrderJunction*> );
	inline int get_FirstNeighbor() {
		return m_neighborID[0];
	}
	inline int get_SecondNeighbor() {
		return m_neighborID[1];
	}
	inline vector<InterfacialElement*> get_vertices() {
		return m_vertices;
	}
	inline Vector3d get_Position() {
		return Vector3d(-1,-1,-1);;
	}
};

class GrainBoundary: public InterfacialElement {
	vector<TripleLine*> m_edges; // saves the indexes of edges in clockwise order
	int m_neighborID;
	Vector3d inclination;
public:
	friend class GrainHull;
	GrainBoundary(int key, GrainHull *owner);
	~GrainBoundary();
	void computeEnergy();
	void computeMobility();
	void findAdjacentTripleLines(vector<TripleLine*> );
	inline vector<TripleLine*> get_edges() {
		return m_edges;
	}
	inline Vector3d get_Position() {
		return Vector3d(-1,-1,-1);
	}
};

#endif
