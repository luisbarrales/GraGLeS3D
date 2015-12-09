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

class InterfacialElement {
protected:
	vector<Triangle> m_Triangles;
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

	void addTriangle(Triangle current) {
		m_Triangles.push_back(current);
	}
	int get_m_Key_NeighborList() {
		return m_Key_NeighborList;
	}
	inline double get_Correction_Weight() {
		return m_energy * m_mobility;
	}
	inline double get_energy() {
		return m_energy;
	}
	inline double get_mobility() {
		return m_mobility;
	}
};
class HighOrderJunction: public InterfacialElement {
	//TODO:
public:
	friend class GrainHull;
	HighOrderJunction(int key, GrainHull *owner);
	~HighOrderJunction();
	void computeEnergy();
	void computeMobility();
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
};

class TripleLine: public InterfacialElement {
	vector<int> m_vertices;
	int m_neighborID[2];
public:
	friend class GrainHull;
	TripleLine(int key, GrainHull *owner);
	TripleLine(int neighbor1, int neighbor2, GrainHull* owner);
	~TripleLine();
	void computeEnergy();
	void computeMobility();
};

class GrainBoundary: public InterfacialElement {
	vector<int> m_edges; // saves the indexes of edges in clockwise order
	int m_neighborID;
public:
	friend class GrainHull;
	GrainBoundary(int key, GrainHull *owner);
	~GrainBoundary();
	void computeEnergy();
	void computeMobility();
};

#endif
