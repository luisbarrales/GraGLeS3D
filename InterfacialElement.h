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
	InterfacialElement(int key, GrainHull owner);
	virtual ~InterfacialElement();
	double computeMobility(double misori);
	double computeReadShockleyEnergy(double misori);
	virtual void computeEnergy();
	virtual void computeMobility();
	void addTriangle(Triangle current) {
		m_Triangles.push_back(current);
	}
	int get_m_Key_NeighborList() {
		return m_Key_NeighborList;
	}
};

class QuadrupleJunction: public InterfacialElement {
	Triangle QuadruplePointTriangles[2];
	Vector3d position;
	int neighborID[3];
public:
	QuadrupleJunction(int key, GrainHull owner);
	~QuadrupleJunction();
};

class TripleLine: public InterfacialElement {
	vector<int> vertices;
	int neighborID[2];
public:
	TripleLine(int key, GrainHull owner);
	~TripleLine();
};

class GrainBoundary: public InterfacialElement {
	vector<int> edges; // saves the indexes of edges in clockwise order
	int neighborID;
public:
	GrainBoundary(int key, GrainHull owner);
	~GrainBoundary();
};

#endif
