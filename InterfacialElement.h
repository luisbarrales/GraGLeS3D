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

using namespace std;

class InterfacialElement {
protected:
	vector<Triangle> Triangles;
	int Key_NeighborList;
	double mobility;
	double energy;
public:
	InterfacialElement();
	InterfacialElement(int key, NeighborList newNeighborList);
	virtual ~InterfacialElement();
	double computeMobility(double misori);
	double computeReadShockleyEnergy(double misori);
	void addTriangle(Triangle current) {
		Triangles.push_back(current);
	}
};

class QuadrupleJunction: public InterfacialElement {
	Triangle QuadruplePointTriangles[2];
	Vector3d position;
	int neighborID[3];
public:
	QuadrupleJunction(int key, NeighborList newQuadrupleJunction, LSbox* owner);
	~QuadrupleJunction();
};

class TripleLine: public InterfacialElement {
	vector<int> vertices;
	int neighborID[2];
public:
	TripleLine(int key, NeighborList newTripleLine, LSbox* owner);
	~TripleLine();
};

class GrainBoundary: public InterfacialElement {
	vector<int> edges; // saves the indexes of edges in clockwise order
	int neighborID;
public:
	GrainBoundary(int key, NeighborList newboundary, LSbox* myBox);
	~GrainBoundary();
};

#endif
