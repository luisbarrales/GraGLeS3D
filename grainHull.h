/*
	GraGLeS 2D A grain growth simulation utilizing level set approaches
    Copyright (C) 2015  Christian Miessen, Nikola Velinov

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef __GRAIN_HULL__
#define __GRAIN_HULL__

#include "triangle.h"
#include "marchingCubes.h"
#include "Eigen/Dense"
#include <vector>

using namespace std;

class QuadrupleJunction{
private:
	Triangle QuadruplePointTriangles[2];
	Vector3d position;
	int neighborID[3];
	double mobility;
	double weight;
public:
	QuadrupleJunction();
	~QuadrupleJunction();

};
class TripleLine{
private:
	vector<QuadrupleJunction> vertices;
	vector<Triangle> TripleLineTriangles;
	int neighborID[2];
	double energy;
	double mobility;
public:
	TripleLine();
	~TripleLine();
};

class GrainBoundary{
private:
	vector<Triangle> GBTriangles;
	vector<TripleLines> edges;
	vector<QuadrupleJunction> vertices;
	int neighborID;
	double mobility;
	double energy;
public:
	GrainBoundary(int id);
	~GrainBoundary();
	//TODO:
	void addTriangle(Triangle);
	void computeGrainBoundaryProperties();
};

class LSbox;

class GrainHull
{
private:
	vector<Triangle>			m_actualHull;
	vector<NeighborList>		m_triangleNeighborLists;
	LSbox*						m_owner;
	vector<unsigned int>		m_neighbors;
	vector<GrainBoundary>		m_Grainboundary;
public:
	GrainHull(LSbox* owner);
	~GrainHull();
	bool 						generateHull();
	const NeighborList&			getNeighborList(const Triangle& triangle);
	double 						computeVolume();
	double 						computeSurface();
	const Triangle&				projectPointToSurface(Vector3d& point);
	const vector<unsigned int>&	getAllNeighbors();
	unsigned int				getAllNeighborsCount() { return m_neighbors.size(); }
	void 						plotContour(bool absoluteCoordinates, int timestep);
	//new:
	void 						computeGrainBoudaryElements();
	void 						recordGB(Triangle current);
	void 						computeGrainBoundaryProperties();
	const Triangle&				projectPointGrainBoundary(Vector3d& point, GrainBoundary* nearestPlane);
};

#endif //__GRAIN_HULL__
