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
#include "InterfacialElement.h"

using namespace std;

class LSbox;

class GrainHull {
private:
	vector<Triangle> m_actualHull;
	vector<NeighborList> m_triangleNeighborLists;
	vector<unsigned int> m_neighbors;
	vector<GrainBoundary*> m_Grainboundary;
	vector<TripleLine*> m_TripleLines;
	vector<QuadrupleJunction*> m_QuadruplePoints;
public:
	friend class GrainBoundary;
	friend class TripleLine;
	friend class QuadrupleJunction;
	LSbox* m_owner;
	GrainHull(LSbox* owner);
	~GrainHull();
	bool generateHull();
	const NeighborList& getNeighborList(const Triangle& triangle);
	double computeGrainVolume();
	double computeSurfaceArea();
	const Triangle& projectPointToSurface(Vector3d& point);
	const vector<unsigned int>& getAllNeighbors();
	unsigned int getAllNeighborsCount() {
		return m_neighbors.size();
	}
	void plotContour(bool absoluteCoordinates, int timestep);
	QuadrupleJunction* findQuadrupleJunction(int key);
	TripleLine* findTripleLine(int key);
	GrainBoundary* findGrainBoundary(int key);
	//new:
	void clearInterfacialElements();
	void computeGrainBoundaryElements();
	void subDivideTrianglesToInterfacialElements();
	void computeInterfacialElementMesh();
	double projectPointToGrainBoundary(Vector3d& point, int id);
	void plotInterfacialElements(bool absoluteCoordinates, int timestep);
};

#endif //__GRAIN_HULL__
