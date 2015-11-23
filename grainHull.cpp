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

#include "grainHull.h"
#include "box.h"
#include "marchingCubes.h"
#include <stdexcept>
#include <map>

using namespace std;


QuadrupleJunction::QuadrupleJunction(){
}
QuadrupleJunction::~QuadrupleJunction(){
}


TripleLine::TripleLine(){
}
TripleLine::~TripleLine(){
}


GrainBoundary::GrainBoundary(int id): neighbor(id){
}
GrainBoundary::~GrainBoundary(){
}

void GrainBoundary::addTriangle(Triangle current){
	GBTriangles.push_back(current);
}

GrainHull::GrainHull(LSbox* owner) :
	m_owner(owner) {
}
GrainHull::~GrainHull() {
}

double GrainHull::computeVolume() {
	double volume = 0;
	for (unsigned int i = 0; i < m_actualHull.size(); i++) {
		Triangle& tri = m_actualHull[i];
		double v321 = tri.points[2][0] * tri.points[1][1] * tri.points[0][2];
		double v231 = tri.points[1][0] * tri.points[2][1] * tri.points[0][2];
		double v312 = tri.points[2][0] * tri.points[0][1] * tri.points[1][2];
		double v132 = tri.points[0][0] * tri.points[2][1] * tri.points[1][2];
		double v213 = tri.points[1][0] * tri.points[0][1] * tri.points[2][2];
		double v123 = tri.points[0][0] * tri.points[1][1] * tri.points[2][2];
		volume += (1.0f / 6.0f) * (-v321 + v231 + v312 - v132 - v213 + v123);
	}
	double h = m_owner->get_h();
	volume = abs(volume) * h * h * h;
	return volume;
}

double GrainHull::computeSurface() {
	double surface = 0;
	for (unsigned int i = 0; i < m_actualHull.size(); i++) {
		Triangle& tri = m_actualHull[i];
		Vector3d AB = tri.points[0] - tri.points[1];
		Vector3d BC = tri.points[0] - tri.points[2];
		surface += AB.cross(BC).norm() / 2.0;
	}
	double h = m_owner->get_h();
	surface = surface * h * h;
	return surface;
}

const NeighborList& GrainHull::getNeighborList(const Triangle& triangle) {
	if (triangle.additionalData < 0 || triangle.additionalData
			>= m_triangleNeighborLists.size()) {
		throw runtime_error(
				"Invalid additional data in triangle. Neighbor List unavailable!");
	}

	return m_triangleNeighborLists[triangle.additionalData];
}

const vector<unsigned int>& GrainHull::getAllNeighbors() {
	return m_neighbors;
}

bool GrainHull::generateHull() {
	MarchingCubesAlgorithm marcher(m_owner->getInputDistanceBuffer(), m_owner);
	marcher.generateHull(m_actualHull, m_triangleNeighborLists);
	m_neighbors = marcher.getIdentifiedNeighbors();

	if (m_actualHull.size() == 0)
		return false;
	else
		return true;
}

double pointToTriangleDistance(Vector3d& point, Triangle& triangle);

const Triangle& GrainHull::projectPointToSurface(Vector3d& point) {
	double minimalDistance = 10000000.0;
	unsigned int minIndex = 0xFFFFFFFF;
	for (unsigned int i = 0; i < m_actualHull.size(); i++) {
		double distance = pointToTriangleDistance(point, m_actualHull[i]);
		if (distance < minimalDistance) {
			minimalDistance = distance;
			minIndex = i;
		}
	}
	return m_actualHull[minIndex];
}

void GrainHull::computeGrainBoudaryElements() {
	m_Grainboundary.clear();
	for (unsigned int i = 0; i < m_actualHull.size(); i++) {
		switch (m_actualHull[i].adjacentGrainIDs.size())case 1: {
			recordGB();
			//triangle has only one adjacent grain
			break;
		}
		case 2: {
			//triangle is part of tripleLine
			break;
		}
		case 3: {
			//triangle contains to QuadrupleJunction
			break;
		}
		default: {
			//probably error case
		}

	}
}

void GrainHull::recordGB(Triangle current){
	bool found = false;
	int id = current.adjacentGrainIDs[0];
	for(const auto it : m_Grainboundary)
	{
		if(it->neighborID == id)
		{
			found = true;
			break;
		}
	}
	if (!found){
		m_Grainboundary.push_back(GrainBoundary(id));
		m_Grainboundary.end()->addTriangle(current);
		m_Grainboundary.end()->computeGrainBoundaryProperties();
	}
	else {
		// add triangle too local list
		it->addTriangle(current);
	}
}


const Triangle& GrainHull::projectPointGrainBoundary(Vector3d& point,
		GrainBoundary* nearestPlane) {
	double minimalDistance = 10000000.0;
	unsigned int minIndex = 0xFFFFFFFF;
	for (unsigned int i = 0; i < m_actualHull.size(); i++) {
		double distance = pointToTriangleDistance(point, m_actualHull[i]);
		if (distance < minimalDistance) {
			minimalDistance = distance;
			minIndex = i;
		}
	}
	return m_actualHull[minIndex];
}

double pointToTriangleDistance(Vector3d& point, Triangle& triangle) {
	//From Real-Time Collision Detection by Christer Ericson, published by Morgan Kaufmann Publishers,  2005 Elsevier Inc
	// Check if point in vertex region outside triangle.points[0]
	Vector3d ab = triangle.points[1] - triangle.points[0];
	Vector3d ac = triangle.points[2] - triangle.points[0];
	Vector3d ap = point - triangle.points[0];
	double d1 = ab.dot(ap);
	double d2 = ac.dot(ap);
	if (d1 <= 0.0f && d2 <= 0.0f)
		return (triangle.points[0] - point).norm(); // barycentric coordinates (1,0,0)
	// Check if point in vertex region outside triangle.points[1]
	Vector3d bp = point - triangle.points[1];
	double d3 = ab.dot(bp);
	double d4 = ac.dot(bp);
	if (d3 >= 0.0f && d4 <= d3)
		return (triangle.points[1] - point).norm(); // barycentric coordinates (0,1,0)
	// Check if point in edge region of AB, if so return projection of point onto AB
	double vc = d1 * d4 - d3 * d2;
	if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
		double v = d1 / (d1 - d3);
		return ((triangle.points[0] + ab * v) - point).norm(); // barycentric coordinates (1-v,v,0)
	}
	// Check if point in vertex region outside triangle.points[2]
	Vector3d cp = point - triangle.points[2];
	double d5 = ab.dot(cp);
	double d6 = ac.dot(cp);
	if (d6 >= 0.0f && d5 <= d6)
		return (triangle.points[2] - point).norm(); // barycentric coordinates (0,0,1)
	// Check if point in edge region of AC, if so return projection of point onto AC
	double vb = d5 * d2 - d1 * d6;
	if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
		double w = d2 / (d2 - d6);
		return ((triangle.points[0] + ac * w) - point).norm(); // barycentric coordinates (1-w,0,w)
	}
	// Check if point in edge region of BC, if so return projection of point onto BC
	double va = d3 * d6 - d5 * d4;
	if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
		double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		return ((triangle.points[1] + (triangle.points[2] - triangle.points[1]))
				* w - point).norm(); // barycentric coordinates (0,1-w,w)
	}
	// point inside face region. Compute Q through its barycentric coordinates (u,v,w)
	double denom = 1.0f / (va + vb + vc);
	double v = vb * denom;
	double w = vc * denom;
	return ((triangle.points[0] + ab * v + ac * w) - point).norm(); // = u*triangle.points[0] + v*triangle.points[1] + w*triangle.points[2], u = va * denom = 1.0f - v - w
}

struct vectorComparator {
	bool operator()(const Vector3d& a, const Vector3d& b) const {
		return a[0] < b[0] ? true : (a[0] > b[0] ? false : (a[1] < b[1] ? true
				: (a[1] > b[1] ? false : (a[2] < b[2] ? true : false))));
	}
};

void GrainHull::plotContour(bool absoluteCoordinates, int timestep) {
	string filename = string("GrainHull_") + to_string(
			(unsigned long long) m_owner->getID()) + string("Timestep_")
			+ to_string((unsigned long long) timestep) + string(".vtk");
	FILE* output = fopen(filename.c_str(), "wt");
	if (output == NULL) {
		throw runtime_error("Unable to save box hull!");
	}

	fprintf(output, "%s\n", "# vtk DataFile Version 3.0\n"
		"vtk output\n"
		"ASCII\n"
		"DATASET POLYDATA\n");

	int counter = 0;
	map<Vector3d, int, vectorComparator> mymap;
	map<int, Vector3d> orderedPoints;
	for (unsigned int i = 0; i < m_actualHull.size(); i++) {
		if (mymap.find(m_actualHull[i].points[0]) == mymap.end()) {
			mymap.insert(
					pair<Vector3d, int> (m_actualHull[i].points[0], counter));
			counter++;
		}
		if (mymap.find(m_actualHull[i].points[1]) == mymap.end()) {
			mymap.insert(
					pair<Vector3d, int> (m_actualHull[i].points[1], counter));
			counter++;
		}
		if (mymap.find(m_actualHull[i].points[2]) == mymap.end()) {
			mymap.insert(
					pair<Vector3d, int> (m_actualHull[i].points[2], counter));
			counter++;
		}
	}for ( const auto &myPair : mymap )
	{
		orderedPoints.insert(pair<int, Vector3d>(myPair.second, myPair.first));
	}

	fprintf(output, "POINTS %lu float\n", orderedPoints.size());

	for ( const auto &myPair : orderedPoints )
	{
		fprintf(output, "%f %f %f\n", myPair.second[0], myPair.second[1], myPair.second[2]);
	}

	fprintf(output, "POLYGONS %lu %lu\n", m_actualHull.size(),
			m_actualHull.size() * 4);
	for (unsigned int i = 0; i < m_actualHull.size(); i++) {

		fprintf(output, "3 %d %d %d \n",
				(*(mymap.find(m_actualHull[i].points[2]))).second,
				(*(mymap.find(m_actualHull[i].points[1]))).second,
				(*(mymap.find(m_actualHull[i].points[0]))).second);
	}

	fprintf(output, "POINT_DATA %lu\n", orderedPoints.size());
	fprintf(output, "FIELD FieldData 1\n");
	fprintf(output, "Number of adjacent grains 1 %lu int\n",
			orderedPoints.size());

	for ( const auto &myPair : orderedPoints )
	{
		const Vector3d& point = myPair.second;
		int interestingness = 0;
		for(unsigned int i=0; i<m_actualHull.size(); i++)
		{
			if( point == m_actualHull[i].points[0] || point == m_actualHull[i].points[1] ||
					point == m_actualHull[i].points[2])
			{
				const NeighborList& list = m_triangleNeighborLists[m_actualHull[i].additionalData];
				int interactingGrains=0;
				for(int j=0; j<NEIGHBOR_LIST_SIZE; j++)
				interactingGrains += (list.neighbors[j] == 0xFFFFFFFF ? 0 : 1);
				interestingness = max(interestingness, interactingGrains);
			}
		}
		fprintf(output, "%d ", interestingness);
	}
	fclose(output);
}
