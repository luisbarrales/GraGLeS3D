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
#include "triangle.h"
#include "grainHull.h"
#include "box.h"
#include "marchingCubes.h"
#include <stdexcept>
#include <map>
#include <vector>
#include "spoint.h"
#include "grahamScan.h"
#include "grainhdl.h"
#include <fstream>
#include "Structs.h"

using namespace std;

GrainHull::GrainHull(LSbox* owner) :
		m_owner(owner) {
}
GrainHull::~GrainHull() {
}

double GrainHull::computeGrainVolume() {
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

double GrainHull::computeSurfaceArea() {
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
	if (triangle.additionalData < 0
			|| triangle.additionalData >= m_triangleNeighborLists.size()) {
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

void GrainHull::clearInterfacialElements() {
	for (auto it : m_Grainboundary)
		delete &(*it);
	m_Grainboundary.clear();
	for (auto it : m_TripleLines)
		delete &(*it);
	m_TripleLines.clear();
	for (auto it : m_QuadruplePoints)
		delete &(*it);
	m_QuadruplePoints.clear();
	for (auto it : m_HighOrderJunctions)
		delete &(*it);
	m_HighOrderJunctions.clear();
}

void GrainHull::computeGrainBoundaryElements() {
	clearInterfacialElements();
	for (unsigned int i = 0; i < m_triangleNeighborLists.size(); i++) {
		int junctionType = m_triangleNeighborLists[i].getNeighborsListCount();
		switch (junctionType) {
		case 2: {
			GrainBoundary* newGB = new GrainBoundary(i, this);
			m_Grainboundary.push_back(newGB);
			//triangle has only one adjacent grain
			break;
		}
		case 3: {
			TripleLine* newTL = new TripleLine(i, this);
			m_TripleLines.push_back(newTL);
			//triangle is part of tripleLine
			break;
		}
		case 4: {
			QuadrupleJunction* newQJ = new QuadrupleJunction(i, this);
			m_QuadruplePoints.push_back(newQJ);
			//triangle contains to QuadrupleJunction
			break;
		}
		default: {
			HighOrderJunction* newHJ = new HighOrderJunction(i, this);
			m_HighOrderJunctions.push_back(newHJ);
			// high order junction is found
			//TODO:
			break;
		}
		}
	}
}

void GrainHull::subDivideTrianglesToInterfacialElements() {
	for (unsigned int i = 0; i < m_actualHull.size(); i++) {
		int key = m_actualHull[i].additionalData;
		unsigned int type =
				m_triangleNeighborLists[key].getNeighborsListCount();
		switch (type) {
		case 2: {
			GrainBoundary* GB = findGrainBoundary(key);
			Vector3d current = m_actualHull[i].computeBarycenter();
			GB->addBaryCenter(current);
			GB->addTriangle(m_actualHull[i]);
			break;
		}
		case 3: {
			TripleLine* TL = findTripleLine(key);
			Vector3d current = m_actualHull[i].computeBarycenter();
			TL->addBaryCenter(current);
			TL->addTriangle(m_actualHull[i]);
			break;
		}
		case 4: {
			QuadrupleJunction* QJ = findQuadrupleJunction(key);
			Vector3d current = m_actualHull[i].computeBarycenter();
			QJ->addBaryCenter(current);
			QJ->addTriangle(m_actualHull[i]);
			break;
		}
		default: {
			HighOrderJunction* HJ = findHighOrderJunction(key);
			Vector3d current = m_actualHull[i].computeBarycenter();
			HJ->addBaryCenter(current);
			//cout << "high order junction found" << endl;
			break;
		}
		}
	}
	//m_actualHull.clear();
}

void GrainHull::computeJunctionPosition() {
	for (const auto it : m_QuadruplePoints) {
		it->computePosition();
	}
	for (const auto it : m_HighOrderJunctions) {
		it->computePosition();
	}
	mergeJunction();
}

Vector3d GrainHull::findClosestJunctionTo(Vector3d myposition) {
	double distanceH = 1000000000;
	double distanceQ = 1000000000;
	double distanceTo;
	int posQuadruple;
	int posHigherOrder;

	for (int i = 0; i < m_QuadruplePoints.size(); i++) {
		distanceTo =
				(m_QuadruplePoints[i]->get_Position() - myposition).squaredNorm();
		if (distanceTo < distanceQ) {
			distanceQ = distanceTo;
			posQuadruple = i;
		}
	}

	for (int i = 0; i < m_HighOrderJunctions.size(); i++) {
		distanceTo =
				(m_HighOrderJunctions[i]->get_Position() - myposition).squaredNorm();
		if (distanceTo < distanceH) {
			distanceH = distanceTo;
			posHigherOrder = i;
		}
	}
	if (distanceH > distanceQ) {
		return m_QuadruplePoints[posQuadruple]->get_Position();
	} else {
		return m_HighOrderJunctions[posHigherOrder]->get_Position();
	}
}

void GrainHull::correctJunctionPositionWithNeighborInformation() {
	/*
	 * some debugging
	 */

	vector<double> distances;
	double h = m_owner->get_grainHandler()->get_h();
	string filenamePos;
	filenamePos = "Position_";
	filenamePos += to_string(m_owner->getID());
	filenamePos += "Timestep_";
	filenamePos += to_string(m_owner->get_grainHandler()->get_loop());
	ofstream CorrectingPositions;
//	CorrectingPositions.open(filenamePos.c_str());
//
//	CorrectingPositions << "xcoord ycoord zcoord scalar" << endl;

	/*
	 * end of debugging
	 */

	for (const auto it : m_QuadruplePoints) {
		vector<int> neighbors = it->get_NeighborIDs();
		Vector3d correspondingJunctions(it->get_Position());
		int N = 0;
		for (int i = 0; i < neighbors.size(); i++) {
			if (neighbors[i] != 0) {

				/*
				 * some debugging
				 */
				Vector3d posClosestJunction;
				posClosestJunction = m_owner->get_grainHandler()->getGrainByID(
						neighbors[i])->findClosestJunctionTo(
						it->get_Position());
				distances.push_back(
						(posClosestJunction - it->get_Position()).norm());

				if (m_owner->getID() == 1)
					if (distances.back() > 7) {
						CorrectingPositions << (it->get_Position())(0) << " "
								<< (it->get_Position())(1) << " "
								<< (it->get_Position())(2) << " "
								<< distances.back() << endl;
						CorrectingPositions << posClosestJunction(0) << " "
								<< posClosestJunction(1) << " "
								<< posClosestJunction(2) << " "
								<< distances.back() << endl;

						cout << (it->get_Position())(0) << " "
								<< (it->get_Position())(1) << " "
								<< (it->get_Position())(2) << " "
								<< distances.back() << endl;
						cout << posClosestJunction(0) << " "
								<< posClosestJunction(1) << " "
								<< posClosestJunction(2) << " "
								<< distances.back() << endl;
					}
				/*
				 * end of debugging
				 */
				if (distances.back() < 3 * h) {
					correspondingJunctions += posClosestJunction;
					N++;
				}

			}
		}
		correspondingJunctions /= double(N + 1);
		it->set_BufferPosition(correspondingJunctions);
	}
	for (const auto it : m_HighOrderJunctions) {
		vector<int> neighbors = it->get_NeighborIDs();
		Vector3d correspondingJunctions(it->get_Position());
		int N = 0;
		for (int i = 0; i < neighbors.size(); i++) {
			if (neighbors[i] != 0) {
				/*
				 * some debugging
				 */
				Vector3d posClosestJunction;
				posClosestJunction = m_owner->get_grainHandler()->getGrainByID(
						neighbors[i])->findClosestJunctionTo(
						it->get_Position());
				distances.push_back(
						(posClosestJunction - it->get_Position()).norm());

				if (m_owner->getID() == 1)
					if (distances.back() > 7) {
						CorrectingPositions << (it->get_Position())(0) << " "
								<< (it->get_Position())(1) << " "
								<< (it->get_Position())(2) << " "
								<< distances.back() << endl;
						CorrectingPositions << posClosestJunction(0) << " "
								<< posClosestJunction(1) << " "
								<< posClosestJunction(2) << " "
								<< distances.back() << endl;

						cout << (it->get_Position())(0) << " "
								<< (it->get_Position())(1) << " "
								<< (it->get_Position())(2) << " "
								<< distances.back() << endl;
						cout << posClosestJunction(0) << " "
								<< posClosestJunction(1) << " "
								<< posClosestJunction(2) << " "
								<< distances.back() << endl;
					}

				/*
				 * end of debugging
				 */
				if (distances.back() < 3 * h) {
					correspondingJunctions += posClosestJunction;
					N++;
				}

			}
		}
		correspondingJunctions /= (double) (N + 1);
		it->set_BufferPosition(correspondingJunctions);
	}

	/*
	 * some debugging
	 */

//	if (m_owner->getID() == 1)
//		CorrectingPositions.close();
//	if (m_owner->getID() == 1) {
//		string filename;
//		filename = "DistanceHistogram_Grain_";
//		filename += to_string(m_owner->getID());
//		filename += "Timestep_";
//		filename += to_string(m_owner->get_grainHandler()->get_loop());
//		ofstream DistanceHistogram;
//		DistanceHistogram.open(filename.c_str());
//		for (int i = 0; i < distances.size(); i++) {
////			DistanceHistogram << distances[i] << endl;
//		}
//		DistanceHistogram.close();
//	}
	/*
	 * end of debugging
	 */

}

void GrainHull::switchBufferPositions() {

	for (const auto it : m_QuadruplePoints)
		it->switch_BufferPosition();
	for (const auto it : m_HighOrderJunctions)
		it->switch_BufferPosition();
}

void GrainHull::computeInterfacialElementMesh() {

	// check for infinitisimal short Triplelines, which could cause stability problmes in computation

	for (const auto it : m_TripleLines) {
		//TODO:
		//is it possible to merge Triplelines with less than 1 vertice?
		it->findAdjacentJunctions(m_QuadruplePoints, m_HighOrderJunctions);
	}
	for (const auto it : m_Grainboundary) {
		it->findAdjacentTripleLines(m_TripleLines);
	}
	if (m_owner->get_grainHandler()->get_loop() == 150
			|| m_owner->get_grainHandler()->get_loop() == 100
			|| m_owner->get_grainHandler()->get_loop() == 250
			|| m_owner->get_grainHandler()->get_loop() == 500) {
		meanWidth();
		computeTriplelineLength();
	}
}
void GrainHull::mergeJunction() {
	for (int i = 0; i < m_QuadruplePoints.size(); i++) {
		for (int j = i + 1; j < m_QuadruplePoints.size(); j++) {
			if ((m_QuadruplePoints[i]->get_Position()
					- m_QuadruplePoints[j]->get_Position()).norm() < 2) {
				//create new high order junction / delete old
				HighOrderJunction* newHJ = new HighOrderJunction(
						m_QuadruplePoints[i], m_QuadruplePoints[j], this);
				m_HighOrderJunctions.push_back(newHJ);
				delete m_QuadruplePoints[i];
				delete m_QuadruplePoints[j];
				m_QuadruplePoints.erase(m_QuadruplePoints.begin() + i);
				m_QuadruplePoints.erase(m_QuadruplePoints.begin() + j - 1);
				i--;
				j--;
				break;
			}
		}
	}

	for (int k = 0; k < m_HighOrderJunctions.size(); k++) {
		for (int i = 0; i < m_QuadruplePoints.size(); i++) {
//			cout << "Quadruple point position" << endl;
//			cout << m_QuadruplePoints[i]->get_Position() << endl;
//			cout << "High order junction position" << endl;
//			cout << m_HighOrderJunctions[k]->get_Position() << endl;
			if ((m_QuadruplePoints[i]->get_Position()
					- m_HighOrderJunctions[k]->get_Position()).norm() < 2) {
				// add Quadruple Junction to existing high order junction
				m_HighOrderJunctions[k]->mergeWith(m_QuadruplePoints[i]);
				delete m_QuadruplePoints[i];
				m_QuadruplePoints.erase(m_QuadruplePoints.begin() + i);
				i--;
			}
		}
	}
	//Changed Quadruple Points to HighOrderJunctions
	for (int i = 0; i < m_HighOrderJunctions.size(); i++) {
		for (int k = i + 1; k < m_HighOrderJunctions.size(); k++) {
			if ((m_HighOrderJunctions[i]->get_Position()
					- m_HighOrderJunctions[k]->get_Position()).norm() < 2) {
				// add Quadruple Junction to existing high order junction
				m_HighOrderJunctions[i]->mergeWith((m_HighOrderJunctions[k]));
				delete m_HighOrderJunctions[k];
				m_HighOrderJunctions.erase(m_HighOrderJunctions.begin() + k);
				k--;
			}
		}
	}
}
GrainBoundary* GrainHull::findGrainBoundary(int key) {
	for (const auto it : m_Grainboundary) {
		if (it->get_m_Key_NeighborList() == key) {
			return &(*it);
		}
	}
	return NULL;
}

TripleLine* GrainHull::findTripleLine(int key) {
	for (const auto it : m_TripleLines) {
		if (it->get_m_Key_NeighborList() == key) {
			return &(*it);
		}
	}
	return NULL;
}

QuadrupleJunction* GrainHull::findQuadrupleJunction(int key) {

	for (const auto it : m_QuadruplePoints) {
		if (it->get_m_Key_NeighborList() == key) {
			return &(*it);
		}
	}
	return NULL;
}

HighOrderJunction* GrainHull::findHighOrderJunction(int key) {

	for (const auto it : m_HighOrderJunctions) {
		if (it->get_m_Key_NeighborList() == key) {
			return &(*it);
		}
	}
	return NULL;
}

GBInfo GrainHull::projectPointToGrainBoundary(Vector3d& point, int id) {
	double minimalDistance = 10000000.0;
	GBInfo weight(1., 1.);

	for (int j = 0; j < m_HighOrderJunctions.size(); j++) {
		bool found = false;
		for (int i = 0; i < m_HighOrderJunctions[j]->m_neighborIDs.size();
				i++) {
			if (m_HighOrderJunctions[j]->m_neighborIDs[i] == id)
				found = true;
		}
		if (found) {
			//			for (unsigned int i = 0; i
			//					< m_HighOrderJunctions[j]->m_Triangles.size(); i++) {
			for (unsigned int i = 0;
					i < m_HighOrderJunctions[j]->m_barycenterTriangles.size();
					i++) {
				//				double distance = pointToTriangleDistance(point,
				//						m_HighOrderJunctions[j]->m_Triangles[i]);

				double distance =
						(point
								- m_HighOrderJunctions[j]->m_barycenterTriangles[i]).squaredNorm();
				if (distance < minimalDistance) {
					minimalDistance = distance;
					weight = m_HighOrderJunctions[j]->get_GBInfo();
				}
			}
		}
	}

//search in QuadrupleJunctions
	for (int j = 0; j < m_QuadruplePoints.size(); j++) {
		if (m_QuadruplePoints[j]->m_neighborIDs[0] == id
				|| m_QuadruplePoints[j]->m_neighborIDs[1] == id
				|| m_QuadruplePoints[j]->m_neighborIDs[2] == id) {
			//			for (unsigned int i = 0;
			//					i < m_QuadruplePoints[j]->m_Triangles.size(); i++) {
			for (unsigned int i = 0;
					i < m_QuadruplePoints[j]->m_barycenterTriangles.size();
					i++) {
				//				double distance = pointToTriangleDistance(point,
				//						m_QuadruplePoints[j]->m_Triangles[i]);
				double distance =
						(point - m_QuadruplePoints[j]->m_barycenterTriangles[i]).squaredNorm();
				if (distance < minimalDistance) {
					minimalDistance = distance;
					weight = m_QuadruplePoints[j]->get_GBInfo();
				}
			}
		}
	}
//search in TripleJunctions
	for (int j = 0; j < m_TripleLines.size(); j++) {
		if (m_TripleLines[j]->m_neighborIDs[0] == id
				|| m_TripleLines[j]->m_neighborIDs[1] == id) {
			//			for (unsigned int i = 0; i < m_TripleLines[j]->m_Triangles.size();
			//					i++) {
			for (unsigned int i = 0;
					i < m_TripleLines[j]->m_barycenterTriangles.size(); i++) {
				//				double distance = pointToTriangleDistance(point,
				//						m_TripleLines[j]->m_Triangles[i]);
				double distance =
						(point - m_TripleLines[j]->m_barycenterTriangles[i]).squaredNorm();
				if (distance < minimalDistance) {
					minimalDistance = distance;
					weight = m_TripleLines[j]->get_GBInfo();
				}
			}
		}
	}
//search in GrainBoundaries:
	for (int j = 0; j < m_Grainboundary.size(); j++) {
		if (m_Grainboundary[j]->m_neighborIDs[0] == id) {
			//			for (unsigned int i = 0; i < m_Grainboundary[j]->m_Triangles.size();
			//					i++) {
			for (unsigned int i = 0;
					i < m_Grainboundary[j]->m_barycenterTriangles.size(); i++) {
				//				double distance = pointToTriangleDistance(point,
				//						m_Grainboundary[j]->m_Triangles[i]);
				double distance =
						(point - m_Grainboundary[j]->m_barycenterTriangles[i]).squaredNorm();
				if (distance < minimalDistance) {
					minimalDistance = distance;
					weight = m_Grainboundary[j]->get_GBInfo();
				}
			}
		}
	}
	return weight;
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
		return a[0] < b[0] ?
				true :
				(a[0] > b[0] ?
						false :
						(a[1] < b[1] ?
								true :
								(a[1] > b[1] ?
										false : (a[2] < b[2] ? true : false))));
	}
};

//void GrainHull::plotContour(bool absoluteCoordinates, int timestep) {
//	string filename = string("GrainHull_") + to_string(
//			(unsigned long long) m_owner->getID()) + string("Timestep_")
//			+ to_string((unsigned long long) timestep) + string(".vtk");
//	FILE* output = fopen(filename.c_str(), "wt");
//	if (output == NULL) {
//		throw runtime_error("Unable to save box hull!");
//	}
//
//	fprintf(output, "%s\n", "# vtk DataFile Version 3.0\n"
//		"vtk output\n"
//		"ASCII\n"
//		"DATASET POLYDATA\n");
//
//	int counter = 0;
//	map<Vector3d, int, vectorComparator> mymap;
//	map<int, Vector3d> orderedPoints;
//
//	for (unsigned int i = 0; i < m_actualHull.size(); i++) {
//		if (mymap.find(m_actualHull[i].points[0]) == mymap.end()) {
//			mymap.insert(
//					pair<Vector3d, int> (m_actualHull[i].points[0], counter));
//			counter++;
//		}
//		if (mymap.find(m_actualHull[i].points[1]) == mymap.end()) {
//			mymap.insert(
//					pair<Vector3d, int> (m_actualHull[i].points[1], counter));
//			counter++;
//		}
//		if (mymap.find(m_actualHull[i].points[2]) == mymap.end()) {
//			mymap.insert(
//					pair<Vector3d, int> (m_actualHull[i].points[2], counter));
//			counter++;
//		}
//	}for ( const auto &myPair : mymap )
//	{
//		orderedPoints.insert(pair<int, Vector3d>(myPair.second, myPair.first));
//	}
//
//	fprintf(output, "POINTS %lu float\n", orderedPoints.size());
//
//	for ( const auto &myPair : orderedPoints )
//	{
//		fprintf(output, "%f %f %f\n", myPair.second[0], myPair.second[1], myPair.second[2]);
//	}
//
//	fprintf(output, "POLYGONS %lu %lu\n", m_actualHull.size(),
//			m_actualHull.size() * 4);
//	for (unsigned int i = 0; i < m_actualHull.size(); i++) {
//		fprintf(output, "3 %d %d %d \n",
//				(*(mymap.find(m_actualHull[i].points[2]))).second,
//				(*(mymap.find(m_actualHull[i].points[1]))).second,
//				(*(mymap.find(m_actualHull[i].points[0]))).second);
//	}
//
//	fprintf(output, "POINT_DATA %lu\n", orderedPoints.size());
//	fprintf(output, "FIELD FieldData 1\n");
//	fprintf(output, "Interestingness 1 %lu int\n", orderedPoints.size());
//
//for ( const auto &myPair : orderedPoints )
//{
//	const Vector3d& point = myPair.second;
//	int interestingness = 0;
//	for(unsigned int i=0; i<m_actualHull.size(); i++)
//	{
//		if( point == m_actualHull[i].points[0] || point == m_actualHull[i].points[1] ||
//				point == m_actualHull[i].points[2])
//		{
//			//const NeighborList& list = m_triangleNeighborLists[m_actualHull[i].additionalData];
//			int interactingGrains = m_triangleNeighborLists[m_actualHull[i].additionalData].getNeighborsListCount();
//			//				for(int j=0; j<NEIGHBOR_LIST_SIZE; j++)
//			//				interactingGrains += (list.neighbors[j] == 0xFFFFFFFF ? 0 : 1);
//			interestingness = max(interestingness, interactingGrains);
//		}
//	}
//	fprintf(output, "%d ", interestingness);
//
//}
//fclose( output);
//}

void GrainHull::plotContour(bool absoluteCoordinates, int timestep) {
	string filename = string("GrainHull_")
			+ to_string((unsigned long long) m_owner->getID())
			+ string("Timestep_") + to_string((unsigned long long) timestep)
			+ string(".vtk");
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
					pair<Vector3d, int>(m_actualHull[i].points[0], counter));
			counter++;
		}
		if (mymap.find(m_actualHull[i].points[1]) == mymap.end()) {
			mymap.insert(
					pair<Vector3d, int>(m_actualHull[i].points[1], counter));
			counter++;
		}
		if (mymap.find(m_actualHull[i].points[2]) == mymap.end()) {
			mymap.insert(
					pair<Vector3d, int>(m_actualHull[i].points[2], counter));
			counter++;
		}
	}
	for (const auto &myPair : mymap) {
		orderedPoints.insert(pair<int, Vector3d>(myPair.second, myPair.first));
	}

	fprintf(output, "POINTS %lu float\n", orderedPoints.size());

	for (const auto &myPair : orderedPoints) {
		fprintf(output, "%f %f %f\n", myPair.second[0], myPair.second[1],
				myPair.second[2]);
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
	fprintf(output, "Interestingness 1 %lu int\n", orderedPoints.size());

	for (const auto &myPair : orderedPoints) {
		const Vector3d& point = myPair.second;
		int interestingness = 0;
		int key = 0;
		for (unsigned int i = 0; i < m_actualHull.size(); i++) {
			if (point == m_actualHull[i].points[0]
					|| point == m_actualHull[i].points[1]
					|| point == m_actualHull[i].points[2]) {
				//const NeighborList& list = m_triangleNeighborLists[m_actualHull[i].additionalData];
				int interactingGrains =
						m_triangleNeighborLists[m_actualHull[i].additionalData].getNeighborsListCount();
				key = m_actualHull[i].additionalData;
				//				for(int j=0; j<NEIGHBOR_LIST_SIZE; j++)
				//				interactingGrains += (list.neighbors[j] == 0xFFFFFFFF ? 0 : 1);
				interestingness = max(interestingness, interactingGrains);
			}
		}
		interestingness = 100 * interestingness + key;
		fprintf(output, "%d ", interestingness);

	}
	fclose(output);
}
void GrainHull::plotInterfacialElements(bool absoluteCoordinates,
		int timestep) {
	int ID = 0;
	string filename = string("InterfacialElements_")
			+ to_string((unsigned long long) m_owner->getID())
			+ string("Timestep_") + to_string((unsigned long long) timestep)
			+ string(".vtk");
	FILE* output = fopen(filename.c_str(), "wt");
	if (output == NULL) {
		throw runtime_error("Unable to save Interfacialelements hull!");
	}

	for (const auto it : m_QuadruplePoints) {
		for (const auto it2 : it->m_barycenterTriangles) {
			fprintf(output, "%lf \t %lf \t %lf \n ", it2[0], it2[1], it2[2]);
		}
		Vector3d position = it->get_Position();
		fprintf(output, "Position of Junction: %lf \t %lf \t %lf \n \n",
				position[0], position[1], position[2]);
	}
	for (const auto it : m_HighOrderJunctions) {
		for (const auto it2 : it->m_barycenterTriangles) {
			fprintf(output, "%lf \t %lf \t %lf \n ", it2[0], it2[1], it2[2]);
		}
		Vector3d position = it->get_Position();
		fprintf(output, "Position of Junction: %lf \t %lf \t %lf \n \n",
				position[0], position[1], position[2]);
	}
	fclose(output);
}

//fprintf(output, "\n\nGRAINBOUNDARY %lu\n", m_Grainboundary.size());
//
//for (const auto it : m_Grainboundary) {
//fprintf(output, "\n\nGRAINBOUNDARY \n");
//for (const auto it2 : it->m_barycenterTriangles) {
//	fprintf(output, "%lf \t %lf \t %lf \t %d \n ", it2[0], it2[1],
//			it2[2], ID);
//}
//ID++;
//fprintf(output, "\n\nEDGES %lu\n", it->m_edges.size());
//for (const auto it3 : it->m_edges) {
//	fprintf(output, "\n\nEDGES \n");
//	for (const auto it4 : it3->m_barycenterTriangles) {
//		fprintf(output, "%lf \t %lf \t %lf \t %d \n ", it4[0], it4[1],
//				it4[2], ID);
//	}
//	ID++;
//
//	for (const auto it5 : it3->m_vertices) {
//		fprintf(output, "\n\nVERTICES \n");
//		Vector3d point = it5->get_Position();
//		fprintf(output, "%lf \t %lf \t %lf \t %d \n ", point[0],
//				point[1], point[2], ID);
//		ID++;
//	}
//
//}
//
//}
//}

/*		fprintf(output, "\n\nTRIPLELINES %lu\n", m_TripleLines.size());
 for(const auto it : m_TripleLines)
 {
 for(const auto it3 : it->m_barycenterTriangles) {
 fprintf(output, "%lf \t %lf \t %lf \t %d \n ",it3[0],it3[1],it3[2], ID);
 }
 ID++;
 }

 fprintf(output, "\n\nQUADRUPLEJUNCTIONS %lu\n", m_QuadruplePoints.size());
 for(const auto it : m_QuadruplePoints)
 {
 for(const auto it4 : it->m_barycenterTriangles){
 fprintf(output, "%lf \t %lf \t %lf \t %d \n ", it4[0],it4[1],it4[2], ID);
 }
 ID++;
 }
 }
 */

void GrainHull::meanWidth() {
	m_normalVectors.resize(m_actualHull.size());
	Vector3d normal_temp;

	for (unsigned int i = 0; i < m_actualHull.size(); i++) {
		normal_temp =
				(m_actualHull[i].points[1] - m_actualHull[i].points[0]).cross(
						(m_actualHull[i].points[2] - m_actualHull[i].points[0]));
		if (normal_temp.norm() != 0)
			normal_temp.normalize();
		else {
			//			cout << "vector has length 0" << endl;
			//			cout << normal_temp.transpose() << endl;
			//			cout << (m_actualHull[i].points[1] - m_actualHull[i].points[0]).transpose() << endl;
			//			cout << (m_actualHull[i].points[2] - m_actualHull[i].points[0]).transpose() << endl;
		}
		m_normalVectors[i] = normal_temp;
	}
	/*
	 * Save for every triangle the IDs of the neighbor-triangles
	 */

	vector<int> NeighborList[m_actualHull.size()];

	for (unsigned int i = 0; i < m_actualHull.size(); i++) {
		if (m_normalVectors[i].norm() == 0)
			continue;
		for (int j = 0; j < i; j++) {
			if (m_actualHull[i].IsNeighbor(m_actualHull[j]))
				if (m_normalVectors[j].norm() != 0)
					NeighborList[i].push_back(j);
		}
	}
	/*
	 * Calculate the mean width of the grain
	 */

	m_LD = 0;
	for (unsigned int i = 0; i < m_actualHull.size(); i++) {
		vector<int>* neighbor = &NeighborList[i];
		//		if(neighbor->size()!=3)
		//			cout << "Triangle has less or more than 3 neighbors " << neighbor->size() << endl;
		for (unsigned int j = 0; j < neighbor->size(); j++) {
			m_LD += m_actualHull[i].calculateMeanWidthComponent(
					m_actualHull[(*neighbor)[j]], m_normalVectors[i],
					m_normalVectors[(*neighbor)[j]]);
		}
	}
	m_LD /= 2 * M_PI;
	if (m_LD < 0)
		m_LD *= -1;
}

void GrainHull::computeTriplelineLength() {
	double old_length = m_TripleLineLength;
	m_TripleLineLength = 0;
	for (vector<TripleLine*>::iterator iter = m_TripleLines.begin();
			iter != m_TripleLines.end(); ++iter) {
		vector<InterfacialElement*> vertices_temp = (*iter)->get_vertices();
		if (vertices_temp[0] == NULL || vertices_temp[1] == NULL) {
			continue;
		}
		m_TripleLineLength += (vertices_temp[0]->get_Position()
				- vertices_temp[1]->get_Position()).norm();
	}
	if (Settings::DecoupleGrains == 0) {
		vector<SPoint> ProjectedPoints;
		if (m_owner->getID() == 1) {

			for (const auto it2 : m_TripleLines) {
				for (const auto itBary : it2->m_barycenterTriangles) {
					SPoint newPoint((itBary)[0], (itBary)[2], (itBary)[1]);
					SPoint projection = newPoint.projectPointToPlane(
							SPoint(0, 0, 0.5), SPoint(1, 0, 0),
							SPoint(0, 1, 0));
					ProjectedPoints.push_back(projection);
				}
			}
			GrahamScan scanner(ProjectedPoints);
			vector<SPoint> ConvexHull;
			scanner.generateCovnexHull(ConvexHull);
			m_TripleLineLength = 0;

			for (vector<SPoint>::iterator it = ConvexHull.begin();
					it != (--ConvexHull.end()); it++) {
				//			for (int it =0 ; it < ConvexHull.size()-1; it++) {
				vector<SPoint>::iterator it2 = it + 1;
				m_TripleLineLength += (*it).DistanceTo(*(it2));
				//			m_TripleLineLength += ConvexHull[it].DistanceTo(ConvexHull[it+1]);
			}
			/*
			 unsigned long long timestep =
			 (unsigned long long) m_owner->get_grainHandler()->get_loop();
			 if (((timestep - Settings::StartTime)
			 % int(Settings::AnalysisTimestep * Settings::PlotInterval)) == 0
			 || timestep == Settings::NumberOfTimesteps) {

			 string filename = string("ConvexHull_") + to_string(timestep)
			 + string(".gnu");
			 FILE* output = fopen(filename.c_str(), "wt");
			 for (const auto it : ConvexHull) {
			 std::fprintf(output, "%lf \t %lf  \n", it.x, it.y);
			 }
			 fclose(output);

			 string filename1 = string("ProjectedPoints_")
			 + to_string((unsigned long long) timestep) + string(".gnu");
			 FILE* output1 = fopen(filename1.c_str(), "wt");

			 for (const auto it : ProjectedPoints) {
			 fprintf(output, "%lf \t %lf \n", it.x, it.y);
			 }
			 fclose(output1);
			 }
			 */
		}
	}
}

vector<Face>* GrainHull::get_Faces() {
	vector<Face>* myfaces = new vector<Face>;
	for (auto it : m_Grainboundary) {
		if (m_owner->getID() > it->m_neighborIDs[0])
			myfaces->push_back(
					Face(it->get_area(), m_owner->getID(),
							it->m_neighborIDs[0]));
	}
	return myfaces;
}
