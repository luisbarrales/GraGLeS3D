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
			// high order junction is found
			//TODO:
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
			//GB->addTriangle(m_actualHull[i]);
			break;
		}
		case 3: {
			TripleLine* TL = findTripleLine(key);
			Vector3d current = m_actualHull[i].computeBarycenter();
			TL->addBaryCenter(current);
			//TL->addTriangle(m_actualHull[i]);
			break;
		}
		case 4: {
			QuadrupleJunction* QJ = findQuadrupleJunction(key);
			Vector3d current = m_actualHull[i].computeBarycenter();
			QJ->addBaryCenter(current);
			//QJ->addTriangle(m_actualHull[i]);
			break;
		}
		default: {
			//cout << "high order junction found" << endl;
			break;
		}
		}
	}
	//m_actualHull.clear();
}
/*
bool IsNeighbor(Triangle &T1, Triangle &T2);

void calculateMeanWidthComponent(Triangle &T1, Triangle &T2, Vector3d normal1, Vector3d normal2, double &meanWidth);
*/
void GrainHull::computeInterfacialElementMesh() {
	//TODO:

	for (const auto it : m_TripleLines) {
		it->findAdjacentQuadrupleJunctions(m_QuadruplePoints);
	}
	for (const auto it : m_Grainboundary) {
		it->findAdjacentTripleLines(m_TripleLines);
	}
	//for the purpose of analyzing the geometric objects of the surface of the grain have to be explicitely computed
	// find the center of mass of a QuadruplePoint
	// attach the quadruplePoints to TripleLines
	// interpolate through the set of points describing the TL
	// attach these objects to the GrainBoundary
	// optimize the search routine to find nearest Triangle
	// extend the classes interfacial elements etc. to capture the analytic descriptions

	/*
	 * LD hier
	 * NeighborList local
	 * computeNormal
	 * nur LD speichern als Eigenschaft der GrainHull
	 */

	/*
	 * Save for every triangle the IDs of the neighbor-triangles
	 */
/*
	vector<int>* NeighborList = new vector<int>[m_actualHull.size()];

	for(int i = 0; i < m_actualHull.size(); i++){
		for(int j=0; j < i ; j++){
			if(IsNeighbor(m_actualHull[i],m_actualHull[j]))
				NeighborList[i].push_back(j);
		}
	}
*/
	/*
	 * Calculate the normal vector of every triangle
	 */
/*
	m_normalVectors.clear();
	Vector3d normal_temp;

	for(int i=0; i< m_actualHull.size(); i++){
		normal_temp = (m_actualHull[i].points[1]-m_actualHull[i].points[0]).cross((m_actualHull[i].points[2]-m_actualHull[i].points[0]));
		normal_temp /= normal_temp.norm();
		m_normalVectors.push_back(normal_temp);
	}
*/
	/*
	 * Calculate the mean width of the grain
	 */
/*
	for(int i=0; i<m_actualHull.size(); i++){
		for(int j=0; j<NeighborList[i].size(); j++){
			calculateMeanWidthComponent(m_actualHull[i],m_actualHull[(NeighborList[i])[j]],m_normalVectors[i],m_normalVectors[(NeighborList[i])[j]],m_LD);
		}
	}
	m_LD /= 2*M_PI;
*/
	/*
	 * Calculate the length of the triple lines
	 */ 
/*
	m_TripleLineLength=0;
 	vector<QuadrupleJunction*> vertices_temp;
	TripleLine* TripleLine_temp;
 	for(vector<TripleLine*>::iterator iter = m_TripleLines.begin(); iter != m_TripleLines.end(); ++iter){
		TripleLine_temp = *iter;
		vertices_temp = TripleLine_temp->get_vertices();
		m_TripleLineLength += (vertices_temp[0]->get_Position()-vertices_temp[0]->get_Position()).norm();
	}*/
}
/*
bool IsNeighbor(Triangle &T1, Triangle &T2){
	int NNeighbors=0;
	int t1[2];
	int t2[2];
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			if(T1.points[i] == T2.points[j]) {
				t1[NNeighbors] = i;
				t2[NNeighbors] = j;
				NNeighbors++;
			}
		}
	}
	if(NNeighbors == 2){
*/
		/*
 		 * Check if the triangles have the correct orientation
 		 */
/*
		if((t1[1]-t1[0]+3)%3-(t2[1]-t2[0]+3)%3==0){
			Vector3d temp;
			temp = T2.points[t2[1]];
			T2.points[t2[1]] = T2.points[t2[0]];
			T2.points[t2[0]] = temp;
		}
		return true;
	}
	else
		return false;
}

void calculateMeanWidthComponent(Triangle &T1, Triangle &T2, Vector3d normal1, Vector3d normal2, double &meanWidth){
	double scalarProduct;
	double beta, epsilon;

	scalarProduct = normal1.transpose() * normal2;
*/
	/*
	 * It is possible that through an numerical error the scalar product of two parallel vectors
	 * becomes greater than one. In this case the value of the scalar product is set to one.
	 */
/*
	if(scalarProduct>1){
		scalarProduct = 1.;
	}

	beta = acos(scalarProduct);
*/
	/*
	 * Find the points that are shared by both triangles
	 */
/*
	int sharedPoints[2];
	int temp=0;
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			if(T1.points[i]==T2.points[j]){
				sharedPoints[temp]=i;
				temp++;
			}
		}
	}
*/
	/*
	 * Orient the shared line as it is oriented in triangle T1
	 */
/*	Vector3d sharedLine;
	if((sharedPoints[1]-sharedPoints[0]+3)%3==1){
		sharedLine=T1.points[sharedPoints[1]]-T1.points[sharedPoints[0]];
	}
	else
	{
		sharedLine=T1.points[sharedPoints[0]]-T1.points[sharedPoints[1]];
	}

	epsilon = sharedLine.norm();
*/
	/*
	 * If the cross product of the two normal vectors is parallel to the shared line the angle is positive
	 * if it is antiparallel the angle is negative this is accounted for in the next if-statement
	 */
/*
	if((normal1.cross(normal2).transpose() * sharedLine) < 0)
	{
		beta*=-1;
	}
	if(beta != beta) {
		cout << "beta" << endl;
		cout << normal1.transpose() * normal2 << endl;
	}
	if(epsilon != epsilon) cout << "epsilon" << endl;
	meanWidth += beta*epsilon;
}
*/

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

double GrainHull::projectPointToGrainBoundary(Vector3d& point, int id) {
	double minimalDistance = 10000000.0;
	double weight = 1.0;

	//search in HighOrderJunctions
	//TODO:

	//search in QuadrupleJunctions
	for (int j = 0; j < m_QuadruplePoints.size(); j++) {
		if (m_QuadruplePoints[j]->m_neighborID[0] == id
				|| m_QuadruplePoints[j]->m_neighborID[1] == id
				|| m_QuadruplePoints[j]->m_neighborID[2] == id) {
			for (unsigned int i = 0;
					i < m_QuadruplePoints[j]->m_Triangles.size(); i++) {
				double distance = pointToTriangleDistance(point,
						m_QuadruplePoints[j]->m_Triangles[i]);
				if (distance < minimalDistance) {
					minimalDistance = distance;
					weight = m_QuadruplePoints[j]->get_Correction_Weight();
				}
			}
		}
	}
	//search in TripleJunctions
	for (int j = 0; j < m_TripleLines.size(); j++) {
		if (m_TripleLines[j]->m_neighborID[0] == id
				|| m_TripleLines[j]->m_neighborID[1] == id) {
			for (unsigned int i = 0; i < m_TripleLines[j]->m_Triangles.size();
					i++) {
				double distance = pointToTriangleDistance(point,
						m_TripleLines[j]->m_Triangles[i]);
				if (distance < minimalDistance) {
					minimalDistance = distance;
					weight = m_TripleLines[j]->get_Correction_Weight();
				}
			}
		}
	}
	//search in GrainBoundaries:
	for (int j = 0; j < m_Grainboundary.size(); j++) {
		if (m_Grainboundary[j]->m_neighborID == id) {
			for (unsigned int i = 0; i < m_Grainboundary[j]->m_Triangles.size();
					i++) {
				double distance = pointToTriangleDistance(point,
						m_Grainboundary[j]->m_Triangles[i]);
				if (distance < minimalDistance) {
					minimalDistance = distance;
					weight = m_Grainboundary[j]->get_Correction_Weight();
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
		throw runtime_error("Unable to save box hull!");
	}

	fprintf(output, "\n\nGRAINBOUNDARY %lu\n", m_Grainboundary.size());
	for (const auto it : m_Grainboundary) {
		for (const auto it2 : it->m_barycenterTriangles) {
			fprintf(output, "%lf \t %lf \t %lf \t %d \n ", it2[0], it2[1],
					it2[2], ID);
			fprintf(output, "\n\nEDGES %lu\n", it->m_edges.size());
			for (const auto it3 : it->m_edges) {
				for (const auto it4 : it3->m_barycenterTriangles) {
					fprintf(output, "%lf \t %lf \t %lf \t %d \n ", it4[0],
							it4[1], it4[2], ID);

					fprintf(output, "\n\nVERTICES %lu\n",
							it3->m_vertices.size());
					for (const auto it5 : it3->m_vertices) {
						Vector3d point = it5->get_Position();
						fprintf(output, "%lf \t %lf \t %lf \t %d \n ", point[0],
								point[1], point[2], ID);

						ID++;
					}
				}
				ID++;
			}
		}
		ID++;
	}
}
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
