/*
 * InterfacialElement.cpp
 *
 *  Created on: Nov 24, 2015
 *      Author: cm654063
 */

#include "InterfacialElement.h"
#include "box.h"
#include "grainhdl.h"
#include "marchingCubes.h"
#include "Settings.h"
#include "myQuaternion.h"
#include "grainHull.h"

InterfacialElement::InterfacialElement(int key, GrainHull* owner) :
	m_Key_NeighborList(key), m_owner(owner) {
}

InterfacialElement::~InterfacialElement() {
}
double InterfacialElement::computeMobilityMisori(double misori) {
	if (Settings::UseMobilityModel == 0)
		return 1.0;
	if (Settings::UseMobilityModel == 1)
		return 1 - (1.0 * exp(-5. * (pow(misori / (15 * PI / 180), 4.))));
	else
		return 1.0;
}

double InterfacialElement::computeReadShockleyEnergy(double misori) {
	if (Settings::UniqueGBEnergies == 1)
		return 1.0;
	double gamma_hagb = 1.0;
	double theta_ref = 15 * PI / 180;
	double gamma;
	if (misori > theta_ref)
		gamma = gamma_hagb;
	else
		gamma = gamma_hagb * (misori / theta_ref) * (1.0 - log(
				misori / theta_ref));
	return gamma;

}

GrainBoundary::GrainBoundary(int key, GrainHull* owner) :
	InterfacialElement(key, owner) {
	int currentID = owner->m_owner->getID();
	if (currentID
			!= m_owner->m_triangleNeighborLists[key].neighbors[0]) {
		m_neighborID = m_owner->m_triangleNeighborLists[key].neighbors[0];
	} else
		m_neighborID = m_owner->m_triangleNeighborLists[key].neighbors[1];
	computeMobility();
	computeEnergy();
}
GrainBoundary::~GrainBoundary() {
}

void GrainBoundary::computeEnergy() {
	LSbox* neighbor = m_owner->m_owner->get_grainHandler()->grains[m_neighborID];
	double misori = m_owner->m_owner->computeMisorientation(neighbor);
	m_energy = computeReadShockleyEnergy(misori);
}

void GrainBoundary::computeMobility() {
	LSbox* neighbor = m_owner->m_owner->get_grainHandler()->grains[m_neighborID];
	double misori = m_owner->m_owner->computeMisorientation(neighbor);
	m_mobility = computeMobilityMisori(misori);
}

TripleLine::TripleLine(int key, GrainHull *owner) :
	InterfacialElement(key, owner) {
	int currentID = owner->m_owner->getID();
	int j = 0;
	for (int i=0; i < 3; i++) {
		if (currentID
				!= m_owner->m_triangleNeighborLists[key].neighbors[i]) {
			m_neighborID[j] = m_owner->m_triangleNeighborLists[key].neighbors[i];
			j++;
		}
	}
	computeMobility();
	computeEnergy();
}
TripleLine::TripleLine(int neighbor1, int neighbor2, GrainHull* owner) :
	InterfacialElement(-1, owner) {
	m_neighborID[0] = neighbor1;
	m_neighborID[1] = neighbor2;
	computeMobility();
	computeEnergy();
}

TripleLine::~TripleLine() {
}

void TripleLine::computeEnergy() {
	double sigma;
	double gamma[3] = { 0.0, 0.0, 0.0 };
	double theta_mis;
	LSbox* me = m_owner->m_owner;
	grainhdl* handler = m_owner->m_owner->get_grainHandler();
	LSbox* neighborGrains[2] = { handler->getGrainByID(m_neighborID[0]),
			handler->getGrainByID(m_neighborID[1]) };
	theta_mis = me->computeMisorientation(neighborGrains[0]);
	gamma[0] = computeReadShockleyEnergy(theta_mis);

	theta_mis = neighborGrains[0]->computeMisorientation(neighborGrains[1]);
	gamma[1] = computeReadShockleyEnergy(theta_mis);

	theta_mis = me->computeMisorientation(neighborGrains[1]);
	gamma[2] = computeReadShockleyEnergy(theta_mis);

	sigma = gamma[0] - gamma[1] + gamma[2];

	if (sigma < 0.01) {
		//cout << "negative sigma " << endl;
		sigma = 0.01;
	}
	m_energy = sigma;
}

void TripleLine::computeMobility() {
	double averageMobility = 0;
	double theta_mis;
	LSbox* me = m_owner->m_owner;
	grainhdl* handler = m_owner->m_owner->get_grainHandler();
	LSbox* neighborGrains[2] = { handler->getGrainByID(m_neighborID[0]),
			handler->getGrainByID(m_neighborID[1]) };
	theta_mis = me->computeMisorientation(neighborGrains[0]);
	averageMobility += computeMobilityMisori(theta_mis);
	theta_mis = neighborGrains[0]->computeMisorientation(neighborGrains[1]);
	averageMobility += computeMobilityMisori(theta_mis);
	theta_mis = me->computeMisorientation(neighborGrains[1]);
	averageMobility += computeMobilityMisori(theta_mis);
	averageMobility /= 3;
	//TODO:
	double ds = 3 * sqrt(3 * handler->get_ds() * handler->get_ds());
	// ds is the extension of the Tripleline - maximum 3 times the diagonal of a grid cell
	m_mobility = 1 / ((1 / (ds * Settings::TripleLineDrag)) + 1

	/ averageMobility);
}

void GrainBoundary::findAdjacentTripleLines(vector<TripleLine*> Junctions){
	int i= 0;
	for(const auto it : Junctions)
	{
		if (it->get_FirstNeighbor() == m_neighborID || it->get_SecondNeighbor() == m_neighborID){
				m_edges[i]=&(*it); i++;
				if(i==1) return ;
			}
	}
}
void TripleLine::findAdjacentQuadrupleJunctions(vector<QuadrupleJunction*> Junctions){
	int i= 0;
	for(const auto it : Junctions)
	{
		if (it->get_FirstNeighbor() == m_neighborID[0] || it->get_SecondNeighbor() == m_neighborID[0] || it->get_ThirdNeighbor() == m_neighborID[0])
			if (it->get_FirstNeighbor() == m_neighborID[1] || it->get_SecondNeighbor() == m_neighborID[1] || it->get_ThirdNeighbor() == m_neighborID[1]){
				m_vertices[i]=&(*it); i++;
				if(i==2) return ;
			}
	}
}

QuadrupleJunction::QuadrupleJunction(int key, GrainHull* owner) :
	InterfacialElement(key, owner) {
	int currentID = owner->m_owner->getID();
	int j = 0;
	for (int i=0; i < 4; i++) {
		if (currentID
				!= m_owner->m_triangleNeighborLists[key].neighbors[i]) {
			m_neighborID[j] = m_owner->m_triangleNeighborLists[key].neighbors[i];
			j++;
		}
	}
	computeMobility();
	computeEnergy();
}

QuadrupleJunction::~QuadrupleJunction() {
}

void QuadrupleJunction::computeEnergy() {
	//compute average weight for all possible Triplets
	TripleLine T1(m_neighborID[0], m_neighborID[1], m_owner);
	TripleLine T2(m_neighborID[1], m_neighborID[2], m_owner);
	TripleLine T3(m_neighborID[0], m_neighborID[2], m_owner);
	m_energy = (T1.get_energy() + T2.get_energy() + T3.get_energy()) / 3;
}

void QuadrupleJunction::computeMobility() {
	//TODO::
	m_mobility = 1;
}

