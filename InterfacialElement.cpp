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

InterfacialElement::InterfacialElement(int key, GrainHull owner) :
	m_Key_NeighborList(key), m_owner(owner) {
}

InterfacialElement::~InterfacialElement() {
}
double InterfacialElement::computeMobility(double misori) {
	//	if (Settings::UseMobilityModel == 1) {
	//check for twin boundary
	// 8.66025 = Theta (15)* 1/ sqrt(SIGMA) here SiGMA is the number od coincidence points (3 for the Twinboundary)
	// thus 8.66 is the permissible deviation from the perfect Twinboundary SIGMA 3
	//		if (Settings::IdentifyTwins) {
	//							if (MisoriToTwinBoundary(candidate) < 8.66025 * PI / 180.0)
	//							return 0.01;
	//							else
	//			return 1 - (1.0 * exp(-5. * (pow(thetaMis / (15 * PI / 180), 4.))));
	//		} else
	return 1 - (1.0 * exp(-5. * (pow(misori / (15 * PI / 180), 4.))));
	//	} else
	//		return 1;
}

double InterfacialElement::computeReadShockleyEnergy(double misori) {
	double gamma_hagb = Settings::HAGB;
	double theta_ref = 15 * PI / 180;
	double gamma;
	if (misori > theta_ref)
		gamma = gamma_hagb;
	else
		gamma = gamma_hagb * (misori / theta_ref) * (1.0 - log(
				misori / theta_ref));
	return gamma;
}

GrainBoundary::GrainBoundary(int key, GrainHull owner) :
	InterfacialElement(key, owner) {
	neighborID = m_owner->m_triangleNeighborLists[key].neighbors[0];
	computeMobility();
	computeEnergy();
}
GrainBoundary::~GrainBoundary() {
}

void GrainBoundary::computeEnergy() {
	LSbox* neighbor = m_owner->m_owner->get_grainHandler()->grains[neighborID];
	double misori = m_owner->m_owner->computeMisorientation(neighbor);
	energy = computeReadShockleyEnergy(misori);
}

void GrainBoundary::computeMobility() {
	mobility = computeMobility(misori);
}

TripleLine::TripleLine(int key, GrainHull owner) :
	InterfacialElement(key, owner) {
	neighborID[0] = m_owner->m_triangleNeighborLists[key].neighbors[0];
	neighborID[1] = m_owner->m_triangleNeighborLists[key].neighbors[1];
	computeMobility();
	computeEnergy();
}
TripleLine::~TripleLine() {
}

void TripleLine::computeEnergy() {
	double averageMobility = 0;
	double sigma;
	double gamma[3] = { 0.0, 0.0, 0.0 };
	double gamma_hagb = Settings::HAGB;
	double theta_ref = 15.0 * PI / 180;
	double theta_mis;
	LSbox* me = m_owner->m_owner;
	grainhdl* handler = m_owner->m_owner->m_grainHandler;
	LSbox* neighborGrains[2] = { handler->getGrainByID(neighborID[0]),
			handler->getGrainByID(neighborID[1]) };

	if (Settings::ResearchMode != 1) {

		theta_mis = me->computeMisorientation(neighborGrains[0]);
		averageMobility += computeMobility(theta_mis);
		gamma[0] = computeReadShockleyEnergy(theta_mis);

		theta_mis = neighborGrains[0]->computeMisorientation(neighborGrains[1]);
		averageMobility += computeMobility(theta_mis);
		gamma[1] = computeReadShockleyEnergy(theta_mis);

		theta_mis = me->computeMisorientation(neighborGrains[1]);
		averageMobility += computeMobility(theta_mis);
		gamma[2] = computeReadShockleyEnergy(theta_mis);

	} else if (Settings::ResearchMode == 1) {
		gamma[0] = 1.0;
		gamma[1] = 1.0;
		gamma[2] = 1.0;
	}

	// find the asociated weight
	sigma = gamma[0] - gamma[1] + gamma[2];

	if (Settings::UseMobilityModel > 0 && Settings::TripleLineDrag > 0) {
		averageMobility /= 3;
		double ds = handler->get_ds();
		double drag = 1 / ((1 / (ds * Settings::TripleLineDrag)) + 1
				/ averageMobility);
		sigma *= drag;
	}
	if (sigma < 0.01) {
		//cout << "negative sigma " << endl;
		sigma = 0.01;
	}
	m_energy = sigma;
}

void TripleLine::computeMobility() {
	m_mobility = 1;
	//the triple Line will move with an effective mobility computed above
}

QuadrupleJunction::QuadrupleJunction(int key,
		NeighborList newQuadrupleJunction, LSbox* myBox) :
	InterfacialElement(key, newQuadrupleJunction) {
	neighborID[0] = m_owner->m_triangleNeighborLists[key].neighbors[0];
	neighborID[1] = m_owner->m_triangleNeighborLists[key].neighbors[1];
	neighborID[2] = m_owner->m_triangleNeighborLists[key].neighbors[2];
	computeMobility();
	computeEnergy();
}

QuadrupleJunction::~QuadrupleJunction() {
}

void QuadrupleJunction::computeEnergy() {
	//TODO::
	m_energy =1;
}

void QuadrupleJunction::computeMObility() {
	//TODO::
	m_mobility =1;
}

