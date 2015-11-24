/*
 * InterfacialElement.cpp
 *
 *  Created on: Nov 24, 2015
 *      Author: cm654063
 */

#include "InterfacialElement.h"

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
		gamma = gamma_hagb * (misori / theta_ref) * (1.0
				- log(misori / theta_ref));
	return gamma;
}

GrainBoundary::GrainBoundary(int key, NeighborList newboundary, LSbox* myBox) :
	Key_NeighborList(key) {
	int neighborID = newboundary.neighbors[0];
	double misori = myBox->computeMisorientation(
			myBox->m_grainHandler->grains[neighborID]);
	mobility = computeMobility(misori);
	energy = computeReadShockleyEnergy(misori);
}
GrainBoundary::~GrainBoundary(){
}

TripleLine::TripleLine(int key, NeighborList newTripleLine, LSbox* myBox) :
	Key_NeighborList(key) {
}
TripleLine::~TripleLine(){
}

QuadrupleJunction::QuadrupleJunction(int key,
		NeighborList newQuadrupleJunction, LSbox* myBox) :
	Key_NeighborList(key) {
}

QuadrupleJunction::~QuadrupleJunction() {
}
