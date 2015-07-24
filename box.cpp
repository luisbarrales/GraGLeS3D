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
#include "box.h"
#include "Settings.h"
#include "dimensionalBufferReal.h"
#include "pooledDimensionalBufferReal.h"
#include "contourSector.h"
#include "minimalisticBoundary.h"
#include "utilities.h"
#include "grainhdl.h"
#include "mymath.h"
#include "marchingCubes.h"
#include <algorithm>
#include <cstdio>
#include <stdexcept>

#define PERIODIC(x, f) (((x)+f)%f)

LSbox::LSbox(int id, double phi1, double PHI, double phi2, grainhdl* owner) :
	m_ID(id), m_exists(true),m_grainHandler(owner), m_explicitBoundary(this),
	m_isMotionRegular(true), m_intersectsBoundaryGrain(false),
	m_volume(0), m_energy(0), m_perimeter(0)
{
	m_orientationQuat = new double[4];
	double euler[3] = { phi1, PHI, phi2 };
	(*(m_grainHandler->mymath)).euler2quaternion(euler, m_orientationQuat);

	m_inputDistance = new DimensionalBufferReal(0, 0, 0, 0, 0, 0);
	m_outputDistance = new DimensionalBufferReal(0, 0, 0, 0, 0, 0);

	m_volumeEvolution = rnd() * 100; //Zufallszahl zwischen 0 und 100
}

LSbox::LSbox(int id, vector<Vector3d>& hull, grainhdl* owner) :
	m_ID(id), m_exists(true),m_grainHandler(owner), m_explicitBoundary(this),
	m_isMotionRegular(true), m_intersectsBoundaryGrain(false),
	m_volume(0), m_energy(0), m_perimeter(0)
{
	int grid_blowup = owner->get_grid_blowup();
	m_volumeEvolution = rnd() * 100; //Zufallszahl zwischen 0 und 100
	double h = owner->get_h();
	// determine size of grain
	m_orientationQuat = new double[4];
	#pragma omp critical
		{
			if (Settings::UseTexture)
			{
				double newOri[3];
				(*(m_grainHandler->mymath)).newOrientationFromReference(m_grainHandler->bunge,
						m_grainHandler->deviation, newOri);
				(*(m_grainHandler->mymath)).euler2quaternion(newOri, m_orientationQuat);
			}
			else
				(*(m_grainHandler->mymath)).randomOriShoemakeQuat(m_orientationQuat);
		}
		int xmax = 0;
		int xmin = m_grainHandler->get_ngridpoints();
		int ymax = 0;
		int ymin = xmin;
		int zmin = ymin;
		int zmax = 0;

		double x, y, z;
		for (unsigned int k = 0; k < hull.size(); k++) {
			x = hull.at(k)[0];
			y = hull.at(k)[1];
			z = hull.at(k)[2];

			if (y / h < ymin)
				ymin = grid_blowup + y / h;
			if (y / h > ymax)
				ymax = grid_blowup +  y / h + 1;
			if (x / h < xmin)
				xmin = grid_blowup +  x / h;
			if (x / h > xmax)
				xmax = grid_blowup +  x / h + 1;
			if (z / h < zmin)
				zmin = grid_blowup +  z / h;
			if (z / h > zmax)
				zmax = grid_blowup +  z / h + 1;
		}

		xmax += grid_blowup;
		xmin -= grid_blowup;
		ymax += grid_blowup;
		ymin -= grid_blowup;
		zmax += grid_blowup;
		zmin -= grid_blowup;

		m_inputDistance = new DimensionalBufferReal(xmin, ymin, zmin, xmax, ymax, zmax);
		m_outputDistance = new DimensionalBufferReal(xmin, ymin, zmin, xmax, ymax, zmax);

		m_inputDistance->resizeToCube(m_grainHandler->get_ngridpoints());
		m_outputDistance->resizeToCube(m_grainHandler->get_ngridpoints());

		reizeIDLocalToDistanceBuffer();
}

LSbox::LSbox(int id, DimensionalBuffer<int>& IDField, grainhdl* owner) :
	m_ID(id), m_exists(true),m_grainHandler(owner), m_explicitBoundary(this),
	m_isMotionRegular(true), m_intersectsBoundaryGrain(false),
	m_volume(0), m_energy(0), m_perimeter(0)
{
	int grid_blowup = owner->get_grid_blowup();
	m_volumeEvolution = rnd() * 100; //Zufallszahl zwischen 0 und 100
	double h = owner->get_h();
	// determine size of grain
	m_orientationQuat = new double[4];
	#pragma omp critical
		{
			if (Settings::UseTexture)
			{
				double newOri[3];
				(*(m_grainHandler->mymath)).newOrientationFromReference(m_grainHandler->bunge,
						m_grainHandler->deviation, newOri);
				(*(m_grainHandler->mymath)).euler2quaternion(newOri, m_orientationQuat);
			}
			else
				(*(m_grainHandler->mymath)).randomOriShoemakeQuat(m_orientationQuat);
		}
		int xmax = 0;
		int xmin = m_grainHandler->get_ngridpoints();
		int ymax = 0;
		int ymin = xmin;
		int zmin = ymin;
		int zmax = 0;

		for(int i = IDField.getMinY(); i< IDField.getMaxY(); i++)
			for(int j = IDField.getMinX(); j< IDField.getMaxX(); j++)
				for(int k = IDField.getMinZ(); k< IDField.getMaxZ(); k++)
				{
					if(m_ID == IDField.getValueAt(i,j,k))
					{
						xmax = max(j, xmax);
						xmin = min(j, xmin);
						ymax = max(i, ymax);
						ymin = min(i, ymin);
						zmax = max(k, zmax);
						zmin = min(k, zmin);
					}
				}

		xmax += grid_blowup;
		xmin -= grid_blowup;
		ymax += grid_blowup;
		ymin -= grid_blowup;
		zmax += grid_blowup;
		zmin -= grid_blowup;

		m_inputDistance = new DimensionalBufferReal(xmin, ymin, zmin, xmax, ymax, zmax);
		m_outputDistance = new DimensionalBufferReal(xmin, ymin, zmin, xmax, ymax, zmax);

		m_inputDistance->resizeToCube(m_grainHandler->get_ngridpoints());
		m_outputDistance->resizeToCube(m_grainHandler->get_ngridpoints());

		reizeIDLocalToDistanceBuffer();
}
LSbox::~LSbox() {
	if (m_orientationQuat != NULL)
		delete[] m_orientationQuat;
	delete m_inputDistance;
	delete m_outputDistance;
}

double LSbox::get_h() {
	return m_grainHandler->get_h();
}

//TODO: Distance function must be properly implemented using point to mesh distance
void LSbox::calculateDistanceFunction(DimensionalBuffer<int>& IDField)
{
	int min = m_grainHandler->get_grid_blowup();
	int max = m_grainHandler->get_ngridpoints() - min - 1;
	for (int k = m_inputDistance->getMinZ(); k < m_inputDistance->getMaxZ(); k++) {
			for (int i = m_inputDistance->getMinY(); i < m_inputDistance->getMaxY(); i++) {
				for (int j = m_inputDistance->getMinX(); j < m_inputDistance->getMaxX(); j++) {
					if(i < min || i > max || j < min || j > max || k < min || k > max)
						m_inputDistance->setValueAt(i, j, k, -m_grainHandler->get_h());
					else if (m_ID == m_grainHandler->IDField.getValueAt(i, j, k))
						m_inputDistance->setValueAt(i, j, k, m_grainHandler->get_h());
					else
						m_inputDistance->setValueAt(i, j, k, -m_grainHandler->get_h());
				}
			}
		}
	executeRedistancing();
}

// Convolution und Helperfunctions 
/**************************************/
/**************************************/

void LSbox::executeConvolution(ExpandingVector<char>& mem_pool) {

	if (grainExists() != true)
		return;
	//  set references for the convolution step
	fftwp_complex *fftTemp = (fftwp_complex*) &mem_pool[0];
	convolutionGeneratorFFTW(fftTemp, m_forwardPlan, m_backwardsPlan);

	reizeIDLocalToDistanceBuffer();
	m_IDLocal.clear();
}

void LSbox::cleanupConvolution() {
	fftw_destroy_planp(m_forwardPlan);
	fftw_destroy_planp(m_backwardsPlan);
}

double LSbox::getGBEnergyTimesGBMobility(int i, int j) {
	return -1;
}

double LSbox::getGBEnergyTimesGBMobility(LSbox* neighbour) {
	return -1;
}

double LSbox::getGBEnergy(LSbox* neighbour) {
	return -1;
}

void LSbox::reizeIDLocalToDistanceBuffer() {
	int xmaxId = m_outputDistance->getMaxX();
	int xminId = m_outputDistance->getMinX();
	int ymaxId = m_outputDistance->getMaxY();
	int yminId = m_outputDistance->getMinY();
	int zminId = m_outputDistance->getMinZ();
	int zmaxId = m_outputDistance->getMaxZ();
	m_IDLocal.resize(xminId, yminId, zminId, xmaxId, ymaxId, zmaxId);
}
void LSbox::makeFFTPlans(double *in, double* out, fftw_complex *fftTemp,
		fftw_plan *fftplan1, fftw_plan *fftplan2) { /* creates plans for FFT and IFFT */
	int n = m_outputDistance->getMaxX() - m_outputDistance->getMinX();
	int m = m_outputDistance->getMaxY() - m_outputDistance->getMinY();
	int l = m_outputDistance->getMaxZ() - m_outputDistance->getMinZ();

	*fftplan1 = fftw_plan_dft_r2c_3d(n, m, l, in, fftTemp, FFTW_ESTIMATE);
	*fftplan2 = fftw_plan_dft_c2r_3d(n, m, l, fftTemp, out, FFTW_ESTIMATE);

	/*
	 The flags argument is usually either FFTW_MEASURE or FFTW_ESTIMATE. FFTW_MEASURE
	 instructs FFTW to run and measure the execution time of several FFTs in order to find the
	 best way to compute the transform of size n. This process takes some time (usually a few
	 seconds), depending on your machine and on the size of the transform. FFTW_ESTIMATE,
	 on the contrary, does not run any computation and just builds a reasonable plan that is
	 probably sub-optimal. In short, if your program performs many transforms of the same size
	 and initialization time is not important, use FFTW_MEASURE; otherwise use the estimate. */
}

void LSbox::makeFFTPlans(float *in, float* out, fftwf_complex *fftTemp,
		fftwf_plan *fftplan1, fftwf_plan *fftplan2) { /* creates plans for FFT and IFFT */
	int n = m_outputDistance->getMaxX() - m_outputDistance->getMinX();
	int m = m_outputDistance->getMaxY() - m_outputDistance->getMinY();
	int l = m_outputDistance->getMaxZ() - m_outputDistance->getMinZ();

	*fftplan1 = fftwf_plan_dft_r2c_3d(n, m, l, in, fftTemp, FFTW_ESTIMATE);
	*fftplan2 = fftwf_plan_dft_c2r_3d(n, m, l, fftTemp, out, FFTW_ESTIMATE);

}

void LSbox::createConvolutionPlans(ExpandingVector<char>& memory_dump) {
	fftwp_complex *fftTemp = (fftwp_complex*) &memory_dump[0];

	makeFFTPlans(m_inputDistance->getRawData(), m_outputDistance->getRawData(),
			fftTemp, &m_forwardPlan, &m_backwardsPlan);
}


void LSbox::initConvoMemory(ExpandingVector<char>& memory_dump) {
	int n = m_outputDistance->getMaxX() - m_outputDistance->getMinX();
	int m = m_outputDistance->getMaxY() - m_outputDistance->getMinY();
	int l = m_outputDistance->getMaxZ() - m_outputDistance->getMinZ();
	int desired_size = n * m * (floor(l / 2) + 1) * sizeof(fftwp_complex);
	memory_dump.expand(desired_size);
}

void LSbox::convolutionGeneratorFFTW(fftwp_complex *fftTemp, fftwp_plan fftplan1,
		fftwp_plan fftplan2) {

	int n = m_outputDistance->getMaxX() - m_outputDistance->getMinX();
	double dt = m_grainHandler->get_dt();
	int n2 = floor(n / 2) + 1;
	int nn = m_grainHandler->get_realDomainSize();
	double nsq = nn * nn *nn;
	double G;
	int j2;
	int i2;

	executeFFTW(fftplan1);
	//	Forward DFT

	switch (Settings::ConvolutionMode)
	{
		case E_GAUSSIAN:
		{
			double n_nsq = n*n*n;
			//			Convolution with Normaldistribution
			for (int k = 0; k < n; k++) {
				int k2 = min(k,n-k);
				for (int i = 0; i < n; i++) {
					i2 = min(i,n-i);
					for (int j = 0; j < n2; j++) {
						j2 = min(j,n-j);
						G = exp(-(i2 * i2 + j2 * j2 + k2 * k2) * 8.0 * dt * nsq / n_nsq * PI * PI * PI) / n_nsq;
						fftTemp[j + n2 * (i + n * k)][0] = fftTemp[j + n2 * (i + n * k)][0] * G;
						fftTemp[j + n2 * (i + n * k)][1] = fftTemp[j + n2 * (i + n * k)][1] * G;
					}
				}
			}
			break;
		}
	default:
		throw runtime_error("Unknown convolution mode!");
	}

	executeFFTW(fftplan2);
	//	Inverse DFT
}


void LSbox::executeFFTW(fftw_plan fftplan) {
	fftw_execute(fftplan);
}

void LSbox::executeFFTW(fftwf_plan fftplan) {
	fftwf_execute(fftplan);
}

/**************************************/
/**************************************/

/**************************************/
/**************************************/

void LSbox::switchInNOut() {
	DimensionalBufferReal* temp;

	temp = m_inputDistance;
	m_inputDistance = m_outputDistance;
	m_outputDistance = temp;
}

// Comparison + Helperfunctions
/**************************************/
/**************************************/

void LSbox::executeSetComparison()
{
	m_newXMin = m_outputDistance->getMaxX();
	m_newXMax = m_outputDistance->getMinX();
	m_newYMin = m_outputDistance->getMaxY();
	m_newYMax = m_outputDistance->getMinY();
	m_newZMin = m_outputDistance->getMaxZ();
	m_newZMax = m_outputDistance->getMinZ();

	for (int k = m_outputDistance->getMinZ(); k < m_outputDistance->getMaxZ(); k++)
	{
		for (int i = m_outputDistance->getMinY(); i < m_outputDistance->getMaxY(); i++)
		{
			for (int j = m_outputDistance->getMinX(); j < m_outputDistance->getMaxX(); j++)
			{

				if (abs(m_inputDistance->getValueAt(i, j, k)) < 0.7 * m_grainHandler->delta)
				{
					m_outputDistance->setValueAt(i,j,k,
							0.5 * (m_inputDistance->getValueAt(i, j, k) - m_outputDistance->getValueAt(i, j, k)));
				}
				else
					m_outputDistance->setValueAt(i, j, k, m_inputDistance->getValueAt(i, j, k));

				if (m_outputDistance->getValueAt(i, j, k) >= 0)
				{
					if (i < m_newYMin)
						m_newYMin = i;
					if (i > m_newYMax)
						m_newYMax = i;
					if (j < m_newXMin)
						m_newXMin = j;
					if (j > m_newXMax)
						m_newXMax = j;
					if(k < m_newZMin)
						m_newZMin = k;
					if(k > m_newZMax)
						m_newZMax = k;
				}
			}
		}
	}

	m_newXMin -= m_grainHandler->get_grid_blowup();
	m_newXMax += m_grainHandler->get_grid_blowup();
	m_newYMin -= m_grainHandler->get_grid_blowup();
	m_newYMax += m_grainHandler->get_grid_blowup();
	m_newZMin -= m_grainHandler->get_grid_blowup();
	m_newZMax += m_grainHandler->get_grid_blowup();
}

bool LSbox::checkIntersection(LSbox* box2) {
	if        (m_inputDistance->getMinX() > box2->m_inputDistance->getMaxX()
			|| m_inputDistance->getMaxX() < box2->m_inputDistance->getMinX()
			|| m_inputDistance->getMinY() > box2->m_inputDistance->getMaxY()
			|| m_inputDistance->getMaxY() < box2->m_inputDistance->getMinY()
			|| m_inputDistance->getMinZ() > box2->m_inputDistance->getMaxZ()
			|| m_inputDistance->getMaxZ() < box2->m_inputDistance->getMinZ())
		return false;
	return true;
}

void LSbox::executeComparison()
{
	if (grainExists() != true)
		return;
	m_outputDistance->clearValues(-1.0);
	m_secondOrderNeighbours = m_comparisonList;

	for (unsigned int neighs = 0; neighs < m_secondOrderNeighbours.size(); neighs++)
	{
		LSbox* neighbor = m_grainHandler->getGrainByID(m_secondOrderNeighbours[neighs]);
		int x_min_new, x_max_new, y_min_new, y_max_new, z_min_new, z_max_new;

		if (m_inputDistance->getMinX() < neighbor->m_inputDistance->getMinX())
			x_min_new = neighbor->m_inputDistance->getMinX();
		else
			x_min_new = m_inputDistance->getMinX();

		if (m_inputDistance->getMaxX() > neighbor->m_inputDistance->getMaxX())
			x_max_new = neighbor->m_inputDistance->getMaxX();
		else
			x_max_new = m_inputDistance->getMaxX();

		if (m_inputDistance->getMinY() < neighbor->m_inputDistance->getMinY())
			y_min_new = neighbor->m_inputDistance->getMinY();
		else
			y_min_new = m_inputDistance->getMinY();

		if (m_inputDistance->getMaxY() > neighbor->m_inputDistance->getMaxY())
			y_max_new = neighbor->m_inputDistance->getMaxY();
		else
			y_max_new = m_inputDistance->getMaxY();

		if (m_inputDistance->getMinZ() < neighbor->m_inputDistance->getMinZ())
			z_min_new = neighbor->m_inputDistance->getMinZ();
		else
			z_min_new = m_inputDistance->getMinZ();

		if (m_inputDistance->getMaxZ() > neighbor->m_inputDistance->getMaxZ())
			z_max_new = neighbor->m_inputDistance->getMaxZ();
		else
			z_max_new = m_inputDistance->getMaxZ();

		for (int k = z_min_new; k < z_max_new; k++)
		{
			for (int i = y_min_new; i < y_max_new; i++)
			{
				for (int j = x_min_new; j < x_max_new; j++)
				{
					double dist = neighbor->getDistanceFromInputBuff(i, j, k);
					if (dist > m_outputDistance->getValueAt(i, j, k))
					{
							m_outputDistance->setValueAt(i, j, k, dist);
							m_IDLocal.getValueAt(i, j, k).grainID = neighbor->getID();
					}
				}
			}
		}
	}

	if (BoundaryIntersection())
	{
		m_intersectsBoundaryGrain = true;
		boundaryCondition();
	}
	else
		m_intersectsBoundaryGrain = false;
}

bool LSbox::BoundaryIntersection() {
	int xMinBoundary = m_grainHandler->get_grid_blowup()
			+ m_grainHandler->getBoundaryGrainTube();
	int yMinBoundary = xMinBoundary;
	int zMinBoundary = xMinBoundary;

	int xMaxBoundary = m_grainHandler->get_ngridpoints() - m_grainHandler->get_grid_blowup()
			- m_grainHandler->getBoundaryGrainTube();
	int yMaxBoundary = xMaxBoundary;
	int zMaxBoundary = xMaxBoundary;

	if (m_outputDistance->getMinX() > xMinBoundary && m_outputDistance->getMaxX()
			< xMaxBoundary && m_outputDistance->getMinY() > yMinBoundary
			&& m_outputDistance->getMaxY() < yMaxBoundary
			&& m_outputDistance->getMinZ() > zMinBoundary
			&& m_outputDistance->getMaxZ() < zMaxBoundary)
		return false;
	else
	{
		return true;
	}
}

double LSbox::getDistanceFromInputBuff(int i, int j, int k)
{
	return m_inputDistance->getValueAt(i,j,k);
}

//Differs largely from latest version of the 2D code! Needs to be checked
void LSbox::boundaryCondition()
{
	int grid_blowup = m_grainHandler->get_grid_blowup();
	double h = m_grainHandler->get_h();
	int m = m_grainHandler->get_ngridpoints();
	double distZMin, distZMax, distXMin, distXMax, distX;
	double distYMin, distYMax, distY, distZ;
	double dist;
	for (int k = m_inputDistance->getMinZ(); k < m_inputDistance->getMaxZ(); k++)
	{
		for (int i = m_inputDistance->getMinY(); i < m_inputDistance->getMaxY(); i++)
		{
			for (int j = m_inputDistance->getMinX(); j < m_inputDistance->getMaxX(); j++)
			{
				distXMin = -(j - grid_blowup);
				distYMin = -(i - grid_blowup);
				distZMin = -(k - grid_blowup);
				distXMax = (j - (m - grid_blowup));
				distYMax = (i - (m - grid_blowup));
				distZMax = (k - (m - grid_blowup));

				if (abs(distXMin) < abs(distXMax))
					distX = distXMin;
				else
					distX = distXMax;

				if (abs(distYMin) < abs(distYMax))
					distY = distYMin;
				else
					distY = distYMax;

				if (abs(distZMin) < abs(distZMax))
					distZ = distZMin;
				else
					distZ = distZMax;

				// the point is outside in one of the 8 corners:
				if (distX > 0 && distY > 0 && distZ > 0)
					dist = sqrt((double) distX * distX + distY * distY + distZ * distZ);

				// the point is inside the domain - > value is maximum of the negative distances:
				else if (distX < 0 && distY < 0 && distZ < 0) {
					dist = max (distX, distY);
					dist = max(dist, distZ);
				}

				// the point is outside in x direction
				else if (distX >= 0)
				{
					if (distY < 0 && distZ < 0) // next to one x-plane
						dist = distX;
					else if (distY < 0 && distZ >= 0)
						dist = sqrt((double) distX * distX + distZ * distZ);
					else if (distY >= 0 && distZ < 0)
						dist = sqrt((double) distX * distX + distY * distY);
				}

				else if (distY >= 0)
				{
					if (distZ < 0 && distX < 0) // next to one y-plane
						dist = distY;
					else if (distX < 0 && distZ >= 0)
						dist = sqrt((double) distY * distY + distZ * distZ);
				}

				else if (distZ >= 0)
				{
					if (distY < 0 && distX < 0) // next to one z-plane
						dist = distZ;
				}

				else if (distX == 0 || distY == 0 || distZ == 0)
					dist = 0; // one or more are zero the other lower

				if (dist > -grid_blowup)
					if (dist * h > m_outputDistance->getValueAt(i, j, k))
					{
						m_outputDistance->setValueAt(i, j, k, dist * h);
						m_IDLocal.getValueAt(i,j,k).grainID = 0;
					}
			}
		}
	}
}

void LSbox::computeSecondOrderNeighbours() {
	vector<unsigned int> neighbourCandidates;
	if (grainExists() != true)
		return;

	bool just_in;
	m_comparisonList.clear();

	for (unsigned int i = 0; i< m_secondOrderNeighbours.size(); i++)
	{
		if (m_grainHandler->getGrainByID(m_secondOrderNeighbours[i])->grainExists() == true)
			m_comparisonList.push_back(m_secondOrderNeighbours[i]);

		LSbox* currentNeighbor = m_grainHandler->getGrainByID(m_secondOrderNeighbours[i]);

		for (unsigned int j = 0; j < currentNeighbor->m_secondOrderNeighbours.size(); j++)
		{
			LSbox* hisNeighbor = m_grainHandler->getGrainByID(currentNeighbor->m_secondOrderNeighbours[j]);
			if (hisNeighbor->grainExists() == true)
				if (checkIntersection(hisNeighbor))
				{
					neighbourCandidates.push_back(hisNeighbor->getID());
				}
		}
	}
	for (unsigned int i = 0; i < neighbourCandidates.size(); i++)
	{
		just_in = false;
		if (neighbourCandidates[i] == getID())
			continue;
		if (m_grainHandler->getGrainByID(neighbourCandidates[i]) == m_grainHandler->boundary)
			continue;

		for (unsigned int j = 0; j < m_comparisonList.size(); j++) {
			if (m_comparisonList[j] == neighbourCandidates[i])
			{
				just_in = true;
				break;
			}
		}
		if ((!just_in))
			m_comparisonList.push_back(neighbourCandidates[i]);
	}
	neighbourCandidates.clear();
}

/**************************************/
// end of Comparison
/**************************************/

// Find Contour operates on inputDistance
/**************************************/
/**************************************/

void LSbox::extractContour() {
	MarchingCubesAlgorithm marcher(*m_inputDistance, this);
	marcher.generateHull(m_grainHull);
	if(m_grainHull.size() == 0)
	{
		m_exists = false;
		return;
	}

	if(m_exists)
	{
		m_outputDistance->resize(m_newXMin, m_newYMin, m_newZMin, m_newXMax, m_newYMax, m_newZMax);
		m_outputDistance->resizeToCube(m_grainHandler->get_ngridpoints());
	}

	return;
}

void LSbox::updateFirstOrderNeigbors() {
	if (grainExists() != true)
		return;
	return;
}
double LSbox::computeVolume() {
	return -1;
}

void LSbox::computeVolumeAndEnergy() {
	if (grainExists() != true || m_isMotionRegular == false)
		return;
}

/**************************************/
//  Redistancing
/**************************************/
void LSbox::executeRedistancing() {
	if (grainExists() != true)
		return;

	double h = m_grainHandler->get_h();
	double candidate, i_slope, distToZero;

	m_outputDistance->clearValues(-1.0);

	// 	resize the outputDistance array. be careful because during this part of algorithm both arrays have not the same size!!
	int intersec_xmin, intersec_xmax, intersec_ymin, intersec_ymax, intersec_zmin, intersec_zmax;

	if (m_inputDistance->getMinX() < m_outputDistance->getMinX())
		intersec_xmin = m_outputDistance->getMinX();
	else
		intersec_xmin = m_inputDistance->getMinX();

	if (m_inputDistance->getMinY() < m_outputDistance->getMinY())
		intersec_ymin = m_outputDistance->getMinY();
	else
		intersec_ymin = m_inputDistance->getMinY();

	if (m_inputDistance->getMinZ() < m_outputDistance->getMinZ())
		intersec_zmin = m_outputDistance->getMinZ();
	else
		intersec_zmin = m_inputDistance->getMinZ();

	if (m_inputDistance->getMaxX() < m_outputDistance->getMaxX())
		intersec_xmax = m_inputDistance->getMaxX();
	else
		intersec_xmax = m_outputDistance->getMaxX();

	if (m_inputDistance->getMaxY() < m_outputDistance->getMaxY())
		intersec_ymax = m_inputDistance->getMaxY();
	else
		intersec_ymax = m_outputDistance->getMaxY();

	if (m_inputDistance->getMaxZ() < m_outputDistance->getMaxZ())
		intersec_zmax = m_inputDistance->getMaxZ();
	else
		intersec_zmax = m_outputDistance->getMaxZ();

	// first to updates layer by layer to take advantage of the order of point in memory - there are aligned layer by layer.
	for (int k = intersec_zmin; k < intersec_zmax - 1; k++)
	{
		for (int i = intersec_ymin; i < m_outputDistance->getMaxY(); i++)
		{
			for (int j = intersec_xmin; j < m_outputDistance->getMaxX() - 1; j++)
			{
				// x-direction forward
				if (j < intersec_xmax - 1 && i < intersec_ymax)
				{
					if (m_inputDistance->getValueAt(i, j, k)* m_inputDistance->getValueAt(i, j + 1, k) <= 0.0)
					{
						// interpolate
						i_slope = (m_inputDistance->getValueAt(i, j + 1, k) - m_inputDistance->getValueAt(i, j, k)) / h;
						distToZero = -m_inputDistance->getValueAt(i, j, k) / i_slope;
						if (abs(m_outputDistance->getValueAt(i, j, k)) > abs(distToZero))
							m_outputDistance->setValueAt(i, j, k,	-distToZero * sgn(i_slope));
					}
					candidate = m_outputDistance->getValueAt(i, j, k) + (sgn(m_inputDistance->getValueAt(i, j + 1, k)) * h);
					if (abs(candidate) < abs(	m_outputDistance->getValueAt(i, j + 1, k)))
						m_outputDistance->setValueAt(i, j + 1, k, candidate);
				}
				else
				{
					candidate = m_outputDistance->getValueAt(i, j, k) + (sgn(m_outputDistance->getValueAt(i, j + 1, k))	* h);
					if (abs(candidate) < abs(m_outputDistance->getValueAt(i, j + 1, k)))
						m_outputDistance->setValueAt(i, j + 1, k, candidate);
				}
			}
		}

		for (int i = intersec_ymin; i < m_outputDistance->getMaxY(); i++)
		{
			for (int j = intersec_xmax - 1; j > m_outputDistance->getMinX(); j--)
			{
				// x-direction outputDistanceward
				//check for sign change
				if (j > intersec_xmin && i < intersec_ymax)
				{
					// calculate new distance candidate and assign if appropriate
					candidate = m_outputDistance->getValueAt(i, j, k) + (sgn(m_inputDistance->getValueAt(i, j - 1, k)) * h);
					if (abs(candidate) < abs(m_outputDistance->getValueAt(i, j - 1, k)))
						m_outputDistance->setValueAt(i, j - 1, k, candidate);
				}
				else
				{
					candidate = m_outputDistance->getValueAt(i, j, k)
							+ sgn(m_outputDistance->getValueAt(i, j - 1, k))* h;
					if (abs(candidate) < abs(m_outputDistance->getValueAt(i, j - 1, k)))
						m_outputDistance->setValueAt(i, j - 1, k, candidate);
				}
			}
		}

		// y-direction forward
		for (int j = intersec_xmin; j < m_outputDistance->getMaxX(); j++)
		{
			for (int i = intersec_ymin; i < m_outputDistance->getMaxY() - 1; i++)
			{
				if (j < intersec_xmax && i < intersec_ymax - 1)
				{
					if (m_inputDistance->getValueAt(i, j, k) * m_inputDistance->getValueAt(i + 1, j, k) <= 0.0)
					{
						// interpolate
						i_slope = (m_inputDistance->getValueAt(i + 1, j, k) - m_inputDistance->getValueAt(i, j, k)) / h;
						distToZero = -m_inputDistance->getValueAt(i, j, k) / i_slope;
						if (abs(m_outputDistance->getValueAt(i, j, k)) > abs(distToZero))
							m_outputDistance->setValueAt(i, j, k,-distToZero * sgn(i_slope));
					}
					// calculate new distance candidate and assign if appropriate
					candidate = m_outputDistance->getValueAt(i, j, k)	+ (sgn(m_inputDistance->getValueAt(i + 1, j,k)) * h);
					if (abs(candidate) < abs(m_outputDistance->getValueAt(i + 1, j, k)))
						m_outputDistance->setValueAt(i + 1, j, k, candidate);
				}
				else
				{
					candidate = m_outputDistance->getValueAt(i, j, k) + (sgn(m_outputDistance->getValueAt(i + 1, j, k))	* h);
					if (abs(candidate) < abs(m_outputDistance->getValueAt(i + 1, j, k)))
						m_outputDistance->setValueAt(i + 1, j, k, candidate);
				}
			}
		}

		for (int j = intersec_xmin; j < m_outputDistance->getMaxX(); j++)
		{
			for (int i = intersec_ymax - 1; i > m_outputDistance->getMinY(); i--)
			{
				if (j < intersec_xmax && i > intersec_ymin)
				{
					// calculate new distance candidate and assign if appropriate
					candidate = m_outputDistance->getValueAt(i, j, k)	+ (sgn(m_inputDistance->getValueAt(i - 1, j, k)) * h);
					if (abs(candidate) < abs(m_outputDistance->getValueAt(i - 1, j, k)))
						m_outputDistance->setValueAt(i - 1, j, k, candidate);
				}
				else
				{
					candidate = m_outputDistance->getValueAt(i, j, k) + (sgn(m_outputDistance->getValueAt(i - 1, j, k)) * h);
					if (abs(candidate) < abs(m_outputDistance->getValueAt(i - 1, j, k)))
						m_outputDistance->setValueAt(i - 1, j, k, candidate);
				}
			}
		}
	}
	// update the point into the third dimensio. the strategy has to change to avoid unneccesary cache loads
	// the idea is to compare all points in one layer first to the next and go on:
	//TODO redist into the third direction
	// z forward:
	for (int k = intersec_zmin; k < m_outputDistance->getMaxZ() - 1; k++)
	{
		for (int i = intersec_ymin; i < m_outputDistance->getMaxY(); i++)
		{
			for (int j = intersec_xmin; j < m_outputDistance->getMaxX(); j++)
			{
				// x-direction forward
				if (k < intersec_zmax - 1 && i < intersec_ymax && j < intersec_xmax)
				{
					if (m_inputDistance->getValueAt(i, j, k) * m_inputDistance->getValueAt(i, j, k + 1) <= 0.0)
					{
						// interpolate
						i_slope = (m_inputDistance->getValueAt(i, j, k + 1) - m_inputDistance->getValueAt(i, j, k)) / h;
						distToZero = -m_inputDistance->getValueAt(i, j, k) / i_slope;
						if (abs(m_outputDistance->getValueAt(i, j, k)) > abs(distToZero))
							m_outputDistance->setValueAt(i, j, k, -distToZero * sgn(i_slope));
					}
					candidate = m_outputDistance->getValueAt(i, j, k) + (sgn(m_inputDistance->getValueAt(i, j, k + 1)) * h);
					if (abs(candidate) < abs(m_outputDistance->getValueAt(i, j, k + 1)))
						m_outputDistance->setValueAt(i, j, k + 1, candidate);
				}
				else
				{
					candidate = m_outputDistance->getValueAt(i, j, k) + (sgn(m_outputDistance->getValueAt(i, j, k + 1))	* h);
					if (abs(candidate) < abs(m_outputDistance->getValueAt(i, j, k + 1)))
						m_outputDistance->setValueAt(i, j, k + 1, candidate);
				}
			}
		}
	}
	// z backward:
	for (int k = intersec_zmax - 1; k > m_outputDistance->getMinZ(); k--)
	{
		for (int i = intersec_ymin; i < m_outputDistance->getMaxY(); i++)
		{
			for (int j = intersec_xmin; j < m_outputDistance->getMaxX(); j++)
			{
				// x-direction forward
				if (k > intersec_zmin && i < intersec_ymax && j < intersec_xmax)
				{
					if (m_inputDistance->getValueAt(i, j, k) * m_inputDistance->getValueAt(i, j, k - 1) <= 0.0)
					{
						// interpolate
						i_slope = (m_inputDistance->getValueAt(i, j, k - 1)- m_inputDistance->getValueAt(i, j, k)) / h;
						distToZero = -m_inputDistance->getValueAt(i, j, k) / i_slope;
						if (abs(m_outputDistance->getValueAt(i, j, k)) > abs(distToZero))
							m_outputDistance->setValueAt(i, j, k,-distToZero * sgn(i_slope));
					}
					candidate = m_outputDistance->getValueAt(i, j, k) + (sgn(m_inputDistance->getValueAt(i, j, k - 1)) * h);
					if (abs(candidate) < abs(m_outputDistance->getValueAt(i, j, k - 1)))
						m_outputDistance->setValueAt(i, j, k - 1, candidate);
				}
				else
				{
					candidate = m_outputDistance->getValueAt(i, j, k) + (sgn(m_outputDistance->getValueAt(i, j, k - 1))	* h);
					if (abs(candidate) < abs(m_outputDistance->getValueAt(i, j, k - 1)))
						m_outputDistance->setValueAt(i, j, k - 1, candidate);
				}
			}
		}
	}

	// for all layers the redist is done - compare into the depth to do
	m_outputDistance->clampValues(-m_grainHandler->delta, m_grainHandler->delta);

	m_inputDistance->resize(m_outputDistance->getMinX(), m_outputDistance->getMinY(),
							m_outputDistance->getMinZ(), m_outputDistance->getMaxX(),
							m_outputDistance->getMaxY(), m_outputDistance->getMaxZ());
}

/**************************************/
// end of redist
/**************************************/

/**************************************/
// plot the box and all its properties
/**************************************/
void LSbox::resizeGrid(int newSize) {

}

void LSbox::recalculateIDLocal() {
	reizeIDLocalToDistanceBuffer();
	executeComparison();
}


double LSbox::computeMisorientation(LSbox* grain_2) {

	double result = (*(m_grainHandler->mymath)).misorientationCubicQxQ(m_orientationQuat[0],
			m_orientationQuat[1], m_orientationQuat[2], m_orientationQuat[3],
			grain_2->m_orientationQuat[0], grain_2->m_orientationQuat[1],
			grain_2->m_orientationQuat[2], grain_2->m_orientationQuat[3]);

	if(result < 3 * PI/180.0)
		result = 3 * PI/180.0;
	return result;
}
double LSbox::computeMisorientation(unsigned int grainID) {
	return computeMisorientation(m_grainHandler->getGrainByID(grainID));
}

double LSbox::GBmobilityModel(double thetaMis) {
	return 1.0;
}

bool LSbox::isNeighbour(LSbox* candidate) {

}


double LSbox::getWeight(int i, int j, bool minimal) {
	return -1;
}

void LSbox::calculateCentroid(SPoint& centroid, vector<GrainJunction> junctions) {

	int nVertices = junctions.size();
	SPoint cent;
	cent.x = 0.0;
	cent.y = 0.0;
	double areaSigned = 0.0;
	double x0 = 0.0; // First vertex' x value
	double y0 = 0.0; // Second vertex' Y value
	double x1 = 0.0; // Next vertex' x value
	double y1 = 0.0; // Next vertex' y value
	double loopArea = 0.0; // Intermediate area

	int i, j = 0;
	for (i = 0, j = nVertices - 1; i < nVertices; j = i++) {

		x0 = junctions[i].coordinates.x;
		y0 = junctions[i].coordinates.y;
		x1 = junctions[j].coordinates.x;
		y1 = junctions[j].coordinates.y;
		loopArea = x0 * y1 - x1 * y0;
		areaSigned += loopArea;
		cent.x += (x0 + x1) * loopArea;
		cent.y += (y0 + y1) * loopArea;
	}

	areaSigned *= 0.5;
	cent.x /= (6.0 * areaSigned);
	cent.y /= (6.0 * areaSigned);

	centroid = cent;
}

double LSbox::getWeigthFromHandler(int i, int j) {
	return m_grainHandler->weightsMatrix[i][j];
}

vector<double> LSbox::linearRegression(vector<SPoint>& points2D) {

	//! the resulting vector contains two elements. The first element is
	//! the slope of the regression line and the second element is the
	//! y-intercept of this line.

	vector<double> linearRegression;
	linearRegression.reserve(2);

	int numberPoints = points2D.size();

	double meanX = 0.0;
	double meanY = 0.0;

	vector<SPoint>::iterator iter;
	for (iter = points2D.begin(); iter != points2D.end(); iter++) {

		meanX += (*iter).x;
		meanY += (*iter).y;
	}

	meanX /= double(numberPoints);
	meanY /= double(numberPoints);

	double covarianceXY = 0.0;
	double varianceX = 0.0;
	for (int i = 0; i < numberPoints; i++) {
		covarianceXY += (points2D[i].x - meanX) * (points2D[i].y - meanY);
		varianceX += (points2D[i].x - meanX) * (points2D[i].x - meanX);
	}

	linearRegression.push_back(covarianceXY / varianceX);
	linearRegression.push_back(meanY - ((covarianceXY / varianceX) * meanX));

	return linearRegression;
}

void LSbox::calculateTriangleCentroid(vector<SPoint>& triangleCentroid,
		vector<SPoint> triangle) {

	int nVertices = triangle.size();
	SPoint cent;
	cent.x = 0.0;
	cent.y = 0.0;
	double areaSigned = 0.0;
	double x0 = 0.0; // First vertex' x value
	double y0 = 0.0; // Second vertex' Y value
	double x1 = 0.0; // Next vertex' x value
	double y1 = 0.0; // Next vertex' y value
	double loopArea = 0.0; // Intermediate area

	int i, j = 0;
	for (i = 0, j = nVertices - 1; i < nVertices; j = i++) {

		x0 = triangle[i].x;
		y0 = triangle[i].y;
		x1 = triangle[j].x;
		y1 = triangle[j].y;
		//cout << "Junction type of i is " << junctions[i].junction_type << " Junction type of j is " << junctions[j].junction_type <<endl;
		loopArea = x0 * y1 - x1 * y0;
		areaSigned += loopArea;
		cent.x += (x0 + x1) * loopArea;
		cent.y += (y0 + y1) * loopArea;
	}

	areaSigned *= 0.5;
	cent.x /= (6.0 * areaSigned);
	cent.y /= (6.0 * areaSigned);

	triangleCentroid.push_back(cent);
}

void LSbox::outputMemoryUsage(ofstream& output) {
	output << m_inputDistance->getTotalMemoryUsed()
			+ m_outputDistance->getTotalMemoryUsed()
			+ m_IDLocal.getTotalMemoryUsed() << endl;

	output << m_outputDistance->getMaxX() - m_outputDistance->getMinX() << " "
			<< m_outputDistance->getMaxY() - m_outputDistance->getMinY() << endl;
}

vector<int> LSbox::getDirectNeighbourIDs() {
	return vector<int>();
}

vector<double> LSbox::getGBLengths() {
	return vector<double>();
}

void LSbox::computeDirectNeighbours(const RTree<unsigned int, int, 3, float>& tree)
{
	int min[3], max[3];
	min[0] = getMinX(); min[1] = getMinY(); min[2] = getMinZ();
	max[0] = getMaxX(); max[1] = getMaxY(); max[2] = getMaxZ();
	vector<unsigned int>	intersectingGrains;
	tree.Search(min, max, intersectingGrains);
	for(unsigned int k=0; k < intersectingGrains.size(); k++)
	{
		if(m_ID != intersectingGrains[k])
		{
			m_secondOrderNeighbours.push_back(intersectingGrains[k]);
			m_explicitBoundary.addDirectNeighbor(intersectingGrains[k]);
		}
	}
}

struct vectorComparator{
	bool operator() (const Vector3d& a, const Vector3d& b) const
	{
			return a[0] < b[0] ? true :
					(a[0] > b[0] ? false :
							(a[1] < b[1] ? true :
									( a[1] > b[1] ? false :
											( a[2] < b[2] ? true : false
											)
									)
							)
					);
	}
};

void LSbox::plotBoxContour(bool absoluteCoordinates)
{
	string filename = string("GrainHull_")+to_string((unsigned long long)m_ID)+string("Timestep_")+
			to_string((unsigned long long)m_grainHandler->loop)+string(".vtk");
	FILE* output = fopen(filename.c_str(), "wt");
	if(output == NULL)
	{
		throw runtime_error("Unable to save box hull!");
	}

	fprintf(output, "%s\n", "# vtk DataFile Version 3.0\n"
							"vtk output\n"
							"ASCII\n"
							"DATASET POLYDATA\n");

	int counter = 0;
	map<Vector3d, int, vectorComparator> mymap;
	map<int, Vector3d> orderedPoints;
	for(unsigned int i=0; i<m_grainHull.size(); i++)
	{
		if(mymap.find(m_grainHull[i].points[0]) == mymap.end())
		{
			mymap.insert(pair<Vector3d, int>(m_grainHull[i].points[0], counter));
			counter++;
		}
		if(mymap.find(m_grainHull[i].points[1]) == mymap.end())
		{
			mymap.insert(pair<Vector3d, int>(m_grainHull[i].points[1], counter));
			counter++;
		}
		if(mymap.find(m_grainHull[i].points[2]) == mymap.end())
		{
			mymap.insert(pair<Vector3d, int>(m_grainHull[i].points[2], counter));
			counter++;
		}
	}
	for ( const auto &myPair : mymap )
	{
		orderedPoints.insert(pair<int, Vector3d>(myPair.second, myPair.first));
	}

	fprintf(output, "POINTS %lu float\n", orderedPoints.size());

	for ( const auto &myPair : orderedPoints )
	{
		fprintf(output, "%f %f %f\n", myPair.second[0], myPair.second[1], myPair.second[2]);
	}

	fprintf(output, "POLYGONS %lu %lu\n", m_grainHull.size(), m_grainHull.size() * 4);
	for(unsigned int i=0; i < m_grainHull.size(); i++)
	{

		fprintf(output, "3 %d %d %d \n", (*(mymap.find(m_grainHull[i].points[2]))).second,
										 (*(mymap.find(m_grainHull[i].points[1]))).second,
										 (*(mymap.find(m_grainHull[i].points[0]))).second);
	}

	fprintf(output, "POINT_DATA %lu\n", orderedPoints.size());
	fprintf(output, "FIELD FieldData 1\n");
	fprintf(output, "Interestingness 1 %lu int\n", orderedPoints.size());

	for ( const auto &myPair : orderedPoints )
	{
		const Vector3d& point = myPair.second;
		int interestingness = 0;
		for(unsigned int i=0; i<m_grainHull.size(); i++)
		{
			if( point == m_grainHull[i].points[0] || point == m_grainHull[i].points[1] ||
				point == m_grainHull[i].points[2])
			{
				interestingness = max(interestingness, m_grainHull[i].additionalData);
			}
		}
		fprintf(output, "%d ", interestingness);
	}
	fclose(output);
}

void LSbox::plotBoxVolumetric(string identifier, E_BUFFER_SELECTION bufferSelection)
{
	string filename = string("GrainVolume_")+to_string((unsigned long long)m_ID)+string("Timestep_")+
				to_string((unsigned long long)m_grainHandler->loop)+identifier+string(".vtk");
	FILE* output = fopen(filename.c_str(), "wt");
	if(output == NULL)
	{
		throw runtime_error("Unable to save box hull!");
	}
	DimensionalBufferReal* distance;
	switch(bufferSelection)
	{
		case E_INPUT_DISTANCE:
			distance = m_inputDistance;
			break;
		case E_OUTPUT_DISTANCE:
			distance = m_outputDistance;
			break;
		default:
			throw runtime_error(string("Invalid buffer selection in plotBoxVolumetric for grain") + to_string((unsigned long long)m_ID)+
					string(" at timestep ") + to_string((unsigned long long)m_grainHandler->loop));
	}

	fprintf(output, "%s\n", "# vtk DataFile Version 3.0\n"
							"vtk output\n"
							"ASCII\n"
							"DATASET STRUCTURED_GRID");
	int totalPoints = 	(distance->getMaxX()-distance->getMinX())*
						(distance->getMaxY()-distance->getMinY())*
						(distance->getMaxZ()-distance->getMinZ());
	fprintf(output, "DIMENSIONS %d %d %d\n",distance->getMaxX()-distance->getMinX(),
											distance->getMaxY()-distance->getMinY(),
											distance->getMaxZ()-distance->getMinZ());
	fprintf(output, "POINTS %d int\n", totalPoints);
	for(int i=distance->getMinY(); i<distance->getMaxY(); i++)
		for(int j=distance->getMinX(); j<distance->getMaxX(); j++)
			for(int k=distance->getMinZ(); k<distance->getMaxZ(); k++)
				fprintf(output, "%d %d %d\n", j,i,k);

	fprintf(output, "\nPOINT_DATA %d\n", totalPoints);
	fprintf(output, "FIELD FieldData 1\n");
	fprintf(output, "Distance 1 %d float\n", totalPoints);
	for(int i=distance->getMinY(); i<distance->getMaxY(); i++)
			for(int j=distance->getMinX(); j<distance->getMaxX(); j++)
				for(int k=distance->getMinZ(); k<distance->getMaxZ(); k++)
					fprintf(output, "%f\n", distance->getValueAt(i,j,k));
	fclose(output);
}

void LSbox::plotBoxIDLocal()
{
	string filename = string("GrainIDLocal_")+to_string((unsigned long long)m_ID)+string("Timestep_")+
				to_string((unsigned long long)m_grainHandler->loop)+string(".vtk");
	FILE* output = fopen(filename.c_str(), "wt");
	if(output == NULL)
	{
		throw runtime_error("Unable to save box hull!");
	}

	fprintf(output, "%s\n", "# vtk DataFile Version 3.0\n"
							"vtk output\n"
							"ASCII\n"
							"DATASET STRUCTURED_GRID");
	int totalPoints = 	(m_IDLocal.getMaxX()-m_IDLocal.getMinX())*
						(m_IDLocal.getMaxY()-m_IDLocal.getMinY())*
						(m_IDLocal.getMaxZ()-m_IDLocal.getMinZ());
	fprintf(output, "DIMENSIONS %d %d %d\n",m_IDLocal.getMaxX()-m_IDLocal.getMinX(),
											m_IDLocal.getMaxY()-m_IDLocal.getMinY(),
											m_IDLocal.getMaxZ()-m_IDLocal.getMinZ());
	fprintf(output, "POINTS %d int\n", totalPoints);
	for(int i=m_IDLocal.getMinY(); i<m_IDLocal.getMaxY(); i++)
		for(int j=m_IDLocal.getMinX(); j<m_IDLocal.getMaxX(); j++)
			for(int k=m_IDLocal.getMinZ(); k<m_IDLocal.getMaxZ(); k++)
				fprintf(output, "%d %d %d\n", j,i,k);

	fprintf(output, "\nPOINT_DATA %d\n", totalPoints);
	fprintf(output, "FIELD FieldData 1\n");
	fprintf(output, "Distance 1 %d float\n", totalPoints);
	for(int i=m_IDLocal.getMinY(); i<m_IDLocal.getMaxY(); i++)
			for(int j=m_IDLocal.getMinX(); j<m_IDLocal.getMaxX(); j++)
				for(int k=m_IDLocal.getMinZ(); k<m_IDLocal.getMaxZ(); k++)
					fprintf(output, "%d\n", m_IDLocal.getValueAt(i,j,k).grainID);
	fclose(output);
}
int LSbox::getNeighbourAt(int i, int j, int k)
{
	if(m_IDLocal.isPointInside(i,j,k))
	{
		return m_IDLocal.getValueAt(i,j,k).grainID;
	}
	else
	{
		return -1;
	}
}
