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
#include "grainhdl.h"
#include "Settings.h"
#include "spoint.h"
#include "RTree.h"
#include "utilities.h"
#include "box.h"
#include "mymath.h"
#include "voro++/include/voro++/voro++.hh"
#include "rapidxml.hpp"
#include "rapidxml_print.hpp"
#include "Eigen/Dense"

#include <sys/time.h>
#include <stdexcept>

using namespace rapidxml;
using namespace Eigen;
using namespace voro;

grainhdl::grainhdl() :
	m_ThreadPoolCount(0) {
}
grainhdl::~grainhdl() {
	delete mymath;
}

void grainhdl::initializeSimulation() {
	initEnvironment();
	mymath = new mathMethods();
	// 	readInit();

	//! The NumberOfParticles passed via parameters.xml is altered
	//! for MicrostructureGenMode 4
	ngrains = Settings::NumberOfParticles;
	currentNrGrains = ngrains;
	realDomainSize = pow(Settings::NumberOfParticles, 1 / 3.0)
			* Settings::NumberOfPointsPerGrain; // half open container of VORO++
	discreteEnergyDistribution.resize(Settings::DiscreteSamplingRate);
	fill(discreteEnergyDistribution.begin(), discreteEnergyDistribution.end(),
			0);

	switch (Settings::ConvolutionMode) {
	case E_LAPLACE: {
		dt = 0.8 / double(realDomainSize * realDomainSize * realDomainSize);
		break;
	}
	case E_LAPLACE_RITCHARDSON: {
		dt = 0.8 / double(realDomainSize * realDomainSize * realDomainSize);
		break;
	}
	case E_GAUSSIAN: {
		dt = 1. / double(realDomainSize * realDomainSize * realDomainSize);
		break;
	}
	default: {
		throw std::runtime_error("Unknown convolution mode!");
	}
	}
	h = 1.0 / double(realDomainSize);

	//Recalculate the setting parameters for the sector radiuses
	Settings::ConstantSectorRadius *= h;
	Settings::InterpolatingSectorRadius *= h;

	delta = Settings::DomainBorderSize * 1 / double(realDomainSize);
	grid_blowup = Settings::DomainBorderSize;
	BoundaryGrainTube = grid_blowup;
	ngridpoints = realDomainSize + (2 * grid_blowup);
	boundary = new LSbox(0, 0, 0, 0, this);
	// 	(*boundary).plot_box(false,2,"no.gnu");
	//!grains.resize(Settings::NumberOfParticles + 1);
	grains.resize(Settings::NumberOfParticles + 1);

	switch (Settings::MicrostructureGenMode) {
	case E_GENERATE_WITH_VORONOY: {
		if (Settings::UseTexture) {
			bunge = new double[3] { PI / 2, PI / 2, PI / 2 };
			deviation = 15 * PI / 180;
		} else {
			bunge = NULL;
			deviation = 0;
		}
		ST = NULL;
		VOROMicrostructure();
		break;
	}
	default: {
		throw std::runtime_error("Unknown microstructure generation mode!");
	}
	}

	// 	construct_boundary();
	//program options:
	cout << endl << "******* PROGRAM OPTIONS: *******" << endl << endl;
	cout << "Number of Grains: " << ngrains << endl;
	cout << "simulated Timesteps: " << Settings::NumberOfTimesteps << endl;
	cout << "DELTA TUBE: " << delta << endl;
	cout << "Timestepwidth " << dt << endl;
	cout << "Number of Gridpoints: " << realDomainSize << " in [0,1]" << endl
			<< endl;

	cout << endl << "******* start simulation: *******" << endl << endl;
}

void grainhdl::VOROMicrostructure() {

	stringstream filename, plotfiles;

	grains.resize(Settings::NumberOfParticles + 1);

	bool randbedingung = false; // bei false ist der container halb offen?! d.h. gitterwert mit 1 werden keinem partikel zugeordnet
	if (randbedingung == false)
		realDomainSize -= 1;

	voronoicell_neighbor c;
	container con(0, 1, 0, 1, 0, 1, 5, 5, 5, randbedingung, randbedingung,
			randbedingung, 2);
	c_loop_all vl(con);

	/**********************************************************/
	// Randomly add particles into the container
	for (int i = 0; i < ngrains; i++) {
		double x = rnd();
		double y = rnd();
		double z = rnd();
		con.put(i, x, y, z);
	}

	/**********************************************************/

	vector<vector<Vector3d> > initialHulls;
	vector<double> cellCoordinates;
	if (vl.start()) {
		initialHulls.resize(ngrains);
		int cellIndex = 0;
		do {
			double cur_x, cur_y, cur_z;
			con.compute_cell(c, vl);
			vl.pos(cur_x, cur_y, cur_z);
			c.vertices(cur_x, cur_y, cur_z, cellCoordinates);
			for (unsigned int i = 0; i < cellCoordinates.size() / 3; i++) {
				initialHulls.at(cellIndex).push_back(
						Vector3d(cellCoordinates.at(3 * i + 1),
								cellCoordinates.at(3 * i),
								cellCoordinates.at(3 * i + 2)));
			}
			cellIndex++;
		} while (vl.inc());

		IDField.resize(0, 0, 0, ngridpoints, ngridpoints, ngridpoints);
		double x, y, z, rx, ry, rz;
		int cell_id;
		for (int k = 0; k < ngridpoints; k++) {
			for (int i = 0; i < ngridpoints; i++) {
				for (int j = 0; j < ngridpoints; j++) {
					x = double((j - grid_blowup) * h);
					y = double((i - grid_blowup) * h);
					z = double((k - grid_blowup) * h);

					if (i < grid_blowup || j < grid_blowup || k < grid_blowup
							|| i >= ngridpoints - grid_blowup || j
							>= ngridpoints - grid_blowup || k >= ngridpoints
							- grid_blowup) {
						IDField.setValueAt(i, j, k, 0);
					} else if (con.find_voronoi_cell(x, y, z, rx, ry, rz,
							cell_id)) {
						int box_id = cell_id + 1;
						IDField.setValueAt(i, j, k, box_id);
					} else {
						IDField.setValueAt(i, j, k, 0);
					}
				}
			}
		}
	} else {
		throw runtime_error("Voronoy container error at start() method!");
	}

	buildBoxVectors(initialHulls);
}

void grainhdl::readMicrostructure() {
}

void grainhdl::readMicrostructureFromVertex() {
}

void grainhdl::distanceInitialisation() {
	unsigned int j;
#pragma omp parallel for private(j)
	for (unsigned int i = 0; i < Settings::MaximumNumberOfThreads; i++)
		for (j = 0; j < Settings::NumberOfParticles
				/ Settings::MaximumNumberOfThreads + 1; j++) {
			if (j * Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num()
					< grains.size()) {
				int id = j * Settings::MaximumNumberOfThreads + 1
						+ omp_get_thread_num();
				if (grains[id] != NULL)
					grains[id]->calculateDistanceFunction(IDField);
			}
		}
}

void grainhdl::convolution(double& planOverhead) {
	unsigned int j;
	double timer = 0;
	timeval time;
#pragma omp parallel for private(j)
	for (unsigned int i = 0; i < Settings::MaximumNumberOfThreads; i++)
		for (j = 0; j < Settings::NumberOfParticles
				/ Settings::MaximumNumberOfThreads + 1; j++) {
			if (j * Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num()
					< grains.size()) {
				int id = j * Settings::MaximumNumberOfThreads + 1
						+ omp_get_thread_num();
				if (grains[id] != NULL)
					grains[id]->initConvoMemory(
							m_ThreadMemPool[omp_get_thread_num()]);
			}
		}
	gettimeofday(&time, NULL);
	timer = time.tv_sec + time.tv_usec / 1000000.0;
	for (unsigned int i = 1; i < grains.size(); i++) {
		if (grains[i] != NULL)
			grains[i]->createConvolutionPlans(
					m_ThreadMemPool[(i - 1) % Settings::MaximumNumberOfThreads]);
	}
	gettimeofday(&time, NULL);
	planOverhead += time.tv_sec + time.tv_usec / 1000000.0 - timer;

#pragma omp parallel for private(j)
	for (unsigned int i = 0; i < Settings::MaximumNumberOfThreads; i++)
		for (j = 0; j < Settings::NumberOfParticles
				/ Settings::MaximumNumberOfThreads + 1; j++) {
			if (j * Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num()
					< grains.size()) {
				int id = j * Settings::MaximumNumberOfThreads + 1
						+ omp_get_thread_num();
				if (grains[id] != NULL)
					grains[id]->executeConvolution(
							m_ThreadMemPool[omp_get_thread_num()]);
			}
		}

	gettimeofday(&time, NULL);
	timer = time.tv_sec + time.tv_usec / 1000000.0;
	for (unsigned int i = 1; i < grains.size(); i++) {
		if (grains[i] != NULL)
			grains[i]->cleanupConvolution();
	}
	gettimeofday(&time, NULL);
	planOverhead += time.tv_sec + time.tv_usec / 1000000.0 - timer;
}

void grainhdl::comparison_box() {
	unsigned int j;
#pragma omp parallel for private(j)
	for (unsigned int i = 0; i < Settings::MaximumNumberOfThreads; i++)
		for (j = 0; j < Settings::NumberOfParticles
				/ Settings::MaximumNumberOfThreads + 1; j++) {
			if (j * Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num()
					< grains.size()) {
				int id = j * Settings::MaximumNumberOfThreads + 1
						+ omp_get_thread_num();
				if (grains[id] != NULL) {
					grains[id]->executeComparison();
					grains[id]->executeSetComparison();
				}
			}
		}
}

void grainhdl::level_set() {
	unsigned int j;
#pragma omp parallel for private(j)
	for (unsigned int i = 0; i < Settings::MaximumNumberOfThreads; i++)
		for (j = 0; j < Settings::NumberOfParticles
				/ Settings::MaximumNumberOfThreads + 1; j++) {
			if (j * Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num()
					< grains.size()) {
				int id = j * Settings::MaximumNumberOfThreads + 1
						+ omp_get_thread_num();
				if (grains[id] == NULL)
					continue;
				if (grains[id]->grainExists() == false) {
					delete grains[id];
					grains[id] = NULL;
				} else
					grains[id]->extractContour();
			}
		}
}

void grainhdl::redistancing() {
	unsigned int j;
#pragma omp parallel for private(j)
	for (unsigned int i = 0; i < Settings::MaximumNumberOfThreads; i++)
		for (j = 0; j < Settings::NumberOfParticles
				/ Settings::MaximumNumberOfThreads + 1; j++) {
			if (j * Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num()
					< grains.size()) {
				int id = j * Settings::MaximumNumberOfThreads + 1
						+ omp_get_thread_num();
				if (grains[id] == NULL)
					continue;
				grains[id]->executeRedistancing();
			}
		}
}

void grainhdl::save_texture() {
	double totalSurface = 0;
	double total_energy = 0.0;
	string filename = string("Texture_") + to_string((unsigned long long) loop)
			+ string(".ori");
	ofstream file;
	file.open(filename.c_str());
	double *euler = new double[3];
	for (int i = 1; i < grains.size(); i++) {
		if (grains[i] != NULL && grains[i]->grainExists()) {
			totalSurface += grains[i]->getSurface() * 0.5;
			total_energy += grains[i]->getEnergy() * 0.5;
			;
			euler = grains[i]->getOrientationQuat()->Quaternion2EulerConst();
			file << grains[i]->getID() << " "
					<< grains[i]->getDirectNeighbourCount() << " "
					<< grains[i]->intersectsBoundaryGrain() << " "
					<< grains[i]->getVolume() << " " << 0 << " "
					<< grains[i]->getSurface() << " " << grains[i]->getEnergy()
					<< " " << euler[0] << " " << euler[1] << " " << euler[2]
					<< "\n";

		}
	}
	file.close();
	nr_grains.push_back(currentNrGrains);
	time.push_back(simulationTime);
	totalenergy.push_back(total_energy);
	cout << "Timestep " << loop << " complete:" << endl;
	cout << "Number of grains remaining in the Network :" << currentNrGrains
			<< endl;
	cout << "Amount of free Energy in the Network :" << total_energy << endl;
	cout << "Total GB Surface in the Network :" << totalSurface << endl << endl
			<< endl;
}

void grainhdl::run_sim() {
	double parallelRest = 0;
	double convo_time = 0;
	double comparison_time = 0;
	double levelset_time = 0;
	double redistancing_time = 0;
	double plan_overhead = 0;

	timeval time;
	double timer;
	distanceInitialisation();
	simulationTime = 0;
	find_neighbors();
	for (loop = Settings::StartTime; loop <= Settings::StartTime
			+ Settings::NumberOfTimesteps; loop++) {
		gettimeofday(&time, NULL);
		timer = time.tv_sec + time.tv_usec / 1000000.0;
		gridCoarsement();
		gettimeofday(&time, NULL);
		parallelRest += time.tv_sec + time.tv_usec / 1000000.0 - timer;

		gettimeofday(&time, NULL);
		timer = time.tv_sec + time.tv_usec / 1000000.0;
		convolution(plan_overhead);
		gettimeofday(&time, NULL);
		convo_time += time.tv_sec + time.tv_usec / 1000000.0 - timer;

		gettimeofday(&time, NULL);
		timer = time.tv_sec + time.tv_usec / 1000000.0;
		switchDistancebuffer();
		updateSecondOrderNeighbors();
		gettimeofday(&time, NULL);
		parallelRest += time.tv_sec + time.tv_usec / 1000000.0 - timer;

		gettimeofday(&time, NULL);
		timer = time.tv_sec + time.tv_usec / 1000000.0;
		comparison_box();
		gettimeofday(&time, NULL);
		comparison_time += time.tv_sec + time.tv_usec / 1000000.0 - timer;

		gettimeofday(&time, NULL);
		timer = time.tv_sec + time.tv_usec / 1000000.0;
		switchDistancebuffer();
		gettimeofday(&time, NULL);
		parallelRest += time.tv_sec + time.tv_usec / 1000000.0 - timer;

		gettimeofday(&time, NULL);
		timer = time.tv_sec + time.tv_usec / 1000000.0;
		level_set();
		gettimeofday(&time, NULL);
		levelset_time += time.tv_sec + time.tv_usec / 1000000.0 - timer;

		gettimeofday(&time, NULL);
		timer = time.tv_sec + time.tv_usec / 1000000.0;
		redistancing();
		gettimeofday(&time, NULL);
		redistancing_time += time.tv_sec + time.tv_usec / 1000000.0 - timer;

		countGrains();
		if (((loop - Settings::StartTime) % int(Settings::AnalysisTimestep))
				== 0 || loop == Settings::NumberOfTimesteps) {
			saveNetworkState();
			save_texture();
			save_sim();
for		(const auto & it : grains)
		{	if (it!= NULL)
			if(it->getID()!=0) {
				//it->plotBoxContour();
				//it->plotBoxVolumetric("",E_OUTPUT_DISTANCE);
			}
		}
	}
	simulationTime += dt;
	if (currentNrGrains < Settings::BreakupNumber) {
		cout << "Network has coarsed to less than specified by Settings::BreakupNumber. "
		<< "Remaining Grains: " << currentNrGrains
		<< ". Break and save." << endl;
		break;
	}
}
// 	utils::CreateMakeGif();

cout << "Simulation complete." << endl;
cout << "Simulation Time: " << simulationTime << endl;
cout << "Detailed timings: " << endl;
cout << "Convolution time: " << convo_time << endl;
cout << "     Of which plan overhead is: " << plan_overhead << endl;
cout << "Comparison time: " << comparison_time << endl;
cout << "Redistancing time: " << redistancing_time << endl;
cout << "Levelset time: " << levelset_time << endl;
cout << "GridCoarse/SwitchBuffer/UpNeigh: " << parallelRest << endl;
cout << "Sum parallel regions: " << convo_time + comparison_time
+ levelset_time + parallelRest + redistancing_time << endl;
}

void grainhdl::countGrains() {
	currentNrGrains = 0;
	for (auto it = ++grains.begin(); it != grains.end(); it++) {
		if (*it == NULL)
			continue;
		currentNrGrains++;
	}
}

void grainhdl::saveMicrostructure() {
}
void grainhdl::createParamsForSim(const char* param_filename,
		const char* vertex_dump_filename) {
	xml_document<> doc_tree;

	xml_node<>* declaration = doc_tree.allocate_node(node_declaration);
	declaration->append_attribute(doc_tree.allocate_attribute("version", "1.0"));
	declaration->append_attribute(
			doc_tree.allocate_attribute("encoding", "utf-8"));
	doc_tree.append_node(declaration);

	doc_tree.append_node(
			Settings::generateXMLParametersNode(&doc_tree,
					vertex_dump_filename, loop, currentNrGrains));
	ofstream output;
	output.open(param_filename);
	output << doc_tree;
	output.close();

}

void grainhdl::save_sim() {
	ofstream myfile;
	myfile.open("NrGrains&EnergyStatistics.txt");
	for (unsigned int i = 0; i < nr_grains.size(); i++) {
		myfile << time[i] << "\t";
		myfile << nr_grains[i] << "\t";
		//myfile << totalenergy[i] << "\t";
		myfile << realDomainSize << endl;
	}
	myfile.close();

}

void grainhdl::updateSecondOrderNeighbors() {
	unsigned int j;
#pragma omp parallel for private(j)
	for (unsigned int i = 0; i < Settings::MaximumNumberOfThreads; i++)
		for (j = 0; j < Settings::NumberOfParticles
				/ Settings::MaximumNumberOfThreads + 1; j++) {
			if (j * Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num()
					< grains.size()) {
				int id = j * Settings::MaximumNumberOfThreads + 1
						+ omp_get_thread_num();
				if (grains[id] == NULL)
					continue;
				grains[id]->computeSecondOrderNeighbours();
			}
		}
}

void grainhdl::find_neighbors() {
	RTree<unsigned int, int, 3, float> tree;
	int min[3], max[3];
	for (unsigned int i = 1; i <= Settings::NumberOfParticles; i++) {
		if (grains[i] == NULL)
			continue;
		min[0] = grains[i]->getMinX();
		min[1] = grains[i]->getMinY();
		min[2] = grains[i]->getMinZ();
		max[0] = grains[i]->getMaxX();
		max[1] = grains[i]->getMaxY();
		max[2] = grains[i]->getMaxZ();
		tree.Insert(min, max, i);
	}

	unsigned int j;
#pragma omp parallel for private(j)
	for (unsigned int i = 0; i < Settings::MaximumNumberOfThreads; i++)
		for (j = 0; j < Settings::NumberOfParticles
				/ Settings::MaximumNumberOfThreads + 1; j++) {
			if (j * Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num()
					< grains.size()) {
				int id = j * Settings::MaximumNumberOfThreads + 1
						+ omp_get_thread_num();
				if (grains[id] == NULL)
					continue;
				grains[id]->computeDirectNeighbours(tree);
			}
		}
}

void grainhdl::saveSpecialContourEnergies(int id) {
}

void grainhdl::saveNetworkState() {
	unsigned int j;
#pragma omp parallel for private(j)
	for (unsigned int i = 0; i < Settings::MaximumNumberOfThreads; i++)
		for (j = 0; j < Settings::NumberOfParticles
				/ Settings::MaximumNumberOfThreads + 1; j++) {
			if (j * Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num()
					< grains.size()) {
				int id = j * Settings::MaximumNumberOfThreads + 1
						+ omp_get_thread_num();
				if (grains[id] == NULL)
					continue;
			}
		}
}

void grainhdl::save_id() {
}

void grainhdl::saveAllContourLines() {

}
void grainhdl::removeGrain(int id) {
	grains[id] = NULL;
}

void grainhdl::switchDistancebuffer() {
	unsigned int j;
#pragma omp parallel for private(j)
	for (unsigned int i = 0; i < Settings::MaximumNumberOfThreads; i++)
		for (j = 0; j < Settings::NumberOfParticles
				/ Settings::MaximumNumberOfThreads + 1; j++) {
			if (j * Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num()
					< grains.size()) {
				int id = j * Settings::MaximumNumberOfThreads + 1
						+ omp_get_thread_num();
				if (grains[id] == NULL)
					continue;
				grains[id]->switchInNOut();
			}
		}
}

void grainhdl::clear_mem() {
	if (ST != NULL) {
		delete[] ST;
	}
}

void grainhdl::initEnvironment() {
	//Set up correct Maximum Number of threads
	if (Settings::ExecuteInParallel) {
		Settings::MaximumNumberOfThreads = omp_get_max_threads();

	} else {
		Settings::MaximumNumberOfThreads = 1;
		omp_set_num_threads(Settings::MaximumNumberOfThreads);
	}

	m_ThreadPoolCount = Settings::MaximumNumberOfThreads;
	m_ThreadMemPool.resize(m_ThreadPoolCount);
	initNUMABindings();

#pragma omp parallel
	{
		double max_size = Settings::NumberOfPointsPerGrain
				* Settings::NumberOfPointsPerGrain * 50;
		int power_of_two = 1 << (int) (ceil(log2(max_size)) + 0.5);
		//!int power_of_two = 1 << (int) (ceil(log2(2<<20)) + 0.5); //!27
		m_ThreadMemPool[omp_get_thread_num()].resize(power_of_two);
	}
}

void grainhdl::buildBoxVectors(vector<vector<Vector3d>>& hulls) {
	unsigned int j;
	bool exceptionHappened = false;
	string error_message;

#pragma omp parallel for private(j)
	for (unsigned int i = 0; i < Settings::MaximumNumberOfThreads; i++)
		for (j = 0; j < Settings::NumberOfParticles
				/ Settings::MaximumNumberOfThreads + 1; j++) {
			unsigned int id = j * Settings::MaximumNumberOfThreads + 1
					+ omp_get_thread_num();
			if (id < grains.size()) {
				try {
					LSbox* grain = new LSbox(id, IDField, this);
					grains[id] = grain;
				} catch (exception& e) {
#pragma omp critical
					{
						exceptionHappened = true;
						error_message += string("Grain ") + to_string(
								(unsigned long long) id) + string(
								" failed at timestep ") + to_string(
								(unsigned long long) loop)
								+ " in its constructor! Reason : " + e.what()
								+ string("\n");
					}
				}
			}
		}
	if (exceptionHappened) {
		throw runtime_error(error_message);
	}
}

void grainhdl::set_h(double hn) {
	h = hn;
}
void grainhdl::set_realDomainSize(int realDomainSizen) {
	realDomainSize = realDomainSizen;
	ngridpoints = realDomainSize + 2 * grid_blowup;
}
/**
 * This function analyzes the input file for MicrostructureGenMode 4.
 * The amount of lines of the input file is determined. This number
 * indicates the number of points specified in the file, i.e the number of
 * grains. This means one pair of x-y-coordinates each in every line. A
 * tabulator is used as the separator.
 *
 * @return the amount of points in the input file
 */
int grainhdl::read_ScenarioPoints() {

	int counter = 0;
	string line;
	ifstream reader(Settings::ReadFromFilename.c_str());
	while (std::getline(reader, line)) {
		counter++;
	}
	return counter;
}

struct NUMANode {
	int num_cpus;
	int numa_cpus[64];
};

void grainhdl::initNUMABindings() {
	vector<NUMANode> nodes;
	nodes.reserve(16);
	numa_available();
	// returns a mask of CPUs on which the current task is allowed to run.
	bitmask* mask = numa_get_run_node_mask();
	bitmask* cpus = numa_allocate_cpumask();
	for (unsigned int j = 0; j < mask->size; j++) {
		if (numa_bitmask_isbitset(mask, j)) {
			printf("We are allowed to used node %d\n", j);
			NUMANode node;
			memset(&node, 0xFF, sizeof(node));
			//converts a node number to a bitmask of CPUs.
			//The user must pass a bitmask structure with a mask buffer long enough to represent all possible cpu's
			numa_node_to_cpus(j, cpus);
			node.num_cpus = numa_bitmask_weight(cpus);
			int cpuCounter = 0;
			for (unsigned int i = 0; i < cpus->size; i++) {

				if (numa_bitmask_isbitset(cpus, i) && numa_bitmask_isbitset(
						numa_all_cpus_ptr, i)) {
					node.numa_cpus[cpuCounter] = i;
					cpuCounter++;
				}
			}
			nodes.push_back(node);
		}
	}
	numa_free_cpumask(cpus);
#pragma omp parallel
	{
		int threadID = omp_get_thread_num();
		for (unsigned int i = 0; i < nodes.size(); i++) {
			if (threadID < nodes.at(i).num_cpus) {
#pragma omp critical
				{
					printf("Will bind thread %d to cpu %d\n",
							omp_get_thread_num(),
							nodes.at(i).numa_cpus[threadID]);
					cpu_set_t set;
					CPU_ZERO(&set);
					CPU_SET(nodes.at(i).numa_cpus[threadID], &set);
					int res = sched_setaffinity(0, sizeof(set), &set);
					printf(res == 0 ? "Successful\n" : "Failed\n");
				}
				break;
			}
			threadID -= nodes.at(i).num_cpus;
		}
	}
}

void grainhdl::gridCoarsement() {
	if ((double) currentNrGrains / (double) ngrains
			< Settings::GridCoarsementGradient && loop != 0
			&& Settings::GridCoarsement) {
		int newSize = pow(currentNrGrains, 1 / 3.0)
				* Settings::NumberOfPointsPerGrain;
		cout << "coarsing the current grid in Timestep: " << loop << endl;
		cout << "newSize :" << newSize << endl << endl;
#pragma omp parallel
		{
			for (unsigned int j = 0; j < Settings::NumberOfParticles
					/ Settings::MaximumNumberOfThreads + 1; j++) {
				if (j * Settings::MaximumNumberOfThreads + 1
						+ omp_get_thread_num() < grains.size()) {
					int id = j * Settings::MaximumNumberOfThreads + 1
							+ omp_get_thread_num();
					if (grains[id] == NULL)
						continue;
					grains[id]->resizeGrid(newSize);
				}

			}
		}
		realDomainSize = newSize;
		delta = Settings::DomainBorderSize * 1 / double(realDomainSize);
		ngridpoints = realDomainSize + 2 * grid_blowup;
		h = 1.0 / realDomainSize;
		//! DISCREPANCY: Compare to the application of dt in the convolution, time decreasing factor 0.8
		switch (Settings::ConvolutionMode) {
		case E_LAPLACE: {
			dt = 0.8 / double(realDomainSize * realDomainSize * realDomainSize);
			break;
		}
		case E_LAPLACE_RITCHARDSON: {
			dt = 0.8 / double(realDomainSize * realDomainSize * realDomainSize);
			break;
		}
		case E_GAUSSIAN: {
			dt = 1. / double(realDomainSize * realDomainSize * realDomainSize);
			break;
		}
		default: {
			throw std::runtime_error("Unknown convolution mode!");
		}
		}
		ngrains = currentNrGrains;
#pragma omp parallel
		{
			for (unsigned int j = 0; j < Settings::NumberOfParticles
					/ Settings::MaximumNumberOfThreads + 1; j++) {
				if (j * Settings::MaximumNumberOfThreads + 1
						+ omp_get_thread_num() < grains.size()) {
					int id = j * Settings::MaximumNumberOfThreads + 1
							+ omp_get_thread_num();
					if (grains[id] == NULL)
						continue;
					grains[id]->recalculateIDLocal();
				}
			}
		}
#pragma omp parallel
		{
			for (unsigned int j = 0; j < Settings::NumberOfParticles
					/ Settings::MaximumNumberOfThreads + 1; j++) {
				if (j * Settings::MaximumNumberOfThreads + 1
						+ omp_get_thread_num() < grains.size()) {
					int id = j * Settings::MaximumNumberOfThreads + 1
							+ omp_get_thread_num();
					if (grains[id] == NULL)
						continue;
					grains[id]->extractContour();
				}
			}
		}

	} else {
		switchDistancebuffer();
	}
}
