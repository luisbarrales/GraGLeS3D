#include "grainhdl.h"
#include "Settings.h"
#include <sys/time.h>
#include "rapidxml.hpp"
#include "rapidxml_print.hpp"

using namespace rapidxml;

grainhdl::grainhdl() :
		m_ThreadPoolCount(0)
{}
grainhdl::~grainhdl(){
	delete mymath;
}

void grainhdl::setSimulationParameter(){
	initEnvironment();
	mymath = new mathMethods();
	// 	readInit();
	Mode = (int)Settings::MicrostructureGenMode;
	ngrains = Settings::NumberOfParticles;
	currentNrGrains = ngrains;
	hagb = Settings::HAGB;
	if(Mode==1) realDomainSize= sqrt(ngrains)*Settings::NumberOfPointsPerGrain-1;	// half open container of VORO++
	if(Mode==2 || Mode ==3 ) realDomainSize= sqrt(ngrains)*Settings::NumberOfPointsPerGrain-1;
	discreteEnergyDistribution.resize(Settings::DiscreteSamplingRate);
	fill(discreteEnergyDistribution.begin(),discreteEnergyDistribution.end(),0 );
	
	dt = 0.4/double(realDomainSize*realDomainSize /4);
	h = 1.0/double(realDomainSize);

	delta = Settings::DomainBorderSize * 1/double(realDomainSize);
	tubeRadius = sqrt(2)*1.5*h + 0.001;
	grid_blowup = Settings::DomainBorderSize;
	BoundaryGrainTube=grid_blowup;
	ngridpoints = realDomainSize + (2*grid_blowup); 
	boundary = new LSbox(0, 0, 0, 0, this);
// 	(*boundary).plot_box(false,2,"no.gnu");
	
	switch (Mode) {
		case 1: {
			if(Settings::UseTexture){
				bunge = new double[3]{PI/2, PI/2, PI/2};
				deviation = 15*PI/180;
			}
			else { 
				bunge = NULL; 
				deviation = 0;
			}
			ST = NULL;			
			VOROMicrostructure();
// 			generateRandomEnergy();
			break;
		}
		case 2: {
			bunge = NULL; deviation = 0;
			ST=new double [ngrains*ngrains];
			std::fill_n(ST,ngrains*ngrains,0);
			readMicrostructureFromVertex();
			break;
		}
		case 3:{
			if(Settings::UseTexture){
				bunge = new double[3]{PI/2, PI/2, PI/2};
				deviation = 15*PI/180;
			}
			else {
				bunge = NULL;
				deviation = 0;
			}
			ST = NULL;
			readMicrostructure();
			break;
		}
	}		
// 	construct_boundary();
	//program options:
    cout << endl << "******* PROGRAM OPTIONS: *******" << endl << endl;
    cout << "Number of Grains: " << ngrains << endl;
    cout << "simulated Timesteps: " << Settings::NumberOfTimesteps << endl;
	cout << "DELTA TUBE: " << delta << endl;
    cout << "Timestepwidth " << dt << endl;
    cout << "Number of Gridpoints: " << ngridpoints << endl << endl;
    
    cout << endl << "******* start simulation: *******" << endl << endl;
}



void grainhdl::VOROMicrostructure(){	
	
	stringstream filename, plotfiles;
	int current_cell, cell_id;
	double x,y,z,rx,ry,rz;
	int cell_order[ngrains];
	
	grains.resize(ngrains+1);
	
	vector<LSbox*> local_grains;	

	std::vector<LSbox*>::iterator itg;
	
	// stores the centroids of the cells ; access by (3*Id, 3*ID +1, 3*ID +2)
    part_pos = new double[3*ngrains];
	
	bool randbedingung = false; // bei false ist der container halb offen?! d.h. gitterwert mit 1 werden keinem partikel zugeordnet
    if (randbedingung == false) realDomainSize-=1;	
	
	voronoicell_neighbor c;
	container con(0,1,0,1,0,1,5,5,5,randbedingung,randbedingung,randbedingung,2);
    c_loop_all vl(con);
	
	/**********************************************************/
	// Randomly add particles into the container
	
		for(int i=0;i<ngrains;i++) {
			x=utils::rnd();
			y=utils::rnd();
			z=0;
			con.put(i,x,y,z);
		}
		
	/**********************************************************/
		
    for(int i=0; i < realDomainSize; i++) for(int j= 0; j < realDomainSize; j++){
	x=double(i*h); 
	y=double(j*h); // only point within the domain
    if(con.find_voronoi_cell(x,y,z,rx,ry,rz,cell_id)){
	  cell_id= cell_id++;
	  part_pos[3*(cell_id-1)]=rx;
	  part_pos[3*(cell_id-1)+1]=ry;
	  part_pos[3*(cell_id-1)+2]=rz;

	}
    else fprintf(stderr,"# find_voronoi_cell error for %g %g 0\n",x,y);
    }  
// 	con.draw_cells_gnuplot("particles.gnu");
	int i=0;

	if(vl.start()) 
	do {
		con.compute_cell(c,vl);
		cell_order[ngrains-1-i]=(vl.pid()+1);
		int box_id = vl.pid()+1;
		LSbox* newBox = new LSbox(box_id, c, part_pos,this);
		grains[box_id]= newBox;

	} while(vl.inc());

	delete [] part_pos;
}

void grainhdl::readMicrostructure(){
	FILE * levelset;
	levelset = fopen(Settings::ReadFromFilename.c_str(), "r");
	if ( levelset == NULL )
	{
		cout << "Could not read from specified file !";
		exit(2);
	}
	int id;
	cout << ngrains <<endl;

	double q1, q2, q3, q4, xr, yr, xl, yl;

	grains.resize(ngrains+1);
	int i=0;
	int nvertices;
	double* vertices = new double [1000];

	for(int nn=1; nn<= ngrains; nn++){

		fscanf(levelset, "%d\t %d\t %lf\t %lf\t%lf\t%lf\n", &id, &nvertices, &q1, &q2, &q3, &q4);

		for(unsigned int j=0; j<nvertices; j++){
			fscanf(levelset, "%lf\t %lf\n", &xl, &yl);
			vertices[2*j]   = xl;
			vertices[(2*j)+1] = yl;
		}
		fscanf(levelset, "\n");
		LSbox* newBox = new LSbox(id, nvertices, vertices, q1, q2, q3, q4, this);
		grains[nn]= newBox;
	}
	fclose(levelset);
	delete [] vertices;
}

void grainhdl::readMicrostructureFromVertex(){
	FILE * levelset;
 	levelset = fopen( Settings::ReadFromFilename.c_str(), "r" );

	long id;
	int nedges;
	double phi1, PHI, phi2, xr, yr, xl, yl;
	double* edges;
	
	fscanf(levelset, "%d\n", &ngrains);
	cout << "ngrains : " << ngrains << endl;;
	grains.resize(ngrains+1);
	
	int i=0;
	for(int nn=0; nn< ngrains; nn++){
		
		fscanf(levelset, "%ld\t %d\t %lf\t %lf\t%lf\n", &id, &nedges, &phi1, &PHI, &phi2);
		edges = new double [nedges * 4];
		cout << id << " || " << nedges << " || " << phi1 << " || " << PHI << " || " << phi2<< endl;
		
		for(unsigned int j=0; j<nedges; j++){
			fscanf(levelset, "%lf\t %lf\t %lf\t%lf\n", &xl, &yl, &xr, &yr);	
			cout << xl << " ||\t "<< yl << " ||\t "<< xr << " ||\t "<< yr<< " ||\t " << endl;
			int k = 4*j;
			edges[k]   = xl;
			edges[k+1] = yl;
			edges[k+2] = xr;
			edges[k+3] = yr;
		}
		
		LSbox* newBox = new LSbox(id, nedges, edges, phi1, PHI, phi2, this);
		grains[id]= newBox;
				
	    // calculate distances	    
	    newBox->distancefunctionToEdges(nedges, edges);
		
		delete [] edges;
	}
	
	ST = new double [ngrains*ngrains];			//Create ST array and fill with zeros
	std::fill_n(ST,ngrains*ngrains,0);
	
	for(unsigned int i=0; i<ngrains; i++){
		double buffer;
		fscanf(levelset, "%lf\t", &buffer);		
		for(unsigned int j=0; j<ngrains; j++){
			while(j < i) { 
				fscanf(levelset, "%lf\t", &buffer);
				j++;
			}
			fscanf(levelset, "%lf\t", &buffer);
			ST[j+(ngrains*i)]= (double) buffer;
			ST[i+(ngrains*j)] = ST[j+(ngrains*i)];
		}
		fscanf(levelset, "\n");
	} 
	fclose(levelset);

	for(unsigned int i=0; i<ngrains; i++){
		for(unsigned int j=0; j<ngrains; j++){
			cout << ST[i+(ngrains*j)] << "  \t";
		}
		cout << endl;
	}
}


void grainhdl::distanceInitialisation(){
	for (int i = 1; i < grains.size(); i++){
		grains[i]->distancefunction();
	}
}
 
void grainhdl::convolution(){
	std::vector<LSbox*>::iterator it;
	int i=0;
	for (it = ++grains.begin(); it !=grains.end(); it++,++i){
		if(*it==NULL) continue;
		if((*it)->id==0) (*it)->plot_box(true,2, "error", true);
		(*it)->convolution(m_ThreadMemPool[0]);
	}
}
 


void grainhdl::comparison_box(){
	for (int i = 1; i < grains.size(); i++){
		if(grains[i]==NULL)
			continue;
		grains[i]->comparison(m_ThreadMemPool[0]);
	}
}


void grainhdl::level_set(){
	for (int i = 1; i < grains.size(); i++){
		if(grains[i]==NULL)
			continue;
		if(grains[i]->get_status() == false ) {
			  delete grains[i];
			  removeGrain(i);
		}
		else grains[i]->find_contour();
	}
}


void grainhdl::redistancing(){
	currentNrGrains=0;
	for (int i = 1; i < grains.size(); i++){
		if(grains[i]==NULL)
			continue;
		currentNrGrains +=1;
		grains[i]->redist_box();
	}
}


void grainhdl::save_texture(){
	FILE* myfile;
	FILE* enLenDis;
	stringstream filename;

	int numberGrains=0;
	double totalLength=0;
	double total_energy= 0.0;
	
	filename << "Texture" << "_"<< loop << ".ori";	
	myfile = fopen(filename.str().c_str(), "w");
	
	filename.str("");
	filename << "EnergyLengthDistribution_" << loop<< ".txt";
	enLenDis = fopen(filename.str().c_str(), "w");
	double dh= hagb / (double)Settings::DiscreteSamplingRate;
	double buffer = 0.24;
	double euler[3];
	vector<characteristics> :: iterator it2;
// 	fprintf(myfile, "%d\n", );
	vector<LSbox*> :: iterator it;
	
	std::fill (discreteEnergyDistribution.begin(),discreteEnergyDistribution.end() , 0.0);
	
	for(it = ++grains.begin(); it != grains.end(); it++){
		if(*it!=NULL && (*it)->get_status()==true){
			numberGrains++;
			total_energy += (*it)->energy;

			(*mymath).quaternion2Euler( (*it)->quaternion, euler );
			fprintf(myfile, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", euler[0], euler[1], euler[2], (*it)->volume, (float) (*it)->grainCharacteristics.size(), (float) (*it)->perimeter, (float) (*it)->energy);			

			for(it2=(*it)->grainCharacteristics.begin(); it2!=(*it)->grainCharacteristics.end(); it2++){
				if(!Settings::IsIsotropicNetwork){
				//take into account that every line is twice in the model
				discreteEnergyDistribution[(int)(((*it2).energyDensity)/dh -0.5) ] += 0.5 * (*it2).length;
				}
				totalLength += 0.5 * (*it2).length;
			}
		}
	}
	double sum=0;
	if(!Settings::IsIsotropicNetwork){
		for (int i=0; i < Settings::DiscreteSamplingRate; i++){
				fprintf(enLenDis, "%lf\t%lf\n",(float)(dh*(i+1)),(float)discreteEnergyDistribution[i]);
				printf("%lf\t%lf\n",(float)(dh*(i+1)),(float)discreteEnergyDistribution[i]);
				sum+=(float)discreteEnergyDistribution[i] * (float)(dh*(i+1));
			}
	}
	totalenergy.push_back(0.5*total_energy);
	nr_grains.push_back(numberGrains);
	cout << "Timestep " << loop << " complete:" << endl;
	cout << "Number of grains remaining in the Network :" << numberGrains<< endl;
	cout << "Amount of free Energy in the Network :" << 0.5*total_energy<< "   "<< sum <<endl;
	cout << "Total GB Length in Network :" << totalLength<< endl << endl << endl;
	fclose(myfile);
	fclose(enLenDis);
}
 
  
 
 
void grainhdl::run_sim(){
	distanceInitialisation();
	simulationTime =0;
	find_neighbors();
// 	determineIDs();
	for(loop=Settings::StartTime; loop <= Settings::StartTime+Settings::NumberOfTimesteps; loop++){
		gridCoarsement();
		convolution();
		switchDistancebuffer();
		updateSecondOrderNeighbors();
		comparison_box();
		switchDistancebuffer();
		level_set();
		redistancing();		
		if ( ((loop-Settings::StartTime) % int(Settings::AnalysysTimestep)) == 0 || loop == Settings::NumberOfTimesteps ) {
			saveAllContourEnergies();
			save_texture();
			if(loop == Settings::NumberOfTimesteps) saveMicrostructure();
		}
		simulationTime += dt;
	}
// 	utils::CreateMakeGif();
	cout << "Simulation complete." << endl;
	cout << "Simulation Time: " << simulationTime<< endl;
}  

void grainhdl::saveMicrostructure(){
	stringstream param_xml_name;
	param_xml_name<< "PARAMETERS_SIM_NONAME_TIMESTEP_"<< loop <<"_GRAINS_"<<currentNrGrains<< ".xml";
	stringstream vertex_dump_name;
	vertex_dump_name<< "NETWORK_NONAME_TIMESTEP_"<< loop <<"_GRAINS_"<<currentNrGrains<< ".dat";

	createParamsForSim(param_xml_name.str().c_str(), vertex_dump_name.str().c_str());

	ofstream output;
	output.open(vertex_dump_name.str());
	std::vector<LSbox*>::iterator it;
		for (it = ++grains.begin(); it !=grains.end(); it++){
			if(*it== NULL) continue;
			output << (*it)->id << "\t" << (*it)->contourGrain.size()<< "\t" << (*it)->quaternion[0] << "\t" << (*it)->quaternion[1] << "\t" << (*it)->quaternion[2] << "\t" << (*it)->quaternion[3] << endl;
			(*it)->plot_box_contour(loop, false, &output);
		}
	output.close();
}
void grainhdl::createParamsForSim(const char* param_filename, const char* vertex_dump_filename)
{
	xml_document<> doc_tree;

	xml_node<>* declaration = doc_tree.allocate_node(node_declaration);
	declaration->append_attribute(doc_tree.allocate_attribute("version", "1.0"));
	declaration->append_attribute(doc_tree.allocate_attribute("encoding", "utf-8"));
	doc_tree.append_node(declaration);

	doc_tree.append_node(Settings::generateXMLParametersNode(&doc_tree, vertex_dump_filename,loop, currentNrGrains));
	ofstream output;
	output.open(param_filename);
	output<<doc_tree;
	output.close();

}
 
void grainhdl::save_sim(){
// 	(*my_weights).plot_weightmap(ngridpoints, ID, ST, zeroBox);		
	ofstream myfile;
	myfile.open ("NrGrains&EnergyStatistics.txt");
	for(int i=0; i< nr_grains.size(); i++){
		myfile << nr_grains[i] << "\t";
		myfile << totalenergy[i] << endl;
	}
	myfile.close();

// 	if (SAVEIMAGE)utils::PNGtoGIF("test.mp4");
	//cout << "number of distanzmatrices: "<< domains.size() << endl;
}

void grainhdl::updateSecondOrderNeighbors(){
	for (int i = 1; i < grains.size(); i++){
		if(grains[i]==NULL)
			continue;
		grains[i]->add_n2o_2();
	}
}

void grainhdl::find_neighbors(){
	std::vector<LSbox*>::iterator it,itc;	
	for (it = ++grains.begin(); it !=grains.end(); it++){
		if(*it== NULL) continue;
		for (itc = ++grains.begin(); itc !=grains.end(); itc++){
			if(*itc== NULL) continue;
			if(*it!=*itc) 
				if ((*it)->checkIntersect(*itc))
					(*it)->grainCharacteristics.push_back(characteristics(*itc,0,0,0));
		}
	}				
}

void grainhdl::saveSpecialContourEnergies(int id){
	if (grains[id]==NULL) return;
	grains[id]->plot_box_contour(loop, true);
}

void grainhdl::saveAllContourEnergies(){
	ofstream output;
	stringstream filename;
	filename << "Network_Timestep_"<<loop<<".gnu";
	output.open(filename.str());

	std::vector<LSbox*>::iterator it;
	for (it = ++grains.begin(); it !=grains.end(); it++){
	  if(*it== NULL) continue;
	  (*it)->plot_box_contour(loop, true, &output);
	}
	output.close();
}

void grainhdl::saveAllContourLines(){	
	stringstream filename;
	filename<< "NetworkAtTime_"<< loop << ".gnu";  
	ofstream dateiname;
	dateiname.open(filename.str());
	std::vector<LSbox*>::iterator it;	
	for (it = ++grains.begin(); it !=grains.end(); it++){
		 if(*it== NULL) continue;
		(*it)->plot_box_contour(loop, false);
	}
	dateiname.close();
}
void grainhdl::removeGrain(int id){
    grains[id]=NULL;
}

void grainhdl::switchDistancebuffer(){
	for (int i = 1; i < grains.size(); i++){
		if(grains[i]==NULL)
			continue;
		grains[i]->switchInNOut();
	}
}

void grainhdl::gridCoarsement(){
  if (sqrt(currentNrGrains)*Settings::NumberOfPointsPerGrain/realDomainSize < 0.95 && loop!=0&& Settings::GridCoarsement){
	  double shrink = 1-sqrt(currentNrGrains)*Settings::NumberOfPointsPerGrain/realDomainSize;
	  for (int i = 1; i < grains.size(); i++){
		if(grains[i]==NULL)
			continue;
		  grains[i]->resizeGrid(shrink);
	    }	
	    realDomainSize = realDomainSize * (1-shrink)+1; 
	    ngridpoints = realDomainSize+2*grid_blowup; 
	    h = 1.0/realDomainSize;
	    dt = 1.0/double(realDomainSize*realDomainSize);  
    }
    else {
      switchDistancebuffer();
    }
}
 
void grainhdl::clear_mem() {
	if (ST!=NULL) {delete  [] ST; }
}

void grainhdl::initEnvironment()
{
	//Set up correct Maximum Number of threads
	if(Settings::ExecuteInParallel)
	{
		if (Settings::MaximumNumberOfThreads == 0)
		{
			Settings::MaximumNumberOfThreads = omp_get_max_threads();
		}
		omp_set_num_threads(Settings::MaximumNumberOfThreads);
	}
	else
	{
		Settings::MaximumNumberOfThreads = 1;
	}
	m_ThreadPoolCount = Settings::MaximumNumberOfThreads;
	m_ThreadMemPool.resize(m_ThreadPoolCount);
	for(int i=0; i < m_ThreadMemPool.size(); i++)
	{
		double max_size = Settings::NumberOfPointsPerGrain * Settings::NumberOfPointsPerGrain * 50;
		int power_of_two = 1 << (int)(ceil(log2(max_size))+0.5);
		m_ThreadMemPool[i].resize(power_of_two);
	}

}

void grainhdl::set_h(double hn){
      h =hn;
}
void grainhdl::set_realDomainSize(int realDomainSizen){
     realDomainSize= realDomainSizen;
     ngridpoints = realDomainSize+2*grid_blowup;
}


