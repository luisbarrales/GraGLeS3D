#include "grainhdl.h"
#include "Settings.h"
#include <sys/time.h>
#include "rapidxml.hpp"
#include "rapidxml_print.hpp"

using namespace rapidxml;

grainhdl::grainhdl(){}
grainhdl::~grainhdl(){
	delete mymath;
}



void grainhdl::setSimulationParameter(){


	mymath = new mathMethods();
	// 	readInit();
	Mode = (int)Settings::MicrostructureGenMode;
	ngrains = Settings::NumberOfParticles;
	
	hagb = Settings::HAGB;
	if(Mode==1) realDomainSize= sqrt(ngrains)*Settings::NumberOfPointsPerGrain-1;	// half open container of VORO++
	if(Mode==2) realDomainSize= sqrt(ngrains)*Settings::NumberOfPointsPerGrain;
	discreteEnergyDistribution.resize(Settings::DiscreteSamplingRate);
	dt = 1.0/double(realDomainSize*realDomainSize);
	h = 1.0/double(realDomainSize);
	delta = Settings::DomainBorderSize * 1/double(realDomainSize);
	tubeRadius = sqrt(2)*1.5*h + 0.001;
	grid_blowup = Settings::DomainBorderSize;
	
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
	// iteration over all cells in the container con:
	if(vl.start()) 
	do {
	  // compute the current cell, taken out of the container
		con.compute_cell(c,vl);
		cell_order[ngrains-1-i]=(vl.pid()+1);
		
		// create a new Box for the current cell
		int box_id = vl.pid()+1;
		LSbox* newBox = new LSbox(box_id, c, part_pos,this);
		
		grains[box_id]= newBox;
		newBox->distancefunction(c, part_pos);        

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
	int nvertex;
	double q1, q2, q3, q4, xr, yr, xl, yl;
	double* vertices;

	fscanf(levelset, "%d\n", &ngrains);
	cout << "ngrains : " << ngrains << endl;;
	grains.resize(ngrains+1);

	int i=0;
	for(int nn=0; nn< ngrains; nn++){
		fscanf(levelset, "%d\t %d\t %lf\t %lf\t%lf\t%lf\n", &id, &nvertex, &q1, &q2, &q3, &q4);
		vertices = new double [nvertex * 4];
		cout << id << " || " << nvertex << " || " << q1 << " || " << q2 << " || " << q3<< " || " << q4 << endl;

		for(unsigned int j=0; j<nvertex; j++){
			fscanf(levelset, "%lf\t %lf\t %lf\t%lf\n", &xl, &yl, &xr, &yr);
			cout << xl << " ||\t "<< yl << " ||\t "<< xr << " ||\t "<< yr<< " ||\t " << endl;
			int k = 4*j;
			vertices[k]   = xl;
			vertices[k+1] = yl;
			vertices[k+2] = xr;
			vertices[k+3] = yr;
		}

		LSbox* newBox = new LSbox(id, nvertex, vertices, q1, q2, q3, q4, this);
		grains[id]= newBox;

		// calculate distances
		newBox->distancefunction(nvertex, vertices);

		delete [] vertices;
	}
}

void grainhdl::readMicrostructureFromVertex(){
	FILE * levelset;
 	levelset = fopen( "lsInput_DRAG.dat", "r" );
// 	levelset = fopen( "lsInput_quadrat.dat", "r" );

	long id;
	int nvertex;
	double phi1, PHI, phi2, xr, yr, xl, yl;
	double* vertices;
	
	fscanf(levelset, "%d\n", &ngrains);
	cout << "ngrains : " << ngrains << endl;;
	grains.resize(ngrains+1);
	
	int i=0;
	for(int nn=0; nn< ngrains; nn++){
		
		fscanf(levelset, "%ld\t %d\t %lf\t %lf\t%lf\n", &id, &nvertex, &phi1, &PHI, &phi2);
		vertices = new double [nvertex * 4];
		cout << id << " || " << nvertex << " || " << phi1 << " || " << PHI << " || " << phi2<< endl;
		
		for(unsigned int j=0; j<nvertex; j++){
			fscanf(levelset, "%lf\t %lf\t %lf\t%lf\n", &xl, &yl, &xr, &yr);	
			cout << xl << " ||\t "<< yl << " ||\t "<< xr << " ||\t "<< yr<< " ||\t " << endl;
			int k = 4*j;
			vertices[k]   = xl;
			vertices[k+1] = yl;
			vertices[k+2] = xr;
			vertices[k+3] = yr;
		}
		
		LSbox* newBox = new LSbox(id, nvertex, vertices, phi1, PHI, phi2, this);
		grains[id]= newBox;
				
	    // calculate distances	    
	    newBox->distancefunction(nvertex, vertices); 
		
		delete [] vertices;
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
// 			cout << "buffer " << buffer <<endl ;
// 			fwrite(/*levelset*/, "test");
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
 

 
// void grainhdl::generateRandomEnergy(){	
// 	const double MIN = 0.6;
// 	const double MAX = 1.5;
// 	for(int i=0; i < ngrains; i++){
// 		for(int j=0; j <=i; j++){			
// 			double zahl=(double)(rand() / (((double)RAND_MAX+1)/ (double)(MAX-MIN)))+MIN;
// // 			if (i==6 && j ==5) { 
// // 				ST[i+(PARTICLES*j)] = 0.1;
// // 				ST[j+(PARTICLES*i)] = 0.1;
// // 			}
// 			ST[i+(ngrains*j)] = zahl;
// 			ST[j+(ngrains*i)] = zahl;
// 			if(i==j) ST[j+(ngrains*i)] = 1.0;
// 		}
// 	} 
// }
//  
 
void grainhdl::convolution(){
	std::vector<LSbox*>::iterator it;
	for (it = ++grains.begin(); it !=grains.end(); it++){	
		if(*it==NULL) continue;
		(*it)->convolution();
	}
}
 


void grainhdl::comparison_box(){
	for (int i = 1; i < grains.size(); i++){
		if(grains[i]==NULL)
			continue;
		grains[i]->comparison();
	}
}


void grainhdl::level_set(){
	for (int i = 1; i < grains.size(); i++){
		if(grains[i]==NULL)
			continue;
		grains[i]->find_contour();
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
	double total_energy= 0.0;
	int numberGrains=0;
	double totalLength=0;
	
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
// 		if(*it == NULL){ ??? better / faster
		if(*it!=NULL&&(*it)->get_status()){
			(*mymath).quaternion2Euler( (*it)->quaternion, euler );
			fprintf(myfile, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", euler[0], euler[1], euler[2], (*it)->volume, (float) (*it)->grainCharacteristics.size(), (float) (*it)->perimeter, (float) (*it)->energy);			
			total_energy += (*it)->energy;
			numberGrains+=1;
			for(it2=(*it)->grainCharacteristics.begin(); it2!=(*it)->grainCharacteristics.end(); it2++){
				if(!Settings::IsIsotropicNetwork){
				//take into account that every line is twice in the model
				discreteEnergyDistribution[(int)(((*it2).energyDensity)/dh -0.5) ] += 0.5 * (*it2).length;
				}
				totalLength += 0.5 * (*it2).length;
			}	
			
		}
	}
	if(!Settings::IsIsotropicNetwork){
		for (int i=0; i < Settings::DiscreteSamplingRate; i++){
				fprintf(enLenDis, "%lf\t%lf\n",(float)(dh*(i+1)),(float)discreteEnergyDistribution[i]);
				printf("%lf\t%lf\n",(float)(dh*(i+1)),(float)discreteEnergyDistribution[i]);
			}
	}
	totalenergy.push_back(0.5*total_energy);
	nr_grains.push_back(numberGrains);
	cout << "Timestep " << loop << " complete:" << endl;
	cout << "Number of grains remaining in the Network :" << nr_grains.back()<< endl;
	cout << "Amount of free Energy in the Network :" << totalenergy.back()<< endl;
	cout << "Total GB Length in Network :" << totalLength<< endl << endl << endl;
	fclose(myfile);
	fclose(enLenDis);
}
 
  
 
 
void grainhdl::run_sim(){
	find_neighbors();
// 	determineIDs();
	for(loop=0; loop <= Settings::NumberOfTimesteps; loop++){
		switchDistancebuffer();
		convolution();
		switchDistancebuffer();
		updateSecondOrderNeighbors();
		comparison_box();
		switchDistancebuffer();
		level_set();
		redistancing();
		if ( (loop % int(Settings::AnalysysTimestep)) == 0 || loop == Settings::NumberOfTimesteps ) {
			saveAllContourEnergies();
			save_texture();
			saveMicrostructure();
		}
		
	}
// 	utils::CreateMakeGif();
	cout << "Simulation complete." << endl;
}  

void grainhdl::saveMicrostructure(){
	stringstream param_xml_name;
	param_xml_name<< "PARAMETERS_SIM_NONAME_TIMESTEP_"<< loop <<"_GRAINS_"<<currentNrGrains<< ".xml";
	stringstream vertex_dump_name;
	vertex_dump_name<< "NETWORK_NONAME_TIMESTEP_"<< loop <<"_GRAINS_"<<currentNrGrains<< ".xml";

	createParamsForSim(param_xml_name.str().c_str(), vertex_dump_name.str().c_str());

}
void grainhdl::createParamsForSim(const char* param_filename, const char* vertex_dump_filename)
{
	xml_document<> doc_tree;

	xml_node<>* declaration = doc_tree.allocate_node(node_declaration);
	declaration->append_attribute(doc_tree.allocate_attribute("version", "1.0"));
	declaration->append_attribute(doc_tree.allocate_attribute("encoding", "utf-8"));
	doc_tree.append_node(declaration);

	doc_tree.append_node(Settings::generateXMLParametersNode(&doc_tree, vertex_dump_filename));
	ofstream output;
	output.open(param_filename);
	output<<doc_tree;
	output.close();

	cout<<doc_tree;
	int p =2;
}

 
void grainhdl::save_sim(){
// 	(*my_weights).plot_weightmap(ngridpoints, ID, ST, zeroBox);		
	ofstream myfile;
	myfile.open ("NrGrains&EnergyStatistics.txt");
	for(int i=1; i< nr_grains.size(); i++){
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
 
void grainhdl::clear_mem() {
	if (ST!=NULL) {delete  [] ST; }
}



