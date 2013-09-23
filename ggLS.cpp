#include "levelsetproject.h"

using namespace voro;

int main() {
	
	grainhdl my_sim;
	
	my_sim.setSimulationParameter();
	
	my_sim.run_sim();
	
	my_sim.save_sim();
	
	my_sim.clear_mem();	
	
}