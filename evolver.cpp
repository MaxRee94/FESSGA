#include "evolver.h"
#include <iostream>


void Evolver::do_2d_crossover(uint* parent1, uint* parent2, uint* child1, uint* child2, int dim_x, int dim_y) {
	int domain_size = dim_x * dim_y;
	vector<uint> crosspoints = { fessga::help::get_rand_uint(0, domain_size), fessga::help::get_rand_uint(0, domain_size) };
	uint crosspoint_1 = min(crosspoints[0], crosspoints[1]);
	uint crosspoint_2 = max(crosspoints[0], crosspoints[1]);
	//cout << "crosspoint 1: " << crosspoint_1 << ", crosspoint 2: " << crosspoint_2 << endl;
	for (int x = 0; x < dim_x; x++) {
		for (int y = 0; y < dim_y; y++) {
			int coord = x * dim_x + y;
			if (coord > crosspoint_1 && coord < crosspoint_2) {
				child1[coord] = parent1[coord];
				child2[coord] = parent2[coord];
			}
			else {
				child1[coord] = parent2[coord];
				child2[coord] = parent1[coord];
			}
		}
	}
}

void Evolver::init_population() {
	population = new uint[pop_size * dim_x * dim_y * dim_z];
}

void Evolver::do_evolution() {

}

