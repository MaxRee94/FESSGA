#include "evolver.h"
#include <iostream>


void Evolver::do_2d_crossover(uint* parent1, uint* parent2, uint* child1, uint* child2) {
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

/*
Mutate the given solution according to the set mutation rate
*/
void Evolver::do_2d_mutation(uint* densities, float _mutation_rate = -1) {
	for (int i = 0; i < domain_size; i++) {
		// Probability of a bit flip is equal to the mutation rate
		float rand_val = fessga::help::get_rand_float(0.0, 1.0);
		bool do_flip = rand_val < _mutation_rate;
		if (do_flip) {
			densities[i] = (int)(!densities[i]);
		}
	}
	// TODO: Test if having larger perturbations (flipping a group of multiple adjacent bits with low probability) is desirable
}

/*
Initialize a population of unique density distributions. Each differs slightly from the distribution loaded from file.
*/
void Evolver::init_population() {
	// Create uint buffer to store the binary density distributions of each individual, concatenated into a flat array
	population = new uint[pop_size * domain_size];

	// Fill the buffer with copies of the base densities
	for (int indiv = 0; indiv < pop_size; indiv++) {
		for (int cell = 0; cell < domain_size; cell++) {
			population[indiv * domain_size + cell] = base_densities[cell];
		}

		// Perturb each individual's density distribution through mutation
		do_2d_mutation(population + indiv * domain_size, initial_perturbation_size);

		// Run the genotype-phenotype mapping pipeline on each solution, to ensure feasibility.
	}
}

void Evolver::do_evolution() {

}

