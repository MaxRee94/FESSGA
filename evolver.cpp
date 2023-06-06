#include "evolver.h"
#include <iostream>


/*
Termination condition based on the 'variation' within the population.
Variation is not the same as variance; it is determined by the number of times each solution differs from another solution
where the 'other solution' is randomly chosen from the population (to avoid doing a costly n^2 check of all combinations
of solutions).
*/
bool variation_minimum_passed(uint* population, int pop_size, int no_cells, float threshold) {
	float sum_of_sq_diffs = 0;
	for (int i = 0; i < pop_size; i++) {
		int individual_diff = 0;
		int other_indiv_idx = i;
		while (other_indiv_idx == i) {
			other_indiv_idx = fessga::help::get_rand_uint(0, pop_size - 1);
		}
		for (int c = 0; c < no_cells; c++) {
			individual_diff += (int)population[i * no_cells + c] - (int)population[other_indiv_idx * no_cells + c];
		}
		sum_of_sq_diffs += individual_diff * individual_diff;
	}
	float variation = sqrt(sum_of_sq_diffs) / (float)(pop_size * no_cells);
	if (variation < threshold) return true;
	else return false;
}


void Evolver::do_2d_crossover(uint* parent1, uint* parent2, uint* child1, uint* child2) {
	vector<uint> crosspoints = { fessga::help::get_rand_uint(0, no_cells - 1), fessga::help::get_rand_uint(0, no_cells - 1) };
	uint crosspoint_1 = min(crosspoints[0], crosspoints[1]);
	cout << "cross 1: " << crosspoint_1 << endl;
	uint crosspoint_2 = max(crosspoints[0], crosspoints[1]);
	cout << "cross 2: " << crosspoint_2 << endl;
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
	for (int i = 0; i < no_cells; i++) {
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
	// Create uint buffer to store the binary density distributions of each individual, concatenated into a single array
	population = new uint[pop_size * no_cells];

	// Fill the buffer with (mutated) copies of the base densities
	for (int indiv = 0; indiv < pop_size; indiv++) {
		for (int cell = 0; cell < no_cells; cell++) {
			population[indiv * no_cells + cell] = starting_densities[cell];
		}

		// Perturb each individual's density distribution through mutation. This is done to add variation to the population.
		do_2d_mutation(population + indiv * no_cells, initial_perturbation_size);

		// Run the genotype-phenotype mapping pipeline on each solution, to ensure feasibility.
	}
}

void Evolver::do_evolution() {

}

