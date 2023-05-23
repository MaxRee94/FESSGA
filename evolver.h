#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <map>
#include "helpers.h"
#include <functional>


bool termination_condition_reached(uint* population, int pop_size, int no_cells, float threshold);

class Evolver {
public:
	Evolver() = default;
	Evolver(
		int _pop_size, float _mutation_rate, function<bool(uint*, int, int, float)> _terminate,
		int _tournament_size, uint* _base_densities, int _dim_x, int _dim_y, int _dim_z = 0,
		float _initial_perturbation_size = 0.5, float variance_treshold = 0.5
	) {
		base_densities = _base_densities;
		dim_x = _dim_x;
		dim_y = _dim_y;
		dim_z = _dim_z;
		if (dim_z == 0) {
			dim_z = 1; domain_2d = true;
		}
		no_cells = dim_x * dim_y * dim_z;
		pop_size = _pop_size;
		mutation_rate = _mutation_rate;
		tournament_size = _tournament_size;
		terminate = _terminate;
		initial_perturbation_size = _initial_perturbation_size;
	}
	void do_2d_crossover(uint* parent1, uint* parent2, uint* child1, uint* child2);
	void do_2d_mutation(uint* densities, float _mutation_rate);
	void init_population();
	void do_evolution();
private:
	int dim_x = 1;
	int dim_y = 1;
	int dim_z = 1;
	bool domain_2d = false;
	int pop_size = 1;
	float mutation_rate = 0;
	float initial_perturbation_size = 0;
	float variance_treshold = 0;
	int no_cells = 1;
	int tournament_size = 1;
	uint* base_densities = 0;
	uint* population = 0;
	function<bool(uint*, int, int, float)> terminate = 0;
};

