#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <map>
#include <functional>
#include "helpers.h"
#include "optimizerBase.h"


bool variation_minimum_passed(uint* population, int pop_size, int no_cells, float threshold);


class Evolver : public OptimizerBase {
public:
	Evolver() = default;
	Evolver(
		string _msh_file, string _casefile, mesher::SurfaceMesh _mesh, string _output_folder, int _pop_size, float _mutation_rate,
		function<bool(uint*, int, int, float)> _termination_condition, int _tournament_size, double _max_stress_threshold,
		uint* _starting_densities, mesher::Grid3D _grid, int _max_iterations, float _initial_perturbation_size = 0.5, float variance_treshold = 0.5
	) : OptimizerBase(_msh_file, _casefile, _mesh, _output_folder, _max_stress_threshold, _starting_densities, _grid, _max_iterations) {
		pop_size = _pop_size;
		mutation_rate = _mutation_rate;
		tournament_size = _tournament_size;
		termination_condition = _termination_condition;
		initial_perturbation_size = _initial_perturbation_size;
	}
	void do_2d_crossover(uint* parent1, uint* parent2, uint* child1, uint* child2);
	void do_2d_mutation(uint* densities, float _mutation_rate);
	void init_population();
	void do_evolution();
private:
	int pop_size = 1;
	float mutation_rate = 0;
	float initial_perturbation_size = 0;
	float variance_treshold = 0;
	int tournament_size = 1;
	uint* population = 0;
	function<bool(uint*, int, int, float)> termination_condition = 0;
};
