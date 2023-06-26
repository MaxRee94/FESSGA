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
#include "meshing.h"
#include "individual.h"

using namespace fessga;


bool variation_minimum_passed(vector<evo::Individual2d> population, int pop_size, int no_cells, float threshold);


class Evolver : public OptimizerBase {
public:
	Evolver() = default;
	Evolver(
		string _msh_file, string _fe_case, msh::SurfaceMesh _mesh, string _output_folder, int _pop_size, float _mutation_rate,
		function<bool(vector<evo::Individual2d>, int, int, float)> _termination_condition, int _tournament_size, double _max_stress_threshold,
		evo::Individual2d _starting_densities, int _max_iterations, bool _export_msh = false, bool _verbose = true,
		float _initial_perturbation_size = 0.5, float variance_treshold = 0.5
	) : OptimizerBase(
		_msh_file, _fe_case, _mesh, _output_folder, _max_stress_threshold, _starting_densities, _max_iterations, _export_msh, _verbose)
	{
		pop_size = _pop_size;
		mutation_rate = _mutation_rate;
		tournament_size = _tournament_size;
		termination_condition = _termination_condition;
		initial_perturbation_size = _initial_perturbation_size;
	}
	void do_2d_crossover(evo::Individual2d parent1, evo::Individual2d parent2, evo::Individual2d child1, evo::Individual2d child2);
	void do_2d_mutation(evo::Individual2d densities, float _mutation_rate);
	void init_population();
	void do_evolution();
private:
	int pop_size = 1;
	float mutation_rate = 0;
	float initial_perturbation_size = 0;
	float variance_treshold = 0;
	int tournament_size = 1;
	vector<evo::Individual2d>* population = 0;
	function<bool(vector<evo::Individual2d>, int, int, float)> termination_condition = 0;
};
