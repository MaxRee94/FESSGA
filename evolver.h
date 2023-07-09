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
#include "individual.h"


class Evolver : public OptimizerBase {
public:
	Evolver() = default;
	Evolver(
		string _msh_file, string _fea_case, msh::SurfaceMesh _mesh, string _output_folder, int _pop_size, float _mutation_rate,
		int _tournament_size, double _max_stress_threshold,
		grd::Densities2d _starting_densities, int _max_iterations, bool _export_msh = false, bool _verbose = true,
		float _initial_perturbation_size = 0.5
	) : OptimizerBase(
		_msh_file, _fea_case, _mesh, _output_folder, _max_stress_threshold, _starting_densities, _max_iterations, _export_msh, _verbose)
	{
		pop_size = _pop_size;
		mutation_rate = _mutation_rate;
		tournament_size = _tournament_size;
		initial_perturbation_size = _initial_perturbation_size;
	}
	void do_2d_crossover(evo::Individual2d parent1, evo::Individual2d parent2, evo::Individual2d child1, evo::Individual2d child2);
	void do_2d_mutation(evo::Individual2d& densities, float _mutation_rate);
	void create_valid_child_densities(vector<evo::Individual2d>* parents, vector<evo::Individual2d>& children);
	void init_population(bool verbose = true);
	void evolve();
	void do_setup();
	void generate_children();
	void evaluate_children(bool verbose = false);
	void do_selection();
	void create_iteration_folder_structure(int iteration);
	virtual void write_densities_to_image();
	bool termination_condition_reached();
	void choose_parents(vector<evo::Individual2d>& parents, vector<evo::Individual2d>* _population);
	void create_individual_mesh(evo::Individual2d* individual, bool verbose = false);
	void create_sif_file(evo::Individual2d* individual, bool verbose = false);
	vector<evo::Individual2d> population;
private:
	int pop_size = 1;
	float mutation_rate = 0;
	float initial_perturbation_size = 0;
	int tournament_size = 1;
	bool last_iteration_was_valid = true;
	int final_valid_iteration = 1;
	string final_valid_iteration_folder = "";
	vector<string> individual_folders;
	PairSet fitnesses;
};
