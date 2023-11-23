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
#include <thread>


class Evolver : public OptimizerBase {
public:
	Evolver() = default;
	Evolver(
		phys::FEACaseManager _fea_manager, msh::SurfaceMesh _mesh, string _base_folder, int _pop_size, int _no_static_iterations_trigger,
		float _mutation_rate_level0, float _mutation_rate_level1, grd::Densities2d _starting_densities, double _variation_trigger, int _max_iterations,
		int _max_iterations_without_change, bool _export_msh, bool _verbose, float _initial_perturb_level0, float _initial_perturb_level1,
		string _crossover_method, float _stress_fitness_influence
	) : OptimizerBase(
		_fea_manager, _mesh, _base_folder, _starting_densities, _max_iterations, _export_msh, _verbose)
	{
		pop_size = _pop_size;
		mutation_rate_level0 = _mutation_rate_level0;
		mutation_rate_level1 = _mutation_rate_level1;
		initial_perturb_level0 = _initial_perturb_level0;
		initial_perturb_level1 = _initial_perturb_level1;
		variation_trigger = _variation_trigger;
		no_static_iterations_trigger = _no_static_iterations_trigger;
		max_iterations_without_change = _max_iterations_without_change;
		crossover_method = _crossover_method;
		best_solutions_folder = output_folder + "/best_solutions";
		best_individuals_images_folder = image_folder + "/best_individuals";
		stress_fitness_influence = _stress_fitness_influence;
		IO::create_folder_if_not_exists(best_individuals_images_folder);
		IO::create_folder_if_not_exists(best_solutions_folder);
		img::write_distribution_to_image(densities, image_folder + "/starting_shape.jpg");
	}
	void start_FEA_threads(int pop_offset, vector<thread*> fea_threads);
	void do_2x_crossover(evo::Individual2d parent1, evo::Individual2d parent2, evo::Individual2d child1, evo::Individual2d child2);
	void do_ux_crossover(evo::Individual2d parent1, evo::Individual2d parent2, evo::Individual2d child1, evo::Individual2d child2);
	void do_2d_mutation(evo::Individual2d& densities, float _mutation_rate_level0, float _mutation_rate_level1);
	void create_valid_child_densities(vector<evo::Individual2d>* parents, vector<evo::Individual2d>& children);
	void init_population(bool verbose = true);
	void evolve();
	void do_setup();
	void create_children(bool verbose = true);
	void evaluate_fitnesses(int offset, bool do_FEA = false, bool verbose = true);
	void do_selection();
	void create_iteration_directories(int iteration);
	virtual void write_densities_to_image(bool verbose = false);
	bool termination_condition_reached();
	void choose_parents(vector<evo::Individual2d>& parents, vector<evo::Individual2d>* _population);
	void create_individual_mesh(evo::Individual2d* individual, bool verbose = false);
	void export_individual(evo::Individual2d* individual, string folder);
	void export_stats(string iteration_name, bool verbose = false);
	void collect_stats();
	void cleanup();
	void update_objective_function();
	void finish_FEA(int pop_offset, vector<thread*> fea_threads);
	void FEA_thread(vector<string> individual_folders, phys::FEACaseManager fea_casemanager, int pop_size, int thread_offset, bool verbose, int stepsize);
	void create_single_individual(bool verbose = false);
	virtual void export_meta_parameters(vector<string>* _ = 0) override;
	tuple<double, double, double, double, double, double, double> get_fitness_stats();
	vector<evo::Individual2d> population;
private:
	int pop_size = 1;
	float mutation_rate_level0, mutation_rate_level1, initial_perturb_level0, initial_perturb_level1;
	float initial_perturbation_size = 0;
	vector<float> fitness_time_series;
	float fitness_time_derivative = 0;
	bool last_iteration_was_valid = true;
	int final_valid_iteration = 1;
	int no_static_iterations_trigger = 100;
	string final_valid_iteration_folder = "";
	vector<string> individual_folders;
	map<int, double> fitnesses_map;
	PairSet fitnesses_pairset;
	vector<FILE*> pipes;
	double best_fitness = -INFINITY;
	double minimum_stress = INFINITY;
	float variation = 0;
	float stress_fitness_influence = 0;
	int iterations_since_fitness_change = 0;
	int max_iterations_without_change = 1;
	int no_unproductive_iterations = 0;
	double variation_trigger;
	string best_individual;
	int best_individual_idx = 0;
	string current_best_solution_folder;
	string best_solutions_folder;
	string best_individuals_images_folder;
	string crossover_method;
	double fitness_mean, fitness_stdev, relative_area_mean, relative_area_stdev, relative_max_stress_mean, relative_max_stress_stdev;
	vector<int> iterations_with_fea_failure;
};
