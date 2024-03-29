#pragma once
#include "evolver.h"
#include <iostream>
#include <math.h>

//#define UNIFORM_POPULATION
//#define FEA_IGNORE
#define BEGIN_WITH_FILLED_DOMAIN

int NO_FEA_THREADS = 7;
int NO_RESULTS_THREADS = 2;

//vector<string> individual_folders, phys::FEACaseManager* fea_casemanager, 

/*
* Method to run a batch of FEA jobs on all given output folders
*/
void Evolver::FEA_thread(
	vector<string> individual_folders, phys::FEACaseManager fea_casemanager,
	int pop_size, int thread_offset, bool verbose, int stepsize
) {
	cout << "Starting FEA batch " + to_string(thread_offset + 1) + "\n";
	// Run FEA on all individuals in the population that have not yet been evaluated (usually only the newly generated children)
	int i = 0;
	for (int i = thread_offset; i < pop_size; i += NO_FEA_THREADS) {
		run_FEA_on_single_solution(individual_folders[i]);
		if (pop_size < 10 || ((i + 1) % (pop_size / 5) == 0))
			cout << "* FEA for individual " << i + 1 << " / " << pop_size << " finished\n";
	}
	IO::write_text_to_file(" ", individual_folders[0] + "/FEA_FINISHED.txt");
}

// Obtain FEA results
void load_physics_batch(
	double* max_stresses, vector<evo::Individual2d>* population, int pop_offset, int thread_offset, int pop_size,
	msh::SurfaceMesh mesh, bool verbose = true
) {
	if (verbose) cout << "Starting batch loader " << thread_offset << endl;
	vector<int> times;
	int j = -1;
	for (int i = pop_offset + thread_offset; i < (pop_offset + pop_size); i += NO_RESULTS_THREADS) {
		j++;
		// Load physics
		bool success = load_physics(&population->at(i), &mesh, &times);
		if (!success) { max_stresses[i] = INFINITY; continue; }
		else max_stresses[j] = population->at(i).fea_results.max;
		if (pop_size < 10 || (i + 1) % (pop_size / 5) == 0)
			cout << "- Read stress distribution for individual " << i - pop_offset + 1 << " / " << pop_size << "\n";
	}

	if (verbose && false) {
		cout << "\nSetup times:\n";
		for (int i = 0; i < population->at(i).fea_casemanager->active_cases.size() * 2 * ((pop_size) / NO_RESULTS_THREADS); i += 2) {
			cout << times[i] << endl;
		}
		cout << "\nRead times:\n";
		for (int i = 1; i < population->at(i).fea_casemanager->active_cases.size() * 2 * ((pop_size) / NO_RESULTS_THREADS); i += 2) {
			cout << times[i] << endl;
		}
	}
}

/*
Get variation within the given population.
Variation is not the same as variance; it is determined by the number of times each solution differs from another solution
where the 'other solution' is randomly chosen from the population (to avoid doing a costly n^2 comparison between all combinations
of solutions).
*/
float _get_variation(vector<evo::Individual2d>* population) {
	float sum_of_sq_diffs = 0;
	int cumulative_count = 0;
	for (int i = 0; i < population->size(); i++) {
		int individual_diff = 0;
		for (int c = 0; c < population->at(0).size; c++) {
			int other_indiv_idx = i;
			while (other_indiv_idx == i) {
				other_indiv_idx = fessga::help::get_rand_uint(0, population->size() - 1);
			}
			individual_diff += population->at(i)[c] != population->at(other_indiv_idx)[c];
		}
		sum_of_sq_diffs += individual_diff * individual_diff;
		cumulative_count += population->at(i).count();
	}
	cumulative_count = cumulative_count / population->size();
	float variation = sqrt(sum_of_sq_diffs) / (float)cumulative_count;

	return variation;
}

float get_variation(vector<evo::Individual2d>* population) {
	float variation = 0.0;
	int no_variation_runs = 5;
	for (int i = 0; i < no_variation_runs; i++) {
		variation += _get_variation(population);
	}
	return variation / (float)no_variation_runs;
}

// Get mean and standard deviation of fitnesses in current generation
tuple<double, double, double, double, double, double, double, double, double> Evolver::get_fitness_stats() { // Compute fitness mean
	vector<double> fitnesses;
	for (auto& [_, fitness] : fitnesses_map) fitnesses.push_back(fitness);
	double fit_mean = help::get_mean(&fitnesses);
	
	// Compute fitness stdev
	double fit_stdev = help::get_stdev(&fitnesses);

	// Compute the fitness derivative as a weighted average of inter-generation fitness differences, with
	// older differences receiving less weight.
	fitness_time_series.push_back(best_fitness);
	float weights = 0;
	float fitness_time_derivative = 0;
	int falloff_length = 50;
	float falloff_rate = 1.0 / (float)(falloff_length * falloff_length);
	fitness_time_derivative = 0;
	if (iteration_number > falloff_length) {
		for (int i = 1; i <= falloff_length; i++) {
			float weight = 1.0 - min((float)1.0, (falloff_rate * (float)i * (float)i));
			float fitness_difference = fitness_time_series[iteration_number - i] - fitness_time_series[iteration_number - i - 1];
			float weighted_difference = fitness_difference * weight;
			weights += weight;
			fitness_time_derivative += weighted_difference;
		}
		fitness_time_derivative /= weights;
	}
	else {
		// Set derivative to an arbitrary (positive) value if an insufficient number of fitness datapoints is available to compute a meaningful derivative approximation.
		fitness_time_derivative = 999;
	}

	// Compute mean relative area
	vector<double> relative_areas;
	for (auto& indiv : population) {
		relative_areas.push_back(indiv.get_relative_area());
	}
	double relative_area_mean = help::get_mean(&relative_areas);

	// Compute stdev relative area
	double relative_area_stdev = help::get_stdev(&relative_areas, relative_area_mean);

	// Compute mean maximum stress relative to threshold
	vector<double> relative_maximum_stresses;
	for (auto& indiv : population) {
		double relative_maximum_stress = indiv.fea_results.max / indiv.fea_casemanager->mechanical_threshold;
		relative_maximum_stresses.push_back(relative_maximum_stress);
	}
	double relative_maximum_stress_mean = help::get_mean(&relative_maximum_stresses);

	// Compute maximum stress relative to threshold stdev
	double relative_maximum_stress_stdev = help::get_stdev(&relative_maximum_stresses, relative_maximum_stress_mean);

	// Obtain displacement and maximum value for stress/yield criterion for the fittest individual in the population
	double fittest_yield_criterion = population[best_individual_idx].fea_results.max_yield_criterion;
	double fittest_displacement = population[best_individual_idx].fea_results.max_displacement;

	return { fit_mean, fit_stdev, fitness_time_derivative, relative_area_mean, relative_area_stdev, relative_maximum_stress_mean, relative_maximum_stress_stdev,
	fittest_yield_criterion, fittest_displacement };
}

/*
Termination condition based on the 'variation' within the population.
*/
bool variation_minimum_passed(vector<evo::Individual2d>* population, float threshold) {
	float variation = get_variation(population);
	if (variation < threshold) return true;
	else return false;
}

void Evolver::collect_stats() {
	variation = get_variation(&population);
	auto [
		_fitness_mean, _fitness_stdev, _fitness_time_derivative, _relative_area_mean,
			_relative_area_stdev, _relative_max_stress_mean, _relative_max_stress_stdev,
			_fittest_yield_criterion, _fittest_displacement
	] = get_fitness_stats();
	fitness_mean = _fitness_mean;
	fitness_stdev = _fitness_stdev;
	relative_area_mean = _relative_area_mean;
	relative_area_stdev = _relative_area_stdev;
	relative_max_stress_mean = _relative_max_stress_mean;
	relative_max_stress_stdev = _relative_max_stress_stdev;
	fitness_time_derivative = _fitness_time_derivative;
	fittest_yield_criterion = _fittest_yield_criterion;
	fittest_displacement = _fittest_displacement;
}

void Evolver::export_stats(string iteration_name, bool verbose) {
	string statistics_file = output_folder + "/statistics.csv";
	if (verbose) cout << "Exporting statistics to " << statistics_file << endl;
	if (initialize && !load_existing_population) {
		IO::write_text_to_file(
			"Iteration, Iteration time, Best fitness, Variation, Fitness mean, Fitness stdev, Fitness Derivative, Mean Relative Area, Stdev Relative Area, Mean Relative Max Stress, Stdev Relative Max Stress, Mutation rate (level 0), Mutation rate (level 1), Fittest Yield Criterion, Fittest Displacement, RAM available(GB), Virtual Memory available(GB), Pagefile available(GB), Percent memory used",
			statistics_file
		);
		return;
	}
	stats.push_back(to_string(best_fitness));
	stats.push_back(to_string(variation));
	stats.push_back(to_string(fitness_mean));
	stats.push_back(to_string(fitness_stdev));
	if (fitness_time_derivative == 999) stats.push_back("");
	else stats.push_back(to_string(fitness_time_derivative));
	stats.push_back(to_string(relative_area_mean));
	stats.push_back(to_string(relative_area_stdev));
	stats.push_back(to_string(relative_max_stress_mean));
	stats.push_back(to_string(relative_max_stress_stdev));
	stats.push_back(to_string(mutation_rate_level0));
	stats.push_back(to_string(mutation_rate_level1));
	stats.push_back(to_string(fittest_yield_criterion));
	stats.push_back(to_string(fittest_displacement));
	export_base_stats();
	vector<string> _stats = {
		"Current stats: \n   Variation = " + to_string(variation), "Best fitness = " + to_string(best_fitness),
		"Fitness derivative = " + to_string(fitness_time_derivative)
	};
	cout << help::join(&_stats, ", ") << endl;
}

// Do 2-point crossover
void Evolver::do_2x_crossover(
	evo::Individual2d parent1, evo::Individual2d parent2, evo::Individual2d child1, evo::Individual2d child2
) {
	vector<uint> crosspoints = { fessga::help::get_rand_uint(0, no_cells - 1), fessga::help::get_rand_uint(0, no_cells - 1) };
	uint crosspoint_1 = min(crosspoints[0], crosspoints[1]);
	uint crosspoint_2 = max(crosspoints[0], crosspoints[1]);
	//cout << "crosspoint 1: " << crosspoint_1 << ", crosspoint 2: " << crosspoint_2 << endl;
	for (int x = 0; x < densities.dim_x; x++) {
		for (int y = 0; y < densities.dim_y; y++) {
			int coord = x * densities.dim_y + y;
			if (coord > crosspoint_1 && coord < crosspoint_2) {
				child1.set(coord, parent1[coord]);
				child2.set(coord, parent2[coord]);
			}
			else {
				child1.set(coord, parent2[coord]);
				child2.set(coord, parent1[coord]);
			}
		}
	}
	child1.update_count();
	child2.update_count();
}

// Do uniform crossover
void Evolver::do_ux_crossover(
	evo::Individual2d parent1, evo::Individual2d parent2, evo::Individual2d child1, evo::Individual2d child2
) {
	for (int x = 0; x < densities.dim_x; x++) {
		for (int y = 0; y < densities.dim_y; y++) {
			int coord = x * densities.dim_y + y;
			int choice = help::get_rand_uint(0, 1);
			if (choice) {
				child1.set(coord, parent1[coord]);
				child2.set(coord, parent2[coord]);
			}
			else {
				child2.set(coord, parent1[coord]);
				child1.set(coord, parent2[coord]);
			}
		}
	}
	child1.update_count();
	child2.update_count();
}

/*
Mutate the given solution according to the set mutation rate
*/
void Evolver::do_2d_mutation(evo::Individual2d& individual, float _mutation_rate_level0, float _mutation_rate_level1) {
	// Level 0 mutation (bit-by-bit)
	for (int i = 0; i < no_cells; i++) {
		// Probability of a bit flip is equal to the mutation rate
		float rand_val = fessga::help::get_rand_float(0.0, 1.0);
		bool do_flip = rand_val < _mutation_rate_level0;
		if (do_flip) {
			individual.set(i, (int)(!individual[i]));
		}
	}

	// Level 1 mutation (groups of 4 adjacent bits arranged in a square)
	vector<pair<int, int>> offsets = { pair(0,0), pair(0,1), pair(1,1), pair(1,0) };
	for (int x = 0; x < individual.dim_x - 1; x++) {
		for (int y = 0; y < individual.dim_y - 1; y++) {
			float rand_val = fessga::help::get_rand_float(0.0, 1.0);
			bool do_flip = rand_val < _mutation_rate_level1;

			if (do_flip) {
				for (auto& offset : offsets) {
					int idx = individual.get_idx(x + offset.first, y + offset.second);
					individual.set(idx, (int)(!individual[idx]));
				}
			}
		}
	}

	individual.update_count();
}

void Evolver::create_single_individual(bool verbose) {
	// Make a copy of the base individual
	evo::Individual2d individual(&densities);

	/*cout << "keep cells:\n";
	individual.visualize_keep_cells();
	cout << "cutout cells:\n";
	individual.visualize_cutout_cells();*/

	float min_fraction_cells = 0.8;
	float max_fraction_cells = 1.2;

#ifdef BEGIN_WITH_FILLED_DOMAIN
	individual.fill_all();
	for (int& cutout_cell : fea_casemanager.cutout_cells) {
		individual.del(cutout_cell);
	}
	max_fraction_cells = INFINITY;
#endif

#ifndef UNIFORM_POPULATION:
	// Perturb the individual's density distribution through mutation. This is done to add variation to the population.
	do_2d_mutation(individual, initial_perturb_level0, initial_perturb_level1);
#endif

	// Run the repair pipeline on each individual, to ensure feasibility.
	bool is_valid = individual.repair();
	if (!is_valid) return; // If the repaired shape is not valid, abort (an attempt is then made to generate a replacement individual)

	
	// If the individual has more than the prescribed range of cells, discard it
	if (individual.count() > max_fraction_cells * densities.count()) {
		return;
	}
	
	// Iteratively fill the smallest fenestrae until the shape has the prescribed number of cells (randomly chosen within prescribed range)
	individual.fill_smaller_fenestrae((int)(help::get_rand_float(min_fraction_cells, max_fraction_cells) * (float)densities.count()), verbose);

	// Export the individual's FEA mesh and case.sif file
	export_individual(&individual, individual_folders[population.size()]);

	// Add the individual to the population
	population.push_back(individual);
}

void Evolver::load_individual(string densities_file) {
	// Create the individual and import its density values from a densities file
	evo::Individual2d individual(&densities);
	individual.do_import(densities_file, mesh.diagonal(0));

	// Export the individual's FEA mesh and case.sif file
	export_individual(&individual, individual_folders[population.size()]);

	// Add the individual to the population
	population.push_back(individual);
}

/*
Initialize a population of unique density distributions. Each differs slightly from the distribution loaded from file.
*/
void Evolver::generate_population(bool verbose) {
	cout << "Generating initial population...\n";

	// The first individual is an unperturbed copy of the given densities densities distribution
	evo::Individual2d individual(&densities);
	export_individual(&individual, individual_folders[population.size()]);
	population.push_back(individual);
	
	// Generate the first #no_threads individuals, and then start the FEA batch threads
	while (population.size() < NO_FEA_THREADS * 2) create_single_individual(verbose);
	
	// Run FEA
	thread fea_thread1, fea_thread2, fea_thread3, fea_thread4, fea_thread5, fea_thread6, fea_thread7, fea_thread8;
	vector<thread*> fea_threads = { &fea_thread1, &fea_thread2, &fea_thread3, &fea_thread4, &fea_thread5, &fea_thread6, &fea_thread7 };
	start_FEA_threads(0, fea_threads);

	int i = NO_FEA_THREADS * 2;
	while (population.size() < pop_size) {
		i++;
		if (i > pop_size * 2 && population.size() == 0) {
			throw std::runtime_error("Error: Unable to generate any valid individuals after " + to_string(i) + " attempts.\n");
		}
		create_single_individual();

		if (verbose && (pop_size < 10 || population.size() % (pop_size / 10) == 0))
			cout << "- Generated individual " << population.size() << " / " << pop_size << "\n";
	}
	cout << "Finished generating initial population.\n";
	finish_FEA(0, fea_threads);
}

void Evolver::load_population(bool verbose) {
	cout << "Loading population from folder...\n";

	int i = 0;
	while (population.size() < pop_size) {
		i++;

		string individual_densities_file = existing_population + "/individual_" + help::add_padding("", i) + to_string(i) + ".dens";
		load_individual(individual_densities_file);

		if (verbose && (pop_size < 10 || population.size() % (pop_size / 10) == 0))
			cout << "- Loaded individual " << population.size() << " / " << pop_size << "\n";
	}

	// Run FEA
	thread fea_thread1, fea_thread2, fea_thread3, fea_thread4, fea_thread5, fea_thread6, fea_thread7, fea_thread8;
	vector<thread*> fea_threads = { &fea_thread1, &fea_thread2, &fea_thread3, &fea_thread4, &fea_thread5, &fea_thread6, &fea_thread7 };
	start_FEA_threads(0, fea_threads);

	cout << "Finished loading population.\n";
	finish_FEA(0, fea_threads);
}

void Evolver::init_population(bool verbose) {
	if (load_existing_population) {
		load_population(verbose);
	}
	else {
		generate_population(verbose);
	}
}

void Evolver::start_FEA_threads(int pop_offset, vector<thread*> fea_threads) {
	for (int i = 0; i < NO_FEA_THREADS; i++) {
		*(fea_threads[i]) = thread(&Evolver::FEA_thread, this, individual_folders, fea_casemanager, pop_size, i, verbose, 0);
	}
}

void Evolver::finish_FEA(int pop_offset, vector<thread*> fea_threads) {
	// Create stress buffers for each FEA results loader
	cout << "Starting results loaders...\n";
	double* max_stresses1 = new double[pop_size / NO_RESULTS_THREADS];
	double* max_stresses2 = new double[pop_size / NO_RESULTS_THREADS];
	/*double* max_stresses3 = new double[pop_size / NO_RESULTS_THREADS];
	double* max_stresses4 = new double[pop_size / NO_RESULTS_THREADS];
	double* max_stresses5 = new double[pop_size / NO_RESULTS_THREADS];
	double* max_stresses6 = new double[pop_size / NO_RESULTS_THREADS];
	double* max_stresses7 = new double[pop_size / NO_RESULTS_THREADS];
	double* max_stresses8 = new double[pop_size / NO_RESULTS_THREADS];*/

	// Start results loaders
	thread results1_thread(load_physics_batch, max_stresses1, &population, pop_offset, 0, pop_size, mesh, true);
	/*thread results2_thread(load_physics_batch, max_stresses2, population, pop_offset, 1, pop_size, mesh, true);
	thread results3_thread(load_physics_batch, max_stresses3, population, pop_offset, 2, pop_size, mesh, true);
	thread results4_thread(load_physics_batch, max_stresses4, population, pop_offset, 3, pop_size, mesh, true);
	thread results5_thread(load_physics_batch, max_stresses5, population, pop_offset, 4, pop_size, mesh, true);
	thread results6_thread(load_physics_batch, max_stresses6, population, pop_offset, 5, pop_size, mesh, true);
	thread results7_thread(load_physics_batch, max_stresses7, population, pop_offset, 6, pop_size, mesh, true);*/
	
	// Also start a loader in the main thread
	load_physics_batch(max_stresses2, &population, pop_offset, NO_RESULTS_THREADS-1, pop_size, mesh, true);
	
	// Wait for all FEA threads to finish
	for (int i = 0; i < NO_FEA_THREADS; i++) {
		fea_threads[i]->join();
	}
	cout << "FEA of initial population finished.\n";

	// Wait for all FEA loaders to finish
	results1_thread.join();
	/*results2_thread.join();
	results3_thread.join();
	results4_thread.join();
	results5_thread.join();
	results6_thread.join();
	results7_thread.join();*/

	// Retrieve max stresses
	int j = 0;
	for (int i = 0; i < pop_size / NO_RESULTS_THREADS; i++) {
		population[pop_offset + i * NO_RESULTS_THREADS].fea_results.max = max_stresses1[j];
		population[pop_offset + i * NO_RESULTS_THREADS + 1].fea_results.max = max_stresses2[j];
		/*population[pop_offset + i * NO_RESULTS_THREADS + 2].fea_results.max = max_stresses3[j];
		population[pop_offset + i * NO_RESULTS_THREADS + 3].fea_results.max = max_stresses4[j];
		population[pop_offset + i * NO_RESULTS_THREADS + 4].fea_results.max = max_stresses5[j];
		population[pop_offset + i * NO_RESULTS_THREADS + 5].fea_results.max = max_stresses6[j];
		population[pop_offset + i * NO_RESULTS_THREADS + 6].fea_results.max = max_stresses7[j];
		population[pop_offset + i * NO_RESULTS_THREADS + 7].fea_results.max = max_stresses8[j];*/
		j++;
	}

	// Cleanup
	delete[] max_stresses1;
	delete[] max_stresses2;
	/*delete[] max_stresses3;
	delete[] max_stresses4;
	delete[] max_stresses5;
	delete[] max_stresses6;
	delete[] max_stresses7;
	delete[] max_stresses8;*/
}

void Evolver::write_densities_to_image(bool verbose) {
	string image_iteration_folder = image_folder + "/" + iteration_name;
	if (verbose) cout << "Exporting individuals to images at directory location " << image_iteration_folder << endl;
	IO::create_folder_if_not_exists(image_iteration_folder);
	for (int i = 0; i < population.size(); i++) {
		img::write_distribution_to_image(
			population[i], image_iteration_folder + help::add_padding("/individual_", i+1) + to_string(i+1) + ".jpg");
	}
	img::write_distribution_to_image(
		population[best_individual_idx], best_individuals_images_folder + "/" + iteration_name + ".jpg"
	);
}

void Evolver::create_iteration_directories(int iteration) {
	// Create iteration folder
	iteration_folder = get_iteration_folder(iteration, true, load_existing_population, existing_population);
	
	// Create individual folders
	individual_folders.clear();
	for (int i = 0; i < pop_size; i++) {
		string individual_folder = iteration_folder + help::add_padding("/individual_", i + 1) + to_string(i + 1);
		IO::create_folder_if_not_exists(individual_folder);
		individual_folders.push_back(individual_folder);
	}
}

bool Evolver::termination_condition_reached() {
	bool terminate = false;
	terminate = terminate || iteration_number >= max_iterations;
	if (iteration_number >= max_iterations) {
		terminate = true;
		cout << "\nTerminating emma: Maximum number of iterations (" + to_string(max_iterations) + ") reached.\n";
	}
	else if (iterations_since_fitness_change > max_iterations_without_change) {
		terminate = true;
		cout << "\nTerminating emma: Maximum number of iterations without a change in best fitness ("
			+ to_string(max_iterations_without_change) + ") reached.\n";
	}
	bool valid_solutions_exist = false;
	for (auto& [_, fitness] : fitnesses_map) if (fitness > -INFINITY) { valid_solutions_exist = true; break; }
	if (!valid_solutions_exist) {
		terminate = true;
		cout << "\nTerminating emma: All solutions in population are invalid.\n";
	}

	return terminate;
}

void Evolver::do_setup() {
	cout << "Beginning Evolver run. Saving results to " << output_folder << endl;
	if (load_existing_population) {
		// Retrieve iteration number from latest iteration folder path
		string latest_iteration_dir = IO::get_latest_path(output_folder + "/iteration#", "_", "_", 10000);
		cout << "latest dir: " << latest_iteration_dir << endl;
		vector<string> substrings;
		help::split(latest_iteration_dir, "iteration_", substrings);
		iteration_number = stoi(substrings[1]) + 1;
		cout << "iteration number: " << iteration_number << endl;
	}
	export_meta_parameters();
	create_iteration_directories(iteration_number);
	if (verbose) densities.print();
	
	time_t start = time(0);
	init_population();
	float seconds_since_start = difftime(time(0), start);
	cout << "Time taken to generate initial population: " << seconds_since_start << endl;
	for (auto& indiv : population) {
		for (auto& keepcell : indiv.fea_casemanager->keep_cells) {
			if (!indiv[keepcell]) {
				cout << "Invidiual with output folder " << indiv.output_folder << endl;
				indiv.print();
			}
		}
	}

	evaluate_fitnesses(0);
	help::sort(fitnesses_map, fitnesses_pairset);
	best_individual_idx = (*fitnesses_pairset.begin()).first;
	current_best_solution_folder = best_solutions_folder + "/" + iteration_name;
	IO::create_folder_if_not_exists(current_best_solution_folder);
	copy_solution_files(population[best_individual_idx].output_folder, current_best_solution_folder);
	collect_stats();
	if (!load_existing_population) export_stats(iteration_name);
	cleanup();
	initialize = false;
}

void Evolver::update_objective_function() {
	if (iterations_since_fitness_change >= no_static_iterations_trigger && variation < variation_trigger) {
		fea_casemanager.mechanical_threshold -= 1e5;
		cout << "-- Optimum shift triggered. Updated objective function. Maximum stress threshold changed from ("
			<< fea_casemanager.mechanical_threshold + 1e5 <<
			") to (" << fea_casemanager.mechanical_threshold << ").\n";
		cout << "-- Updating fitnesses according to new objective function.\n";
		fitnesses_map.clear();
		evaluate_fitnesses(0);
	}
}

// Choose a pair of parents randomly
void Evolver::choose_parents(vector<evo::Individual2d>& parents, vector<evo::Individual2d>* previous_population) {
	while (parents.size() < 2) {
		int parent_idx = help::get_rand_uint(0, previous_population->size() - 1);
		parents.push_back(previous_population->at(parent_idx));
		previous_population->erase(previous_population->begin() + parent_idx);
	}
}

// Perform crossover and mutation. If this yields invalid children, retry until valid children are obtained.
void Evolver::create_valid_child_densities(vector<evo::Individual2d>* parents, vector<evo::Individual2d>& children) {
	while (true) {
		evo::Individual2d child1(&densities), child2(&densities);
		if (crossover_method == "2x") do_2x_crossover(parents->at(0), parents->at(1), child1, child2);
		else if (crossover_method == "ux") do_ux_crossover(parents->at(0), parents->at(1), child1, child2);
		children = { child1, child2 };
		bool valid = true;
		for (auto& child : children) {
			do_2d_mutation(child, mutation_rate_level0, mutation_rate_level1);
			valid = child.repair();
			if (!valid) break;
		}
		if (valid) break;
	}
	int i = 0;
	for (auto& child : children) {
		child.output_folder = individual_folders[population.size() - pop_size + i];
		i++;
	}
}

void Evolver::create_individual_mesh(evo::Individual2d* individual, bool verbose) {
	msh::create_FE_mesh(mesh, *individual, individual->fe_mesh);
	msh::export_as_elmer_files(&individual->fe_mesh, individual->output_folder);
	if (export_msh) msh::export_as_msh_file(&individual->fe_mesh, individual->output_folder);
	if (verbose && IO::file_exists(individual->output_folder + "/mesh.header")) cout <<
		"emma: Exported new FE mesh to " << individual->output_folder << endl;
	else if (verbose) cout << "emma: ERROR: Failed to export new FE mesh.\n";
	string densities_file = individual->do_export(individual->output_folder + "/distribution2d.dens");
	string batch_file = msh::create_batch_file(individual->output_folder);
}

void Evolver::export_individual(evo::Individual2d* individual, string folder) {
	individual->output_folder = folder;
	individual->iteration = iteration_number;
	create_individual_mesh(individual);
	create_sif_files(individual, &individual->fe_mesh, verbose);
}

void Evolver::create_children(bool verbose) {
	cout << "Generating children...\n";
	vector<int> parent_indices;
	vector<evo::Individual2d> previous_population = population;

	// First create #NO_FEA_THREADS * 2 children to be able to begin FEA
	for (int i = 0; i < NO_FEA_THREADS; i++) {
		vector<evo::Individual2d> parents;
		choose_parents(parents, &previous_population);
		vector<evo::Individual2d> children;
		create_valid_child_densities(&parents, children);
		for (int j = 0; j < 2; j++) {
			export_individual(&children[j], individual_folders[i * 2 + j]);
			population.push_back(children[j]);
		}
		if (verbose && (population.size() < 20 || (i + 1) % (pop_size / 10) == 0))
			cout << "- Created child " << (i + 1) * 2 << " / " << pop_size << "\n";
	}
#ifndef FEA_IGNORE
	// Run FEA
	thread fea_thread1, fea_thread2, fea_thread3, fea_thread4, fea_thread5, fea_thread6, fea_thread7, fea_thread8, fea_thread9, fea_thread10;
	vector<thread*> fea_threads = { &fea_thread1, &fea_thread2, &fea_thread3, &fea_thread4, &fea_thread5, &fea_thread6, &fea_thread7, &fea_thread8, &fea_thread9, &fea_thread10 };
	start_FEA_threads(pop_size, fea_threads);
#endif
	// Generate rest of children
	for (int i = NO_FEA_THREADS; i < (pop_size / 2); i++) {
		vector<evo::Individual2d> parents;
		choose_parents(parents, &previous_population);
		vector<evo::Individual2d> children;
		create_valid_child_densities(&parents, children);
		for (int j = 0; j < 2; j++) {
			export_individual(&children[j], individual_folders[i * 2 + j]);
			population.push_back(children[j]);
		}
		if (verbose && (population.size() < 20 || (i+1) % (pop_size / 10) == 0))
			cout << "- Created child " << (i+1)*2 << " / " << pop_size << "\n";
	}
	cout << "Finished generating children.\n";
	finish_FEA(pop_size, fea_threads);
}

void Evolver::export_meta_parameters(vector<string>* _) {
	vector<string> additional_metaparameters = {
		"population size = " + to_string(pop_size),
		"initial perturbation size level 0 = " + to_string(initial_perturb_level0),
		"initial perturbation size level 1 = " + to_string(initial_perturb_level1),
		"initial perturbation size = " + to_string(initial_perturbation_size),
		"mutation rate level 0 = " + to_string(mutation_rate_level0),
		"mutation rate level 1 = " + to_string(mutation_rate_level1),
		"max iterations since fitness change = " + to_string(max_iterations_without_change),
		"crossover method = " + crossover_method,
		"mechanical constraint = " + fea_casemanager.mechanical_constraint,
		(help::is_in(fea_casemanager.mechanical_constraint, "Displacement") ? "max displacement = " + to_string(fea_casemanager.max_displacement) : "")
	};
	OptimizerBase::export_meta_parameters(&additional_metaparameters);
}

void Evolver::evaluate_fitnesses(int offset, bool do_FEA, bool verbose) {
	cout << "Evaluating individual fitnesses...\n";
	iterations_since_fitness_change++;

	// Obtain FEA results and compute fitnesses
	for (int i = offset; i < (pop_size + offset); i++) {
		if (verbose && (i % (pop_size / 5) == 0)) cout << "max stress: " << population[i].fea_results.max << endl;

		// Determine fitness
		double fitness;
		if (population[i].fitness == -INFINITY) {
			fitness = -INFINITY;
			if (!help::is_in(&iterations_with_fea_failure, iteration_number)) iterations_with_fea_failure.push_back(iteration_number);
		}
		else if (population[i].fea_results.max > fea_casemanager.mechanical_threshold) {
			fitness = 1.0 - population[i].fea_results.max / fea_casemanager.mechanical_threshold;
			if (help::is_in(fea_casemanager.mechanical_constraint, "Displacement") && fea_casemanager.mechanical_constraint != "Displacement") {
				double relative_displacement = population[i].fea_results.max_displacement / fea_casemanager.max_displacement;
				double relative_yield_criterion = population[i].fea_results.max_yield_criterion / fea_casemanager.mechanical_threshold;
				fitness = 2.0 - max(1.0, relative_displacement) - max(1.0, relative_yield_criterion);
			}
		}
		else {
			//double relative_maximum_stress = (fea_casemanager.mechanical_threshold - population[i].fea_results.max) / fea_casemanager.mechanical_threshold;
			//fitness = (relative_maximum_stress * stress_fitness_influence + 1.0) / (population[i].get_relative_area());
			fitness = 1.0 / population[i].get_relative_area();
		}
		if (verbose && (i % (pop_size/5) == 0)) cout << "fitness: " << fitness << endl;

		// Add fitness to map
		fitnesses_map.insert(pair(i, fitness));

		// Store stress value of strongest (i.e. lowest-stress/displacement) individual 
		if (population[i].fea_results.max < minimum_stress) {
			minimum_stress = population[i].fea_results.max;
		}

		// Update best fitness if improved
		if (fitness > best_fitness) {
			best_fitness = fitness;
			best_individual_idx = i;
			iterations_since_fitness_change = 0;
		}
	}

	if (iterations_since_fitness_change == 0) {
		cout << "EMMA: new best fitness: " << best_fitness << endl;
	}
	else {
		cout << "EMMA: Fitness (" << best_fitness << ") has remained unchanged for " << iterations_since_fitness_change << " iterations. " <<
			"emma will terminate when fitness has not changed for " << max_iterations_without_change << ".\n";
	}
}

void Evolver::do_selection() {
	cout << "Performing truncation selection...\n";

	// Obtain a list of population indices sorted according to fitness values
	help::sort(fitnesses_map, fitnesses_pairset);
	vector<pair<int, double>> individuals_and_fitnesses;
	for (auto& _pair : fitnesses_pairset) individuals_and_fitnesses.push_back(_pair);

	// Reverse the list so that it is ordered from large fitness values to small ones
	reverse(individuals_and_fitnesses.begin(), individuals_and_fitnesses.end());

	// Create a new fitness map (mapping pop indices to fitnesses) and a list of population indices corresponding to individuals
	// that should be removed from the population.
	map<int, double> new_fitnesses_map;
	vector<int> individuals_to_remove;

	// Store indices of individuals to be removed from population. Also store the fitness of individuals which will survive.
	for (auto& [pop_idx, fitness] : individuals_and_fitnesses) {
		if (new_fitnesses_map.size() < pop_size) {
			// Individual is selected to remain in the population
			new_fitnesses_map.insert(pair(pop_idx, fitness));
		}
		else {
			// Store the individual's index. Removal happens later.
			individuals_to_remove.push_back(pop_idx);
		}
	}
	fitnesses_map = new_fitnesses_map;

	// Sort removal vector in reverse order, so that individuals last in population will be removed first (done so that next individual's indices are not affected).
	sort(individuals_to_remove.begin(), individuals_to_remove.end(), greater<int>());
	
	// Erase individuals marked for removal from population
	for (auto& remove_idx : individuals_to_remove) {
		population[remove_idx].delete_arrays();
		population.erase(population.begin() + remove_idx);

		// Decrement all population indices larger than the removed index by one
		// The fitnesses_map is thereby made to track the shifts in the population.
		map<int, double> _fitnesses_map;
		for (auto& [keep_idx, fitness] : fitnesses_map) {
			if (keep_idx > remove_idx) {
				_fitnesses_map[keep_idx - 1] = fitness;
			}
			else _fitnesses_map[keep_idx] = fitness;
		}
		if (best_individual_idx > remove_idx) best_individual_idx--;
		fitnesses_map = _fitnesses_map;
	}

	// If current iteration produced a new best solution, export this solution to the 'best_solutions' folder
	if (iterations_since_fitness_change == 0) {
		population[best_individual_idx].print();

		string target_folder = IO::create_folder_if_not_exists(best_solutions_folder + "/" + iteration_name);
		copy_solution_files(population[best_individual_idx].output_folder, best_solutions_folder + "/" + iteration_name);
		current_best_solution_folder = best_solutions_folder + "/" + iteration_name;
		
		// Also write a superposition of stress values to the target folder as a .vtk file
		vector<string> vtk_paths;
		msh::get_vtk_paths(population[best_individual_idx].fea_casemanager, population[best_individual_idx].output_folder, vtk_paths);
		uint* densities = new uint[population[0].dim_x * population[0].dim_y];
		population[best_individual_idx].copy_to(densities);
		phys::write_results_superposition(
			vtk_paths, population[best_individual_idx].dim_x, population[best_individual_idx].dim_y,
			population[best_individual_idx].cell_size, mesh.offset, target_folder + "/SuperPosition.vtk", fea_casemanager.mechanical_constraint, &population[best_individual_idx].border_nodes,
			fea_casemanager.max_tensile_strength, fea_casemanager.max_compressive_strength
		);
		string superposition_copy_target = output_folder + "/superpositions/SuperPosition_iteration" + help::add_padding("_", iteration_number) + to_string(iteration_number) + ".vtk";
		cout << "Copying superposition to location " << superposition_copy_target << " ...\n";
		IO::copy_file(target_folder + "/SuperPosition.vtk", superposition_copy_target);
		delete[] densities;
	}

	// Print info on whether best solution exceeds MDT and/or MST
	cout << std::setprecision(4) << std::scientific;
	if (population[best_individual_idx].fea_results.max_displacement > fea_casemanager.max_displacement) {
		double relative_displacement = population[best_individual_idx].fea_results.max_displacement / fea_casemanager.max_displacement;
		cout << "Displacement triggered (" << population[best_individual_idx].fea_results.max_displacement << " > " << fea_casemanager.max_displacement << "). More severe than stress/yield criterion ("
			<< population[best_individual_idx].fea_results.max_yield_criterion << " / " << fea_casemanager.mechanical_threshold << ")? " <<
			(((relative_displacement * fea_casemanager.mechanical_threshold) > population[best_individual_idx].fea_results.max_yield_criterion) ? "Yes\n\n" : "No\n");
	}
	else {
		cout << "Displacement (" << population[best_individual_idx].fea_results.max_displacement << ") does not exceed MDT (" << fea_casemanager.max_displacement << ").\n";
	}
	cout << std::fixed;
}

void Evolver::cleanup() {
	if (iteration_number < 2) return;
	//if (help::is_in(&iterations_with_fea_failure, (iteration_number - 1))) return; // Skip removal of iterations with FEA failure
	string _iteration_folder = get_iteration_folder(iteration_number - 15);
	IO::remove_directory_incl_contents(_iteration_folder);
	cout << "Cleanup: Removed iteration directory and all contained files.\n";
}

void Evolver::do_local_search() {
	for (int i = 0; i < pop_size; i++) {
		population[i].do_thickening();
		export_individual(&population[i], individual_folders[i]);
	}
	// Run FEA
	thread fea_thread1, fea_thread2, fea_thread3, fea_thread4, fea_thread5, fea_thread6, fea_thread7, fea_thread8;
	vector<thread*> fea_threads = { &fea_thread1, &fea_thread2, &fea_thread3, &fea_thread4, &fea_thread5, &fea_thread6, &fea_thread7 };
	start_FEA_threads(0, fea_threads);
	finish_FEA(0, fea_threads);
}

void Evolver::write_population_to_disk() {
	// Remove current population
	IO::remove_directory_incl_contents(current_population_folder);
	IO::create_folder_if_not_exists(current_population_folder);

	// Write current population
	for (int i = 0; i < pop_size; i++) {
		string src = population[i].output_folder + "/distribution2d.dens";
		string trgt = current_population_folder + help::add_padding("/individual_", i + 1) + to_string(i + 1) + ".dens";
		IO::copy_file(src, trgt);
	}
}

void Evolver::do_iteration(bool _do_local_search) {
	iteration_number++;
	cout << "\nStarting iteration " << iteration_number << "...\n";
	if (fea_casemanager.dynamic) update_objective_function();
	create_iteration_directories(iteration_number);
	
	// Create new individuals and run FEA on them
	Timer timer; timer.start();
	if (_do_local_search) {
		do_local_search();
	}
	else {
		create_children();
	}
	timer.stop();
	cout << "Time taken to create children: " << timer.elapsedSeconds() << endl;

	// Evaluate fitnesses
	int pop_offset = pop_size;
	if (_do_local_search) pop_offset = 0;
	evaluate_fitnesses(pop_offset);

	// Double population size if variation is becoming too low.
	if (variation < 0.5 && (pop_size < 200) && !_do_local_search && iteration_number > 300) {
		cout << "INFO: INCREASING POPULATION SIZE FROM " << pop_size << " TO " << pop_size * 2 << endl;
		pop_size *= 2;
		collect_stats();
		export_stats(iteration_name);
		cleanup();
		return;
	}

	// Mutation boost
	if (iterations_since_fitness_change > 50 && mutation_boost_wait_time == 0) {
		// If fitness has not increased for half of maximum no iterations, boost the mutation rates to increase
		// population variance, in an attempt to push the algorithm to search outside the current local optimum. 
		mutation_rate_level0 *= mutation_boost_size;
		mutation_rate_level1 *= mutation_boost_size;
		cout << "MUTATION BOOST ON - Attempting to increase population variation.\n";
		cout << "	Increasing mutation rates by a factor of " << to_string(mutation_boost_size) << ".\n";
		mutation_boost = true;
		mutation_boost_wait_time = max_mutation_boost_wait_time; // Wait for a preset number of iterations before attempting another mutation boost.
	}
	else if (mutation_boost && variation > 1.0) {
		// If population variation has been restored, set the mutation rates back to their original values.
		mutation_rate_level0 /= mutation_boost_size;
		mutation_rate_level1 /= mutation_boost_size;
		cout << "MUTATION BOOST OFF - Restored mutation rate to original values.\n";
		mutation_boost = false;
	}
	if (mutation_boost_wait_time > 0) mutation_boost_wait_time--;

	write_population_to_disk();
	do_selection();
	collect_stats();
	export_stats(iteration_name);
	cleanup();
	//if (iterations_since_fitness_change == 0 && best_fitness > 0 && population[best_individual_idx].count() <= initial_count) {
	//	// If the surface area of the best individual is lower than the original input shape, shift the optimum 
	//	cout << std::setprecision(3) << std::scientific;
	//	cout << "SHIFTING OPTIMUM: Changing max stress threshold from " << fea_casemanager.mechanical_threshold << " to " << minimum_stress << endl;
	//	cout << std::fixed;
	//	fea_casemanager.mechanical_threshold = minimum_stress;
	//	export_meta_parameters();
	//}

	return;
	if (fitness_time_derivative < 0.01) {
		if (best_fitness > 0) return;
		if (no_unproductive_iterations < max_no_unproductive_iterations) {
			no_unproductive_iterations++;
			return;
		}
		float new_threshold = ((1.0 - best_fitness) * (float)fea_casemanager.mechanical_threshold) + ((float)fea_casemanager.mechanical_threshold / 1000.0);
		cout << "CHANGING MAX STRESS THRESHOLD: from " << fea_casemanager.mechanical_threshold << " to " << new_threshold << endl;
		fea_casemanager.mechanical_threshold = new_threshold;
		export_meta_parameters();
	}
	else {
		no_unproductive_iterations = 0;
	}
}

void Evolver::evolve() {
	do_setup();
	start_time = time(0);

	// Mutation boost parameters
	mutation_boost = false;
	max_mutation_boost_wait_time = max_iterations_without_change / 3;
	mutation_boost_wait_time = 0;
	mutation_boost_size = 2.0;

	// Other parameters
	max_no_unproductive_iterations = 1;

	// Evolution loop
	bool thickness_enhanced = false;
	while (!termination_condition_reached()) {
		do_iteration();
		continue;
		if (best_fitness < 0 && variation < 0.8 && !thickness_enhanced) {
			thickness_enhanced = true;
			// Do thickening on all individuals if the best fitness is still negative
			cout << "\n ----- Fitness is negative. Applying thickening to all individuals for the next iteration... -----\n";
			do_iteration(true);
		}
	}
}

