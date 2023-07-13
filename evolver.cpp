#include "evolver.h"
#include <iostream>
#include <thread>

int NO_FEA_THREADS = 4; // must be even number


/*
* Method to run a batch of FEA jobs on all given output folders
*/
void run_FEA_batch(vector<string> individual_folders, int pop_size, int thread_offset, bool verbose) {
	cout << "Starting FEA batch thread " + to_string(thread_offset) + "\n";
	// Run FEA on all individuals in the population that have not yet been evaluated (usually only the newly generated children)
	int i = 0;
	for (int i = thread_offset; i < pop_size; i += NO_FEA_THREADS) {
		string elmer_bat_file = individual_folders[i] + "/run_elmer.bat";
		while (!IO::file_exists(elmer_bat_file)) {} // Wait for elmer batfile to appear on disk
		fessga::phys::call_elmer(elmer_bat_file);
		if (verbose && (pop_size < 10 || (i + 1) % (pop_size / 5) == 0))
			cout << "- Finished FEA for individual " << i + 1 << " / " << pop_size << "\n";
	}
}

// Obtain FEA results
void load_physics_batch(vector<evo::Individual2d>* population, int pop_size, msh::SurfaceMesh* mesh, bool verbose = true) {
	cout << "Starting results loader...\n";
	for (int i = pop_size; i < pop_size * 2; i++) {
		string vtk_file_path = population->at(i).output_folder + "/case0001.vtk";
		while (!IO::file_exists(vtk_file_path)) {} // Wait for casefile to appear on disk
		load_physics(&population->at(i), mesh);
		if (verbose && (pop_size < 10 || (i + 1) % (pop_size / 5) == 0))
			cout << "- Read stress distribution for individual " << i - pop_size + 1 << " / " << pop_size << "\n";
	}
}

/*
Get variation within the given population.
Variation is not the same as variance; it is determined by the number of times each solution differs from another solution
where the 'other solution' is randomly chosen from the population (to avoid doing a costly n^2 check of all combinations
of solutions).
*/
float get_variation(vector<evo::Individual2d>* population) {
	float sum_of_sq_diffs = 0;
	int cumulative_count = 0;
	for (int i = 0; i < population->size(); i++) {
		int individual_diff = 0;
		int other_indiv_idx = i;
		while (other_indiv_idx == i) {
			other_indiv_idx = fessga::help::get_rand_uint(0, population->size() - 1);
		}
		for (int c = 0; c < population->at(0).size; c++) {
			individual_diff += population->at(i)[c] != population->at(other_indiv_idx)[c];
		}
		sum_of_sq_diffs += individual_diff * individual_diff;
		cumulative_count += population->at(i).count();
	}
	cumulative_count = cumulative_count / population->size();
	float variation = sqrt(sum_of_sq_diffs) / (float)cumulative_count;

	return variation;
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
}

void Evolver::export_stats(string iteration_name, bool initialize) {
	string statistics_file = output_folder + "/statistics.csv";
	cout << "Exporting statistics to " << statistics_file << endl;
	export_base_stats(iteration_name);
	if (initialize) IO::write_text_to_file(
		"Iteration, Best fitness, Variation, Best individual",
		statistics_file
	);
	IO::append_to_file(
		statistics_file,
		to_string(iteration_number) + ", " + to_string(best_fitness) + ", " + to_string(variation) + ", " + current_best_solution_folder
	);
	vector<string> stats = {
		"Current variation: " + to_string(variation),
		"current best individual (folder location): " + current_best_solution_folder
	};
	cout << help::join(&stats, ", ") << endl;
}

// Do 2-point crossover
void Evolver::do_2x_crossover(evo::Individual2d parent1, evo::Individual2d parent2, evo::Individual2d child1, evo::Individual2d child2) {
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
void Evolver::do_ux_crossover(evo::Individual2d parent1, evo::Individual2d parent2, evo::Individual2d child1, evo::Individual2d child2) {
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
	for (int x = 0; x < individual.dim_x; x+=2) {
		for (int y = 0; y < individual.dim_y; y += 2) {
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

	if (help::have_overlap(&individual.fea_case.cutout_cells, &individual.fea_case.cells_to_keep))
		cout << "Error: some cutout cells are also marked as keep cells.\n";

	/*cout << "keep cells:\n";
	individual.visualize_keep_cells();
	cout << "cutout cells:\n";
	individual.visualize_cutout_cells();*/

	// Perturb the individual's density distribution through mutation. This is done to add variation to the population.
	do_2d_mutation(individual, initial_perturb_level0, initial_perturb_level1);

	// Run the repair pipeline on each individual, to ensure feasibility.
	bool is_valid = individual.repair();
	if (!is_valid) return; // If the repaired shape is not valid, abort (an attempt is then made to generate a replacement individual, outside this function)

	// Export the individual's FEA mesh and case.sif file
	export_individual(&individual, individual_folders[population.size()]);

	// Add the individual to the population
	population.push_back(individual);
}

/*
Initialize a population of unique density distributions. Each differs slightly from the distribution loaded from file.
*/
void Evolver::init_population(bool verbose) {
	// Generate the first individual, and then start the FEA batch thread
	while (population.size() < NO_FEA_THREADS) create_single_individual(verbose);
	thread fea_thread1(run_FEA_batch, individual_folders, pop_size, 0, verbose);
	thread fea_thread2(run_FEA_batch, individual_folders, pop_size, 1, verbose);
	thread fea_thread3(run_FEA_batch, individual_folders, pop_size, 2, verbose);
	thread fea_thread4(run_FEA_batch, individual_folders, pop_size, 3, verbose);

	int i = NO_FEA_THREADS;
	cout << "Generating initial population...\n";
	while (population.size() < pop_size) {
		i++;
		if (i > pop_size * 2 && population.size() == 0) {
			throw std::runtime_error("Error: Unable to generate any valid individuals after " + to_string(i) + " attempts.\n");
		}

		create_single_individual();

		if (verbose && (pop_size < 10 || population.size() % (pop_size / 10) == 0))
			cout << "- Generated individual " << population.size() << " / " << pop_size << "\n";
	}
	cout << "Generating initial population finished.\n";
	fea_thread1.join();
	fea_thread2.join();
	fea_thread3.join();
	fea_thread4.join();
	cout << "FEA of initial population finished.\n";
}

void Evolver::write_densities_to_image() {
	string image_iteration_folder = image_folder + "/" + iteration_name;
	cout << "Exporting individuals to image files at directory location " << image_iteration_folder << endl;
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
	iteration_folder = get_iteration_folder(iteration, true);
	if (last_iteration_was_valid) {
		final_valid_iteration_folder = iteration_folder;
	}
	
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
		cout << "Terminating EVOMA: Maximum number of iterations (" + to_string(max_iterations) + ") reached.\n";
	}
	else if (iterations_since_fitness_change > max_iterations_without_change) {
		terminate = true;
		cout << "Terminating EVOMA: Maximum number of iterations without a change in best fitness (" + to_string(max_iterations_without_change) + ") reached.\n";
	}

	return terminate;
}

void Evolver::do_setup() {
	cout << "Beginning Evolver run. Saving results to " << base_folder << endl;
	export_meta_parameters();
	create_iteration_directories(iteration_number);
	if (verbose) densities.print();
	init_population();
	evaluate_fitnesses(0);
	help::sort(fitnesses_map, fitnesses_pairset);
	best_individual_idx = (*fitnesses_pairset.begin()).first;
	current_best_solution_folder = best_solutions_folder + "/" + iteration_name;
	IO::create_folder_if_not_exists(current_best_solution_folder);
	copy_solution_files(population[best_individual_idx].output_folder, current_best_solution_folder);
	collect_stats();
	export_stats(iteration_name, true);
	cleanup();
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
}

void Evolver::create_individual_mesh(evo::Individual2d* individual, bool verbose) {
	msh::create_FE_mesh(mesh, *individual, individual->fe_mesh);
	msh::export_as_elmer_files(&individual->fe_mesh, individual->output_folder);
	if (export_msh) msh::export_as_msh_file(&individual->fe_mesh, individual->output_folder);
	if (verbose && IO::file_exists(individual->output_folder + "/mesh.header")) cout << "EVOMA: Exported new FE mesh to " << individual->output_folder << endl;
	else if (verbose) cout << "EVOMA: ERROR: Failed to export new FE mesh.\n";
	string densities_file = individual->do_export(individual->output_folder + "/distribution2d.dens");
	string batch_file = msh::create_batch_file(individual->output_folder);
}

// Create and export a new version of the case.sif file by updating the boundary ids to fit the topology of the current FE mesh
void Evolver::create_sif_file(evo::Individual2d* individual, bool verbose) {
	map<string, vector<int>> bound_id_lookup;
	msh::create_bound_id_lookup(&bound_conds, &individual->fe_mesh, bound_id_lookup);
	msh::assemble_fea_case(&individual->fea_case, &bound_id_lookup);
	IO::write_text_to_file(individual->fea_case.content, individual->output_folder + "/case.sif");
	if (verbose) cout << "EVOMA: Exported updated case.sif file.\n";
}

void Evolver::export_individual(evo::Individual2d* individual, string folder) {
	individual->output_folder = folder;
	create_individual_mesh(individual);
	create_sif_file(individual);
}

void Evolver::create_children(bool verbose) {
	cout << "Generating children...\n";
	vector<int> parent_indices;
	vector<evo::Individual2d> previous_population = population;

	// First create #NO_FEA_THREADS children to be able to begin FEA
	for (int i = 0; i < (NO_FEA_THREADS / 2); i++) {
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
	thread fea_thread1(run_FEA_batch, individual_folders, pop_size, 0, verbose);
	thread fea_thread2(run_FEA_batch, individual_folders, pop_size, 1, verbose);
	thread fea_thread3(run_FEA_batch, individual_folders, pop_size, 2, verbose);
	thread fea_thread4(run_FEA_batch, individual_folders, pop_size, 3, verbose);

	// Generate rest of children
	for (int i = NO_FEA_THREADS / 2; i < (pop_size / 2); i++) {
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

	// Start a thread to read the results of the FEA
	thread results_thread(load_physics_batch, &population, pop_size, &mesh, verbose);

	// Join threads
	fea_thread1.join();
	fea_thread2.join();
	fea_thread3.join();
	fea_thread4.join();
	cout << "Finished FEA for all children.\n";
	results_thread.join();
	cout << "Finished reading FEA results for all children.\n";
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
		"crossover method = " + crossover_method
	};
	OptimizerBase::export_meta_parameters(&additional_metaparameters);
}

void Evolver::evaluate_fitnesses(int offset, bool do_FEA, bool verbose) {
	cout << "Evaluating individual fitnesses...\n";
	iterations_since_fitness_change++;

	// Obtain FEA results and compute fitnesses
	for (int i = offset; i < (pop_size + offset); i++) {
		max_stress = population[i].fea_results.max;

		// Compute fitness
		double fitness;
		if (max_stress > max_stress_threshold) fitness = INFINITY;
		else fitness = population[i].get_relative_volume();
		fitnesses_map.insert(pair(i, fitness));

		// Update best fitness if fitness is improved
		if (fitness < best_fitness) {
			best_fitness = fitness;
			iterations_since_fitness_change = 0;
		}
	}
	if (iterations_since_fitness_change == 0) {
		cout << "EVOMA: New best fitness: " << best_fitness << endl;
	}
	else {
		cout << "EVOMA: Fitness (" << best_fitness << ") has remained unchanged for " << iterations_since_fitness_change << " iterations. " <<
			"EVOMA will terminate when fitness has not changed for " << max_iterations_without_change << ".\n";
	}
}

void Evolver::do_selection() {
	cout << "Performing truncation selection...\n";
	help::sort(fitnesses_map, fitnesses_pairset);
	vector<evo::Individual2d> new_population;
	best_individual_idx = 0;
	map<int, double> new_fitnesses_map;
	for (auto& [pop_idx, fitness] : fitnesses_pairset) {
		new_fitnesses_map.insert(pair(new_population.size(), fitness));
		new_population.push_back(population[pop_idx]);
		if (new_population.size() == pop_size) break;
	}
	population = new_population;
	fitnesses_map = new_fitnesses_map;

	// If current iteration produced a new best solution, export this solution to the 'best_solutions' folder
	if (iterations_since_fitness_change == 0) {
		IO::create_folder_if_not_exists(best_solutions_folder + "/" + iteration_name);
		copy_solution_files(population[best_individual_idx].output_folder, best_solutions_folder + "/" + iteration_name);
		current_best_solution_folder = best_solutions_folder + "/" + iteration_name;
	}
}

void Evolver::cleanup() {
	IO::remove_directory_incl_contents(iteration_folder);
	cout << "Cleanup: Removed iteration directory and all contained files.\n";
}

void Evolver::evolve() {
	do_setup();
	while (!termination_condition_reached()) {
		iteration_number++;
		cout << "\nStarting iteration " << iteration_number << "...\n";
		create_iteration_directories(iteration_number);
		create_children();
		evaluate_fitnesses(pop_size);
		do_selection();
		collect_stats();
		export_stats(iteration_name);
		cleanup();
	}
}

