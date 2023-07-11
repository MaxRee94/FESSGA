#include "evolver.h"
#include <iostream>

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
		to_string(iteration_number) + ", " + to_string(best_fitness) + ", " + to_string(variation) + ", " + best_individual
	);
}

void Evolver::do_2d_crossover(evo::Individual2d parent1, evo::Individual2d parent2, evo::Individual2d child1, evo::Individual2d child2) {
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

/*
Mutate the given solution according to the set mutation rate
*/
void Evolver::do_2d_mutation(evo::Individual2d& individual, float _mutation_rate = -1) {
	for (int i = 0; i < no_cells; i++) {
		// Probability of a bit flip is equal to the mutation rate
		float rand_val = fessga::help::get_rand_float(0.0, 1.0);
		bool do_flip = rand_val < _mutation_rate;
		if (do_flip) {
			individual.set(i, (int)(!individual[i]));
		}
	}
	// TODO: Test if having larger perturbations (flipping a group of multiple adjacent bits with low probability) is desirable
}

/*
Initialize a population of unique density distributions. Each differs slightly from the distribution loaded from file.
*/
void Evolver::init_population(bool verbose) {
	int i = 0;
	cout << "Generating initial population...\n";
	while (population.size() < pop_size) {
		i++;
		if (i > pop_size * 2 && population.size() == 0) {
			throw std::runtime_error("Error: Unable to generate any valid individuals after " + to_string(i) + " attempts.\n");
		}

		// Make a copy of the base individual
		evo::Individual2d individual(&densities);

		if (help::have_overlap(&individual.fea_case.cutout_cells, &individual.fea_case.cells_to_keep))
			cout << "Error: some cutout cells are also marked as keep cells.\n";

		/*cout << "keep cells:\n";
		individual.visualize_keep_cells();
		cout << "cutout cells:\n";
		individual.visualize_cutout_cells();*/

		// Perturb the individual's density distribution through mutation. This is done to add variation to the population.
		do_2d_mutation(individual, initial_perturbation_size);

		// Run the repair pipeline on each individual, to ensure feasibility.
		bool is_valid = individual.repair();
		if (!is_valid) continue; // If the repaired shape is not valid, re-try generating an individual

		// Export the individual's FEA mesh and case.sif file
		export_individual(&individual, individual_folders[population.size()]);
		
		// Add the individual to the population
		population.push_back(individual);

		if (verbose && (pop_size < 10 || population.size() % (pop_size / 10) == 0))
			cout << "- Generated individual " << population.size() << " / " << pop_size << "\n";
	}
	cout << "Generating initial population finished.\n";
}

void Evolver::write_densities_to_image() {
	string image_iteration_folder = image_folder + "/" + iteration_name;
	cout << "Exporting individuals to image files at directory location " << image_iteration_folder << endl;
	IO::create_folder_if_not_exists(image_iteration_folder);
	for (int i = 0; i < population.size(); i++) {
		img::write_distribution_to_image(
			population[i], image_iteration_folder + help::add_padding("/individual_", i+1) + to_string(i+1) + ".jpg");
	}
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
	cout << "Beginning Evolver run. Saving results to " << output_folder << endl;
	create_iteration_directories(iteration_number);
	image_folder = IO::create_folder_if_not_exists(output_folder + "/image_output");
	if (verbose) densities.print();
	img::write_distribution_to_image(densities, image_folder + "/starting_shape.jpg");
	init_population();
	evaluate_fitnesses(0);
	collect_stats();
	export_stats(iteration_name, true);
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
		do_2d_crossover(parents->at(0), parents->at(1), child1, child2);
		children = { child1, child2 };
		bool valid = true;
		for (auto& child : children) {
			do_2d_mutation(child, mutation_rate);
			valid = child.repair();
			if (!valid) break;
		}
		if (valid) break;
	}
}

void Evolver::create_individual_mesh(evo::Individual2d* individual, bool verbose) {
	msh::generate_FE_mesh(mesh, *individual, individual->fe_mesh);
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

void Evolver::generate_children(bool verbose) {
	cout << "Generating children...\n";
	vector<int> parent_indices;
	vector<evo::Individual2d> previous_population = population;
	for (int i = 0; i < (pop_size / 2); i++) {
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
}

void Evolver::evaluate_fitnesses(int offset, bool verbose) {
	cout << "Evaluating individual fitnesses...\n";
	iterations_since_fitness_change++;

	// Run FEA on all individuals in the population that have not yet been evaluated (usually only the newly generated children)
	for (int i = offset; i < (pop_size + offset); i++) {
		if (pipes.size() >= 6) {
			// If at least 6 pipes are running, wait for the first one to complete before starting a new one
			std::array<char, 80> buffer;
			while (fgets(buffer.data(), 80, pipes[0]) != NULL) {}
			_pclose(pipes[0]);
			pipes.erase(pipes.begin());
		}
		fessga::phys::call_elmer(population[i].output_folder + "/run_elmer.bat", &pipes, false);
		if (verbose && (population.size() < 20 || (i+1) % (pop_size / 5) == 0))
			cout << "- Finished FEA for individual " << i+1 - offset << " / " << pop_size << "\n";
	}
	for (auto& pipe : pipes) {
		std::array<char, 80> buffer;
		while (fgets(buffer.data(), 80, pipe) != NULL) {}
		_pclose(pipe);
	}

	// Obtain FEA results and compute fitnesses
#pragma omp parallel for
	for (int i = offset; i < (pop_size + offset); i++) {
		// Obtain FEA results
		load_physics(&population[i]);
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

		if (verbose && (population.size() < 20 || (i+1) % (pop_size / 5) == 0))
			cout << "- Evaluated individual " << i - offset + 1 << " / " << pop_size << "\n";
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
	best_individual = population[0].output_folder;
	vector<evo::Individual2d> new_population;
	map<int, double> new_fitnesses_map;
	for (auto& [pop_idx, fitness] : fitnesses_pairset) {
		new_fitnesses_map.insert(pair(new_population.size(), fitness));
		new_population.push_back(population[pop_idx]);
		if (new_population.size() == pop_size) break;
	}
	population = new_population;
	fitnesses_map = new_fitnesses_map;
}

void Evolver::evolve() {
	do_setup();
	while (!termination_condition_reached()) {
		iteration_number++;
		cout << "\nStarting iteration " << iteration_number << "...\n";
		create_iteration_directories(iteration_number);
		generate_children();
		evaluate_fitnesses(pop_size);
		do_selection();
		collect_stats();
		export_stats(iteration_name);
	}
}

