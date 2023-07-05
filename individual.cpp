#include "individual.h"


void fessga::evo::Individual2d::copy_from_individual(Individual2d* source) {
	for (int i = 0; i < size; i++) values[i] = source->at(i);
	_count = source->count();
}

void fessga::evo::Individual2d::update_phenotype() {
	_copy(values, phenotype, _count, _phenotype_count);
	do_ground_element_filtering();
}

bool fessga::evo::Individual2d::repair() {
	do_feasibility_filtering();
	remove_isolated_material();
	bool is_single_piece = pieces.size() == 1;
	bool is_valid = is_single_piece;
	return is_valid;
}

bool fessga::evo::Individual2d::remove_isolated_material() {
	flush_edit_memory();
	init_pieces();
	if (pieces.size() == 1) return true;
	remove_smaller_pieces();
	return is_single_piece();
}

void fessga::evo::Individual2d::fill_voids(int target_no_neighbors) {
	save_internal_snapshot();
	vector<int> x_bounds = { 0, dim_x - 1 };
	vector<int> y_bounds = { 0, dim_y - 1 };
	for (int x = 0; x < dim_x; x++) {
		for (int y = 0; y < dim_y; y++) {
			int cell = x * dim_y + y;

			// If the cell is already filled, continue to the next cell
			if (values[cell]) continue;

			// Fill the cell if the number of neighbors is equal to the required number computed earlier.
			vector<int> neighbors = get_true_neighbors(x, y, snapshot_internal);
			if (neighbors.size() == target_no_neighbors) {
				fill(cell);
			}
		}
	}
}

void fessga::evo::Individual2d::do_feasibility_filtering(bool verbose) {
	grd::Densities2d previous_state(this);
	bool filtering_had_effect = true;
	int i = 1;
	while (filtering_had_effect) {
		do_single_feasibility_filtering_pass();
		filtering_had_effect = !previous_state.is_identical_to(values);
		previous_state.copy_from(this);
		i++;
	}
	if (verbose) cout << "Performed feasibility filtering (" << i << " passes).\n";
}

void fessga::evo::Individual2d::do_single_feasibility_filtering_pass() {
	// Step 1: Fill void elements that are surrounded on all sides by solid elements
	fill_voids(4);

	// Step 2: Remove solid elements that have no true neighbors
	filter(0);

	// Restore boundary cell elements
	for (auto& cell : fea_case->boundary_cells) fill(cell);

	// Step 3: Fill void elements that are surrounded on all but one side by solid elements
	fill_voids(3);

	// Step 4: Remove solid elements that have exactly one true neighbor
	filter(1);

	// Restore boundary cell elements
	for (auto& cell : fea_case->boundary_cells) fill(cell);
}

void fessga::evo::Individual2d::do_ground_element_filtering() {
}

