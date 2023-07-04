#include "individual.h"



bool fessga::evo::Individual2d::update_phenotype() {
	_copy(values, phenotype, _count, _phenotype_count);
	do_feasibility_filtering();
	remove_isolated_material();
	do_ground_element_filtering();
	return true;
}

bool fessga::evo::Individual2d::remove_isolated_material() {
	flush_edit_memory();
	init_pieces();
	if (pieces.size() == 1) return true;
	remove_smaller_pieces();
	return is_single_piece();
}

void fessga::evo::Individual2d::fill_voids(int _target_no_neighbors) {
	int target_no_neighbors = _target_no_neighbors;
	save_internal_snapshot();
	vector<int> x_bounds = { 0, dim_x - 1 };
	vector<int> y_bounds = { 0, dim_y - 1 };
	for (int x = 0; x < dim_x; x++) {
		for (int y = 0; y < dim_y; y++) {
			int cell = x * dim_y + y;

			// If the cell is already filled, continue to the next cell
			if (values[cell]) continue;

			// Determine the required number of neighbors for the cell to qualify for a fill operation.
			// This is dependent on whether the cell lies on the boundary of the design domain or not.
			target_no_neighbors = _target_no_neighbors;
			bool on_x_bound = help::is_in(&x_bounds, x);
			bool on_y_bound = help::is_in(&y_bounds, y);
			if (on_x_bound != on_y_bound) {
				target_no_neighbors -= 1;
			}
			else if (on_x_bound && on_y_bound) {
				// Only allow a minimum of 2 true neighbors
				if (_target_no_neighbors == 4) target_no_neighbors = 2;
				else continue;
			}

			// Fill the cell if the number of neighbors is equal to the required number computed earlier.
			vector<int> neighbors = get_true_neighbors(x, y, snapshot_internal);
			if (neighbors.size() == target_no_neighbors) {
				fill(cell);
			}
		}
	}
}

void fessga::evo::Individual2d::do_feasibility_filtering() {
	// Step 1: Fill all void elements that are surrounded on all sides by solid elements
	fill_voids(4);

	// Step 2: Remove all solid elements that have no true neighbors

	// Restore boundary cell elements

	// Step 3: Fill all void elements that are surrounded on all but one side by solid elements

	// Step 4: Remove all solid elements that have exactly one true neighbor

	// Restore boundary cell elements
}

void fessga::evo::Individual2d::do_ground_element_filtering() {
}

