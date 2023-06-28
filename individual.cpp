#include "individual.h"



void fessga::evo::Individual2d::update_phenotype() {
	_copy(values, phenotype, _count, _phenotype_count);
	do_feasibility_filtering();
	remove_isolated_material();
	do_ground_element_filtering();
}

void fessga::evo::Individual2d::remove_isolated_material() {
	// Get number of pieces
	//vector<int> visited_cells;
	//densities.get_pieces(&visited_cells);

	// If multiple pieces, remove the ones that do not contain boundary cells
}

void fessga::evo::Individual2d::do_feasibility_filtering() {
}

void fessga::evo::Individual2d::do_ground_element_filtering() {
}

