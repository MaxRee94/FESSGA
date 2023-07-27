#pragma once
#include "individual.h"


void fessga::evo::Individual2d::copy_from_individual(Individual2d* source) {
	for (int i = 0; i < size; i++) values[i] = source->at(i);
	_count = source->count();
}

void fessga::evo::Individual2d::update_phenotype() {
	_copy(values, phenotype, _count, _phenotype_count);
	do_ground_element_filtering();
}

void fessga::evo::Individual2d::do_ground_element_filtering() {
}

void fessga::evo::Individual2d::fill_smaller_fenestrae(int target_no_cells, bool verbose) {

	// Obtain vector of all fenestrae
	// -- Create densities2d object with inverted values of current individual
	grd::Densities2d inverted(dim_x, diagonal, output_folder);
	inverted.copy_from(this);
	inverted.invert();

	// -- Run 'init_pieces' on inverted densities
	inverted.init_pieces(fea_casemanager->cutout_cells[0]);

	// -- Obtain vector of inverted pieces (the fenestrae + the space surrounding the skull)
	vector<grd::Piece> fenestrae = inverted.pieces;

	// -- Detect which piece represents the piece surrounding the skull by searching for edge cells
	for (int i = 0; i < fenestrae.size(); i++) {
		bool is_surrounding_space = false;
		for (int j = 0; j < fenestrae[i].cells.size(); j++) {
			int cell = fenestrae[i].cells[j];
			int x = cell / dim_y;
			int y = cell % dim_y;
			is_surrounding_space = x == 0 || y == 0 || x == (dim_x - 1) || y == (dim_y - 1);
		}

		// -- Remove the found piece from the obtained vector. It should then contain only fenestrae.
		if (is_surrounding_space) {
			fenestrae.erase(fenestrae.begin() + i);
			break;
		}
	}
	
	// Create map from vector indices to number of cells per fenestra
	map<int, double> fenestrae_map;
	for (int i = 0; i < fenestrae.size(); i++) {
		fenestrae_map.insert(pair(i, (double)fenestrae[i].cells.size()));
	}

	// Get PairSet of fenestrae sorted according to number of cells per fenestra
	PairSet fenestrae_set;
	help::sort(fenestrae_map, fenestrae_set);

	// Iteratively fill fenestrae until count() >= target_no_cells
	for (auto& [idx, no_cells] : fenestrae_set) {

		// Don't fill fenestra that contain cutout cells
		bool contains_cutout_cells = false;
		for (auto& cell : fenestrae[idx].cells) {
			if (help::is_in(&fea_casemanager->cutout_cells, cell)) { contains_cutout_cells = true; break; }
		}
		if (contains_cutout_cells) break;

		// Fill all cells in the fenestra
		bool stop = false;
		for (auto& cell : fenestrae[idx].cells) {
			fill(cell);
			if (count() >= target_no_cells) {
				stop = true; break;
			}
		}
		if (stop) break;
	}
}

