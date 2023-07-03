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

void fessga::evo::Individual2d::do_feasibility_filtering() {
}

void fessga::evo::Individual2d::do_ground_element_filtering() {
}

