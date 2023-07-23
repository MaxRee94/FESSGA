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

