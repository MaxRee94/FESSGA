#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <map>
#include "helpers.h"


class OptimizerBase {
public:
	OptimizerBase() = default;
	OptimizerBase(uint* _base_densities, int _dim_x, int _dim_y, int _dim_z = 0) {
		base_densities = _base_densities;
		dim_x = _dim_x;
		dim_y = _dim_y;
		dim_z = _dim_z;
		if (dim_z == 0) {
			dim_z = 1; domain_2d = true;
		}
		no_cells = dim_x * dim_y * dim_z;
	}
private:
	int dim_x = 1;
	int dim_y = 1;
	int dim_z = 1;
	bool domain_2d = false;
	int no_cells = 1;
	uint* base_densities = 0;
};
