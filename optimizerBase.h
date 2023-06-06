#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <map>
#include "helpers.h"
#include "io.h"


class OptimizerBase {
public:
	OptimizerBase() = default;
	OptimizerBase(string _msh_file, string _case_file, string _output_folder, double _max_stress_threshold, uint* _starting_densities, int _dim_x, int _dim_y, int _dim_z = 0) {
		msh_file = _msh_file; case_file = _case_file; output_folder = _output_folder;
		starting_densities = _starting_densities;
		max_stress_threshold = _max_stress_threshold;
		dim_x = _dim_x;
		dim_y = _dim_y;
		dim_z = _dim_z;
		if (dim_z == 0) {
			dim_z = 1; domain_2d = true;
		}
		no_cells = dim_x * dim_y * dim_z;
		fessga::IO::create_folder_if_not_exists(output_folder);

		// Check if output folder is empty
		//if (fessga::IO::is_empty(output_folder)) cout << "--- WARNING: The output folder " << output_folder << " is not empty. The existing files may be overwritten." << endl;
	}
	//void call_elmer(string folder);
	int dim_x = 1;
	int dim_y = 1;
	int dim_z = 1;
	bool domain_2d = false;
	int no_cells = 1;
	double max_stress_threshold = 0.0;
	uint* starting_densities = 0;
	string msh_file, case_file, output_folder;
};


//void OptimizerBase::call_elmer(string folder) {
//
//}

