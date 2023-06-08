#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <map>
#include "helpers.h"
#include "io.h"
#include "meshing.h"

using namespace fessga;


class OptimizerBase {
public:
	OptimizerBase() = default;
	OptimizerBase(
		string _msh_file, string _case_file, mesher::SurfaceMesh _mesh, string _output_folder,
		double _max_stress_threshold, uint* _densities, mesher::Grid3D _grid
	) {
		mesh = _mesh;
		msh_file = _msh_file; case_file = _case_file; output_folder = _output_folder;
		densities = _densities;
		max_stress_threshold = _max_stress_threshold;
		grid = _grid;
		domain_2d = true;
		no_cells = grid.x * grid.y * grid.z;
		IO::create_folder_if_not_exists(output_folder);

		// Check if output folder is empty
		//if (fessga::IO::is_empty(output_folder)) cout << "--- WARNING: The output folder " << output_folder << " is not empty. The existing files may be overwritten." << endl;
	}
	//void call_elmer(string folder);
	mesher::Grid3D grid;
	mesher::SurfaceMesh mesh;
	bool domain_2d = false;
	int no_cells = 1;
	double max_stress_threshold = 0.0;
	uint* densities = 0;
	string msh_file, case_file, output_folder;
};
