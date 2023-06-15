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

	// Constructor for 2d optimization
	OptimizerBase(
		string _msh_file, string _casefile, mesher::SurfaceMesh _mesh, string _output_folder,
		double _max_stress_threshold, uint* _densities, mesher::Grid3D _grid, int _max_iterations
	) {
		mesh = _mesh;
		msh_file = _msh_file; output_folder = IO::get_fullpath(_output_folder);
		casefile.path = _casefile;
		densities = _densities;
		max_stress_threshold = _max_stress_threshold;
		grid = _grid;
		domain_2d = true;
		no_cells = grid.x * grid.y * grid.z;
		max_iterations = _max_iterations;
		IO::create_folder_if_not_exists(output_folder);
		mesher::derive_boundary_conditions(densities, bound_conds, grid, mesh, casefile);
	};
	mesher::Grid3D grid;
	mesher::SurfaceMesh mesh;
	bool domain_2d = false;
	int no_cells = 1;
	double max_stress_threshold = 0.0;
	uint* densities = 0;
	int max_iterations = 0;
	map<string, vector<pair<int, int>>> bound_conds;
	physics::CaseFile casefile;
	string msh_file, output_folder;
};
