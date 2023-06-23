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
		string _msh_file, string _fe_case, mesher::SurfaceMesh _mesh, string _output_folder,
		double _max_stress_threshold, uint* _densities, mesher::Grid3D _grid, int _max_iterations, bool _export_msh, bool _verbose
	) {
		mesh = _mesh;
		msh_file = _msh_file; output_folder = IO::get_fullpath(_output_folder);
		fe_case.path = _fe_case;
		densities = _densities;
		export_msh = _export_msh;
		max_stress_threshold = _max_stress_threshold;
		grid = _grid;
		domain_2d = true;
		no_cells = grid.x * grid.y;
		max_iterations = _max_iterations;
		verbose = _verbose;
		IO::create_folder_if_not_exists(output_folder);
		mesher::derive_boundary_conditions(densities, bound_conds, grid, mesh, fe_case);

		//for (auto& [bound_name, pairs] : bound_conds) {
		//	cout << "bound name: " << bound_name << endl;
		//	cout << "pairs: " << help::join_as_string(pairs, " ") << endl;
		//}
	};
	mesher::Grid3D grid;
	mesher::SurfaceMesh mesh;
	bool domain_2d = false;
	int no_cells = 1;
	double max_stress_threshold = 0.0;
	uint* densities = 0;
	bool verbose = true;
	int max_iterations = 0;
	bool export_msh = false;
	map<string, vector<pair<int, int>>> bound_conds;
	mesher::Case fe_case;
	string msh_file, output_folder;

	// Function to get the folder corresponding to the given iteration number. If the folder does not exist yet, it will be created.
	string get_iteration_folder(int iteration, string& cur_iteration_name, bool verbose = false) {
		cur_iteration_name = fessga::help::add_padding("iteration_", iteration) + to_string(iteration);
		string cur_output_folder = output_folder + "/" + cur_iteration_name;
		if (verbose) cout << "FESS: Creating folder " << cur_output_folder << " for current iteration if it does not exist yet...\n";
		IO::create_folder_if_not_exists(cur_output_folder);

		return cur_output_folder;
	}
};


