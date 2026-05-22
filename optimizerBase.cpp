#pragma once
#include "optimizerBase.h"


bool load_physics(grd::Densities2d* densities, msh::SurfaceMesh* mesh, vector<int>* times, bool verbose) {
	// Obtain vtk file paths
	vector<string> vtk_paths;
	msh::get_vtk_paths(densities->fea_casemanager, densities->output_folder, vtk_paths);
	time_t start = time(0);
	string fea_finished_file = densities->output_folder + "/FEA_FINISHED.txt";
	for (auto& vtk_path : vtk_paths) {
		while (!IO::file_exists(fea_finished_file) || !IO::file_exists(vtk_path)) {
			float seconds_since_start = difftime(time(0), start);
			if (seconds_since_start > 10) {
				cout << "WARNING: physics loader has been waiting for >10 seconds for vtk file " << vtk_path << " to appear.\n";
				break;
			}
		}
		while (!IO::file_exists(vtk_path)) {
			if (difftime(time(0), start) > 60) {
				cout << "\nOptimizerBase: ERROR: Elmer did not produce a .vtk file after 60 seconds (expected path " << vtk_path << ")\n";
				return false;
			}
		}
	}
	densities->vtk_paths = vtk_paths;

	// Initialize data map to contain only 0's
	densities->fea_results.data_map.clear();
	for (int i = 0; i < densities->dim_x * densities->dim_y; i++) {
		if (densities->at(i)) densities->fea_results.data_map.insert(pair(i, 0));
	}

	// Load physics
	bool physics_loaded = fessga::phys::load_2d_physics_data(
		vtk_paths, densities->fea_results, densities->fea_casemanager, densities->dim_x, densities->dim_y,
		densities->cell_size, mesh->offset, densities->fea_casemanager->mechanical_constraint, &densities->border_nodes, verbose, times
	);

	// Check if loading was successful
	if (physics_loaded && verbose) cout << "OptimizerBase: Finished reading stress distribution from .vtk files." << endl;
	else if (!physics_loaded) {
		cout << "OptimizerBase: Error: Unable to read physics data." << endl;
	}


	return true;
}


