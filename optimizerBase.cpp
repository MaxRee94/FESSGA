#pragma once
#include "optimizerBase.h"


bool load_physics(grd::Densities2d* densities, msh::SurfaceMesh* mesh, bool verbose) {
	// Obtain vtk file paths
	vector<string> vtk_paths;
	msh::get_vtk_paths(densities->fea_casemanager, densities->output_folder, vtk_paths);
	for (auto& vtk_path : vtk_paths) {
		if (!IO::file_exists(vtk_path)) {
			cout << "\nOptimizerBase: ERROR: Elmer did not produce a .vtk file (expected path " << vtk_path << ")\n";
			return false;
			//cout << "OptimizerBase: Terminating program." << endl;
			//exit(1);
			//cout << "OptimizerBase";
		}
	}
	densities->vtk_paths = vtk_paths;

	// Initialize data map to contain only 0's
	densities->fea_results.data_map.clear();
	for (auto& [coord, stress] : densities->fea_results.data_map) cout << "coord: " << coord << endl;
	for (int i = 0; i < densities->dim_x * densities->dim_y; i++) {
		if (densities->at(i)) densities->fea_results.data_map.insert(pair(i, 0));
	}

	// Load physics
	bool physics_loaded = fessga::phys::load_2d_physics_data(
		vtk_paths, densities->fea_results, densities->fea_casemanager, densities->dim_x, densities->dim_y,
		densities->cell_size, mesh->offset, densities->fea_casemanager->stress_type
	);

	// Check if loading was successful
	if (physics_loaded && verbose) cout << "OptimizerBase: Finished reading stress distribution from .vtk file." << endl;
	else if (!physics_loaded) {
		cout << "OptimizerBase: Error: Unable to read physics data." << endl;
	}


	return true;
}


