#pragma once
#include "optimizerBase.h"


void load_physics(grd::Densities2d* densities, msh::SurfaceMesh* mesh, bool verbose) {
	vector<string> vtk_paths;
	msh::get_vtk_paths(densities, vtk_paths);
	for (auto& vtk_path : vtk_paths) {
		if (!IO::file_exists(vtk_path)) {
			cout << "\nOptimizerBase: ERROR: Elmer did not produce a .vtk file (expected path " << vtk_path << ")\n";
			cout << "OptimizerBase: Terminating program." << endl;
			exit(1);
		}
	}
	bool physics_loaded = fessga::phys::load_2d_physics_data(
		vtk_paths, &densities->fea_results, densities->dim_x,
		densities->dim_y, densities->cell_size, mesh->offset, "Vonmises"
	);
	if (physics_loaded && verbose) cout << "OptimizerBase: Finished reading stress distribution from .vtk file." << endl;
	else if (!physics_loaded) {
		cout << "OptimizerBase: Error: Unable to read physics data." << endl;
	}
}


