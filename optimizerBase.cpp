#include "optimizerBase.h"


void load_physics(grd::Densities2d* densities, msh::SurfaceMesh* mesh, bool verbose) {
	string cur_case_output_file = densities->output_folder + "/case0001.vtk";
	if (!IO::file_exists(cur_case_output_file)) {
		cout << "\nOptimizerBase: ERROR: Elmer did not produce a .vtk file (expected path " << cur_case_output_file << ")\n";
		cout << "OptimizerBase: Terminating program." << endl;
		exit(1);
	}
	bool physics_loaded = fessga::phys::load_2d_physics_data(
		cur_case_output_file, &densities->fea_results, densities->dim_x, densities->dim_y, densities->cell_size, mesh->offset, "Vonmises"
	);
	if (physics_loaded && verbose) cout << "OptimizerBase: Finished reading stress distribution from .vtk file." << endl;
	else if (!physics_loaded) {
		cout << "OptimizerBase: Error: Unable to read physics data from file " << cur_case_output_file << endl;
	}
}


