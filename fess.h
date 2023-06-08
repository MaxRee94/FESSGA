#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <map>
#include <functional>
#include "helpers.h"
#include "optimizerBase.h"


class FESS : public OptimizerBase {
public:
	FESS() = default;
	FESS(
		string _msh_file, string _case_file, mesher::SurfaceMesh _mesh, string _output_folder, double _min_stress_threshold,
		double _max_stress_threshold, uint* _starting_densities, mesher::Grid3D _grid
	) : OptimizerBase(_msh_file, _case_file, _mesh, _output_folder, _max_stress_threshold, _starting_densities, _grid)
	{
		min_stress_threshold = _min_stress_threshold;
	}
	double min_stress_threshold = 1.0;
	void run();
};


void FESS::run() {
	cout << "Beginning FESS run. Saving results to " << output_folder << endl;

	double min_stress = 1.0;
	string msh = msh_file;
	string cur_iteration_name = "iteration_0001";
	string cur_output_folder = output_folder + "/" + cur_iteration_name;
	int i = 0;
	while (min_stress > min_stress_threshold) {
		cout << "--- Starting iteration " << i << ". Previous lowest stress value: " << min_stress << " N/m^2" << endl;

		// Call Elmer to run FEA on .msh file in cur_output_folder, using .sif file
		physics::call_elmer(output_folder + "/" + cur_iteration_name);

		// Wait for Elmer's analysis to complete. This is the case when a new .vtk file has appeared

		// Obtain stress distribution from the .vtk file
		min_stress = 0.0;

		// Set all cells that have corresponding stress values lower than <min_stress_threshold> to density=0.

		// Create new subfolder for output of current iteration
		string cur_iteration_name = fessga::help::add_padding("iteration_", i + 1) + to_string(i + 1);
		string cur_output_folder = output_folder + "/" + cur_iteration_name;
		cout << "Created output folder " << cur_output_folder << " for current iteration.\n";
		IO::create_folder_if_not_exists(cur_output_folder);

		// Copy the case.sif file to the newly created subfolder
		IO::copy_file(case_file, cur_output_folder + "/case.sif");
		if (IO::FileExists(cur_output_folder + "/case.sif")) cout << "Copied case file to subfolder.\n";
		else cout << "ERROR: Failed to copy case file to subfolder.\n";
		
		// Generate new FE mesh using modified density distribution
		cout << "Generating new FE mesh...\n";
		mesher::FEMesh2D fe_mesh;
		mesher::generate_FE_mesh(grid, mesh, densities, fe_mesh);
		cout << "FE mesh generation done.\n";

		// Export newly generated FE mesh
		mesher::export_as_elmer_files(&fe_mesh, cur_output_folder);
		if (IO::FileExists(cur_output_folder + "/mesh.header")) cout << "Exported new FE mesh to subfolder.\n";
		else cout << "ERROR: Failed to export new FE mesh to subfolder.\n";

		cout << endl;

		i++;
	}

	cout << "\nFinished FESS run. Results were saved to " << output_folder << endl;
}
