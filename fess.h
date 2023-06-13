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
		double _max_stress_threshold, uint* _starting_densities, mesher::Grid3D _grid, int _max_iterations
	) : OptimizerBase(_msh_file, _case_file, _mesh, _output_folder, _max_stress_threshold, _starting_densities, _grid, _max_iterations)
	{
		min_stress_threshold = _min_stress_threshold;
	}
	double min_stress_threshold = 1.0;
	void run();
};


void FESS::run() {
	cout << "Beginning FESS run. Saving results to " << output_folder << endl;

	double min_stress, max_stress;
	string msh = msh_file;
	string cur_iteration_name = "";
	string cur_output_folder = output_folder;
	int i = 1;
	while (i - 1 < max_iterations) {
		cout << "\nFESS: Starting iteration " << i << ".\n";

		// Create new subfolder for output of current iteration
		string cur_iteration_name = fessga::help::add_padding("iteration_", i) + to_string(i);
		string cur_output_folder = output_folder + "/" + cur_iteration_name;
		cout << "FESS: Created output folder " << cur_output_folder << " for current iteration.\n";
		IO::create_folder_if_not_exists(cur_output_folder);

		// Copy the case.sif file to the newly created subfolder
		IO::copy_file(case_file, cur_output_folder + "/case.sif");
		if (IO::file_exists(cur_output_folder + "/case.sif")) cout << "FESS: Copied case file to subfolder.\n";
		else cout << "FESS: ERROR: Failed to copy case file to subfolder.\n";
		
		// Generate new FE mesh using modified density distribution
		cout << "FESS: Generating new FE mesh...\n";
		mesher::FEMesh2D fe_mesh;
		mesher::generate_FE_mesh(grid, mesh, densities, fe_mesh);
		cout << "FESS: FE mesh generation done.\n";

		// Export newly generated FE mesh
		mesher::export_as_elmer_files(&fe_mesh, cur_output_folder);
		if (IO::file_exists(cur_output_folder + "/mesh.header")) cout << "FESS: Exported new FE mesh to subfolder.\n";
		else cout << "FESS: ERROR: Failed to export new FE mesh to subfolder.\n";

		// Call Elmer to run FEA on new FE mesh
		string batch_file = mesher::create_batch_file(cur_output_folder);
		cout << "FESS: Calling Elmer .bat file (" << batch_file << ")\n";
		physics::call_elmer(batch_file);

		// Obtain vonmises stress distribution from the .vtk file
		double* vonmises = new double[(grid.x) * (grid.y)]; // Nodes grid has +1 width along each dimension
		string elmer_output_file = cur_output_folder + "/case0001.vtk";
		physics::FEResults2D fe_results(grid);
		physics::load_2d_physics_data(elmer_output_file, fe_results, grid, mesh.offset, "Vonmises");
		cout << "FESS: Finished reading physics data." << endl;

		// Get minimum and maximum stress values
		max_stress = fe_results.max;
		min_stress = fe_results.min;
		cout << "FESS: Current maximum stress: " << std::setprecision(3) << std::scientific << max_stress << endl;
		cout << "FESS: Current minimum stress: " << std::setprecision(3) << std::scientific << min_stress << endl;
		//max_stress = 1e12; // Temporary value that will make FESS terminate after one iteration. TODO: remove

		// Check termination conditions
		bool terminate = false;
		if (max_stress > max_stress_threshold) {
			cout << "\nFESS: highest stress in FE result (" << std::setprecision(3) << std::scientific << max_stress << ") EXCEEDS MAXIMUM THRESHOLD (" << std::setprecision(3) << std::scientific << max_stress_threshold << ")\n";
			terminate = true;
		}
		if (min_stress > min_stress_threshold) {
			cout << "\nFESS: lowest stress in FE result (" << std::setprecision(3) << std::scientific << min_stress << ") EXCEEDS MINIMUM THRESHOLD (" << std::setprecision(3) << std::scientific << min_stress_threshold << ")\n";
			terminate = true;
		}
		if (terminate) {
			cout << "Terminating FESS algorithm after " << i << " iterations. Results were saved to " << output_folder << endl;
			break;
		}

		// If termination conditions were not reached, prepare density distribution for next iteration by removing 
		// elements below minimum stress threshold
		int no_cells_removed = physics::remove_low_stress_cells(fe_results.values, densities, min_stress_threshold, grid.x, grid.y);
		cout << "FESS: Removed low-stress cells. Relative volume decreased by " << (float)no_cells_removed / (float)grid.size2d << "\n";

		i++;
	}
}
