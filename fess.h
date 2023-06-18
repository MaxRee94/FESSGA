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
		string _msh_file, string _fe_case, mesher::SurfaceMesh _mesh, string _output_folder, double _min_stress_threshold,
		double _max_stress_threshold, uint* _starting_densities, mesher::Grid3D _grid, int _max_iterations, float _greediness, bool _export_msh = false
	) : OptimizerBase(_msh_file, _fe_case, _mesh, _output_folder, _max_stress_threshold, _starting_densities, _grid, _max_iterations, _export_msh)
	{
		min_stress_threshold = _min_stress_threshold;
		greediness = _greediness;
	}
	double min_stress_threshold = 1.0;
	float greediness;
	void run();
};


void FESS::run() {
	cout << "Beginning FESS run. Saving results to " << output_folder << endl;

	// From the fe mesh, get a map<int, int> containing:
	//		* The names of each bound condition as the keys
	//		* Values which are themselves maps, containing [coord : bound_number] key-value pairs

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

		// Generate new FE mesh using modified density distribution
		cout << "FESS: Generating new FE mesh...\n";
		mesher::FEMesh2D fe_mesh;
		mesher::generate_FE_mesh(grid, mesh, densities, fe_mesh);
		cout << "FESS: FE mesh generation done.\n";
		
		// Create and export a new version of the case.sif file by updating the boundary ids to fit the topology of the current FE mesh
		map<string, vector<int>> bound_id_lookup;
		mesher::create_bound_id_lookup(&bound_conds, &fe_mesh, bound_id_lookup);
		mesher::assemble_fe_case(&fe_case, &bound_id_lookup);
		IO::write_text_to_file(fe_case.content, cur_output_folder + "/case.sif");
		cout << "FESS: Exported updated case.sif file.\n";

		// Export newly generated FE mesh
		mesher::export_as_elmer_files(&fe_mesh, cur_output_folder);
		if (export_msh) mesher::export_as_msh_file(&fe_mesh, cur_output_folder);
		if (IO::file_exists(cur_output_folder + "/mesh.header")) cout << "FESS: Exported new FE mesh to subfolder.\n";
		else cout << "FESS: ERROR: Failed to export new FE mesh to subfolder.\n";

		// Export density distribution
		string densities_file = mesher::export_density_distrib(cur_output_folder, densities, grid.x, grid.y);
		cout << "FESS: Exported current density distribution to file: " << densities_file << endl;

		// Call Elmer to run FEA on new FE mesh
		string batch_file = mesher::create_batch_file(cur_output_folder);
		cout << "FESS: Calling Elmer .bat file (" << batch_file << ")\n";
		physics::call_elmer(batch_file);

		// Obtain vonmises stress distribution from the .vtk file
		double* vonmises = new double[(grid.x) * (grid.y)]; // Nodes grid has +1 width along each dimension
		string cur_case_output_file = cur_output_folder + "/case0001.vtk";
		if (!IO::file_exists(cur_case_output_file)) {
			cout << "\nFESS: ERROR: Elmer did not produce a .vtk file (expected path " << cur_case_output_file << ")\n";
			cout << "FESS: Terminating program." << endl;
			exit(1);
		}
		physics::FEResults2D fe_results(grid);
		physics::load_2d_physics_data(cur_case_output_file, fe_results, grid, mesh.offset, "Vonmises");
		cout << "FESS: Finished reading physics data." << endl;

		// Get minimum and maximum stress values
		max_stress = fe_results.max;
		min_stress = fe_results.min;
		cout << "FESS: Current maximum stress: " << std::setprecision(3) << std::scientific << max_stress << endl;
		cout << "FESS: Current minimum stress: " << std::setprecision(3) << std::scientific << min_stress << endl;

		// Check termination conditions
		bool terminate = false;
		if (max_stress > max_stress_threshold) {
			cout << "\nFESS: highest stress in FE result (" << std::setprecision(3) << std::scientific << max_stress << ") EXCEEDS MAXIMUM THRESHOLD (" << std::setprecision(3) << std::scientific << max_stress_threshold << ")\n";
			terminate = true;
		}
		if (terminate) {
			cout << "Terminating FESS algorithm after " << i << " iterations. Results were saved to " << output_folder << endl;
			break;
		}

		// If termination conditions were not reached, prepare density distribution for next iteration by removing 
		// elements below minimum stress threshold. TODO: Prevent removal of elements to which boundary conditions were applied.
		int no_cells_to_remove = max(3, (int)round(greediness * (float)fe_mesh.surfaces.size()));
		physics::remove_low_stress_cells(&fe_results.data, densities, min_stress_threshold, no_cells_to_remove);
		cout << "FESS: Removed " << no_cells_to_remove << " low - stress cells. Relative volume decreased by " << std::fixed
			<< (float)no_cells_to_remove / (float)grid.size2d << ", to "
			<< (float)(fe_mesh.surfaces.size() - no_cells_to_remove) / (float)(grid.size2d) << "\n";

		i++;
	}
}
