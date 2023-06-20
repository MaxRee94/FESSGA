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
	string final_valid_iteration_folder = cur_output_folder;
	int final_valid_iteration = 0;
	int i = 1;
	bool terminate = false;
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
		if (IO::file_exists(cur_output_folder + "/mesh.header")) cout << "FESS: Exported new FE mesh.\n";
		else cout << "FESS: ERROR: Failed to export new FE mesh.\n";

		// Export density distribution
		string densities_file = mesher::export_density_distrib(cur_output_folder, densities, grid.x, grid.y);
		cout << "FESS: Exported current density distribution.\n";
		mesher::print_density_distrib(densities, grid.x, grid.y);

		// Call Elmer to run FEA on new FE mesh
		string batch_file = mesher::create_batch_file(cur_output_folder);
		cout << "FESS: Calling Elmer .bat file...\n";
		physics::call_elmer(batch_file);
		cout << "FESS: ElmerSolver finished. Attempting to read .vtk file...\n";

		// Obtain vonmises stress distribution from the .vtk file
		double* vonmises = new double[(grid.x) * (grid.y)]; // Nodes grid has +1 width along each dimension
		string cur_case_output_file = cur_output_folder + "/case0001.vtk";
		if (!IO::file_exists(cur_case_output_file)) {
			cout << "\nFESS: ERROR: Elmer did not produce a .vtk file (expected path " << cur_case_output_file << ")\n";
			cout << "FESS: Terminating program." << endl;
			exit(1);
		}
		physics::FEResults2D fe_results(grid);
		bool physics_loaded = physics::load_2d_physics_data(cur_case_output_file, fe_results, grid, mesh.offset, "Vonmises");
		if (physics_loaded) cout << "FESS: Finished reading stress distribution from .vtk file." << endl;
		else {
			cout << "FESS: Error: Unable to read physics data from file " << cur_case_output_file << endl;
		}

		// Get minimum and maximum stress values
		max_stress = fe_results.max;
		min_stress = fe_results.min;
		cout << "FESS: Current maximum stress: " << std::setprecision(3) << std::scientific << max_stress << endl;
		cout << "FESS: Current minimum stress: " << std::setprecision(3) << std::scientific << min_stress << endl;

		// Check termination conditions
		if (max_stress > max_stress_threshold) {
			cout << "FESS: highest stress in FE result (" << std::setprecision(3) << std::scientific << max_stress
				<< ") EXCEEDS MAXIMUM THRESHOLD (" << std::setprecision(3) << std::scientific << max_stress_threshold << ")\n";
			terminate = true;
		}
		if (terminate) {
			cout << "\nTerminating FESS algorithm after " << final_valid_iteration << " iterations. Final results were saved to "
				<< final_valid_iteration_folder << endl;
			break;
		}

		// If termination conditions were not reached, prepare density distribution for next iteration by removing 
		// elements below minimum stress threshold. TODO: Prevent removal of elements to which boundary conditions were applied.
		int no_cells_to_remove = max(1, (int)round(greediness * (float)fe_mesh.surfaces.size()));
		vector<int> removed_cells;
		int no_cells_removed = physics::remove_low_stress_cells(
			&fe_results.data, densities, &fe_case, grid, no_cells_to_remove, removed_cells, max_stress_threshold
		);
		int total_no_cells = fe_mesh.surfaces.size() - no_cells_removed;
		cout << "total no cells (in fess): " << total_no_cells << endl;

		// Check if the resulting shape consists of exactly one piece. If not, the shape has become invalid and a repair operation is necessary.
		int no_pieces = 1;
		vector<mesher::Piece> pieces;
		vector<int> visited_cells;
		mesher::get_pieces(densities, grid, &fe_case, &pieces, &visited_cells, total_no_cells, &removed_cells, no_pieces);
		if (no_pieces != 1) {
			cout << "FESS: Element removal resulted in a shape consisting of multiple pieces. Attempting to remove smaller piece(s)..." << endl;
			vector<int> unremoved_piece_indices = physics::remove_smaller_pieces(
				densities, grid, &fe_case, total_no_cells, &pieces, &removed_cells, max_stress_threshold, &fe_results
			);
			if (unremoved_piece_indices.size() == 0) cout << "FESS: All floating pieces successfully removed." << endl;
			else if (unremoved_piece_indices.size() == pieces.size()) {
				cout << "FESS: Smaller pieces could not be removed. Restoring removed cells and re-trying cell removal in 'careful mode'...\n";
				// If none of the pieces could be removed, restore the removed cells and re-try cell removal in so-called 'careful mode'.
				// This means that cells are only removed if they will not result in the splitting of the shape into multiple pieces.
				mesher::restore_removed_cells(densities, grid, &removed_cells);
				removed_cells.clear();
				physics::remove_low_stress_cells(
					&fe_results.data, densities, &fe_case, grid, no_cells_to_remove, removed_cells, max_stress_threshold,
					&pieces[0], fe_mesh.surfaces.size()
				);
				int total_no_cells = fe_mesh.surfaces.size() - no_cells_removed;
			}
		}

		// Check if at least one cell was actually removed. If not, there must be insufficient cells left to continue optimization.
		terminate = false;
		if (no_cells_removed == 0) {
			cout << "FESS: Termination condition reached: Unable to remove any more cells." << endl;
			terminate = true;
		}
		else {
			cout << "FESS: Removed " << no_cells_removed << " low - stress cells. Relative volume decreased by " << std::fixed
				<< (float)no_cells_removed / (float)grid.size2d << ", to "
				<< (float)(total_no_cells) / (float)grid.size2d << "\n";
			final_valid_iteration_folder = cur_output_folder;
			final_valid_iteration++;
		}


		i++;
	}
}
