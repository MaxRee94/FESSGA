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
#include "images.h"



class FESS : public OptimizerBase {
public:
	FESS() = default;
	FESS(
		string _msh_file, string _fe_case, mesher::SurfaceMesh _mesh, string _output_folder, double _min_stress_threshold,
		double _max_stress_threshold, uint* _starting_densities, mesher::Grid3D _grid, int _max_iterations, float _greediness,
		bool _maintain_boundary_cells, bool _export_msh = false, bool _verbose = true
	) : OptimizerBase(_msh_file, _fe_case, _mesh, _output_folder, _max_stress_threshold, _starting_densities, _grid, _max_iterations, _export_msh, _verbose)
	{
		min_stress_threshold = _min_stress_threshold;
		greediness = _greediness;
		maintain_boundary_cells = _maintain_boundary_cells;
	}
	double min_stress_threshold = 1.0;
	float greediness;
	bool maintain_boundary_cells = true;
	void run();
	void log_termination(string final_valid_iteration_folder, int final_valid_iteration);
	void FESS::handle_floating_pieces(
		uint* densities, mesher::Grid3D grid, mesher::Case* fe_case, int& total_no_cells, vector<int>* removed_cells, physics::FEResults2D* fe_results,
		int& no_cells_removed, mesher::FEMesh2D* fe_mesh, int& no_cells_to_remove, bool recurse = true
	);
};


void FESS::log_termination(string final_valid_iteration_folder, int final_valid_iteration) {
	cout << "\nTerminating FESS algorithm after " << final_valid_iteration << " iterations. Final results were saved to "
		<< final_valid_iteration_folder << endl;
}

void flush_whitelist(mesher::Case* fe_case, bool& whitelist_flushed) {
	cout << "FESS: Re-trying optimization after flushing whitelist." << endl;
	fe_case->whitelisted_cells.clear();
	whitelist_flushed = true;
}

void FESS::handle_floating_pieces(
	uint* densities, mesher::Grid3D grid, mesher::Case* fe_case, int& total_no_cells, vector<int>* removed_cells, physics::FEResults2D* fe_results,
	int& no_cells_removed, mesher::FEMesh2D* fe_mesh, int& no_cells_to_remove, bool recurse
) {
	int no_pieces = 1;
	vector<mesher::Piece> pieces;
	vector<int> visited_cells;
	mesher::get_pieces(densities, grid, fe_case, &pieces, &visited_cells, total_no_cells, removed_cells, no_pieces);
	if (no_pieces != 1) {
		// Attempt removal of smaller floating pieces
		cout << "FESS: Element removal resulted in a shape consisting of " << no_pieces << " pieces. Attempting to remove smaller piece(s)..." << endl;
		cout << "piece sizes: "; for (auto& piece : pieces) cout << piece.cells.size() << "  "; cout << endl;
		vector<int> removed_piece_indices = physics::remove_smaller_pieces(
			densities, grid, fe_case, total_no_cells, &pieces, removed_cells, max_stress_threshold, fe_results, maintain_boundary_cells
		);

		// Count no cells that were removed as the result of removing the floating pieces
		int no_cells_in_removed_pieces = 0;
		for (int j = 0; j < removed_piece_indices.size(); j++) no_cells_in_removed_pieces += pieces[removed_piece_indices[j]].cells.size();

		// Check if all pieces were removed
		if (removed_piece_indices.size() == pieces.size()) {
			cout << "FESS: All floating pieces successfully removed." << endl;
		}
		else {
			// If not all the floating pieces could be removed, restore the removed cells and re-try cell removal in so-called 'careful mode'.
			// This means that cells are only removed if they will not result in the splitting of the shape into multiple pieces.
			cout << "FESS: Not all smaller pieces could be removed. Restoring removed cells and re-trying cell removal in 'careful mode'...\n";
			mesher::restore_removed_cells(densities, grid, removed_cells);
			bool is_single_piece = mesher::is_single_piece(densities, grid, fe_case, total_no_cells, removed_cells);
			no_cells_removed = no_cells_in_removed_pieces;
			if (is_single_piece) {
				cout << "Also restoring removed pieces.\n";
				mesher::restore_removed_pieces(densities, &pieces);
				no_cells_in_removed_pieces = 0;
				no_cells_removed = 0;
			}
			total_no_cells = fe_mesh->surfaces.size() - no_cells_removed;
			vector<mesher::Piece> pieces_to_be_removed;
			if (removed_piece_indices.size() > 0) {
				// At least one piece was succesfully removed previously. These should again be removed.
				for (int j = 0; j < removed_piece_indices.size(); j++) pieces_to_be_removed.push_back(pieces[removed_piece_indices[j]]);
				removed_piece_indices = physics::remove_smaller_pieces(
					densities, grid, fe_case, total_no_cells, &pieces_to_be_removed, removed_cells, max_stress_threshold, fe_results, maintain_boundary_cells, false, true
				);

				// Count no cells that were removed as the result of removing the floating pieces
				no_cells_in_removed_pieces = 0;
				for (int j = 0; j < removed_piece_indices.size(); j++) no_cells_in_removed_pieces += pieces[removed_piece_indices[j]].cells.size();
				no_cells_removed = no_cells_in_removed_pieces;

				total_no_cells = fe_mesh->surfaces.size() - no_cells_in_removed_pieces;
				bool is_single_piece = mesher::is_single_piece(densities, grid, fe_case, total_no_cells, removed_cells);
				if (!is_single_piece) {
					cout << "FESS: Something went wrong with piece removal (shape is composed of multiple pieces again). Undoing removal.." << endl;
					mesher::restore_removed_pieces(densities, &pieces_to_be_removed); // Undo piece removal if something went wrong
					no_cells_removed = 0; total_no_cells = fe_mesh->surfaces.size();
				}
				cout << "FESS: " << removed_piece_indices.size() << " pieces were succesfully removed." << endl;
			}
			removed_cells->clear();
			no_cells_to_remove -= no_cells_removed;
			physics::remove_low_stress_cells(
				&fe_results->data, densities, fe_case, grid, no_cells_to_remove, *removed_cells, max_stress_threshold,
				&pieces[0], total_no_cells
			);
		}
		//cout << "densities after attempting to handle multiple pieces:\n";
		//mesher::print_density_distrib(densities, grid.x, grid.y);
		no_cells_removed = removed_cells->size() + no_cells_in_removed_pieces;
		total_no_cells = fe_mesh->surfaces.size() - no_cells_removed;

		// Recurse once to remove any floating pieces still remaining
		if (recurse) {
			no_cells_to_remove = 0;
			handle_floating_pieces(
				densities, grid, fe_case, total_no_cells, removed_cells, fe_results, no_cells_removed, fe_mesh, no_cells_to_remove, false
			);
		}
	}
}


void FESS::run() {
	cout << "Beginning FESS run. Saving results to " << output_folder << endl;

	// From the fe mesh, get a map<int, int> containing:
	//		* The names of each bound condition as the keys
	//		* Values which are themselves maps, containing [coord : bound_number] key-value pairs

	double min_stress, max_stress;
	string msh = msh_file;
	string cur_output_folder = output_folder;
	string final_valid_iteration_folder = cur_output_folder;
	int final_valid_iteration = 1;
	int i = 1;
	bool last_iteration_was_valid = true;
	bool whitelist_flushed = false;
	string image_folder = IO::create_folder_if_not_exists(output_folder + "/image_output");
	while (i - 1 < max_iterations) {
		cout << "\nFESS: Starting iteration " << i << ".\n";

		// Create new subfolder for output of current iteration
		string cur_iteration_name;
		string cur_output_folder = get_iteration_folder(i, cur_iteration_name, true);
		if (last_iteration_was_valid) {
			final_valid_iteration_folder = cur_output_folder;
			img::write_distribution_to_image(densities, grid, image_folder + "/" + cur_iteration_name + ".jpg", 1000, 1000, true);
		}

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
		if (last_iteration_was_valid && verbose) mesher::print_density_distrib(densities, grid.x, grid.y);

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

		// Check if maximum stress exceeds threshold
		if (max_stress > max_stress_threshold) {
			cout << "FESS: highest stress in FE result (" << std::setprecision(3) << std::scientific << max_stress
				<< ") EXCEEDS MAXIMUM THRESHOLD (" << std::setprecision(3) << std::scientific << max_stress_threshold << ")\n";
			final_valid_iteration_folder = get_iteration_folder(i - 1, cur_iteration_name);
			log_termination(final_valid_iteration_folder, final_valid_iteration);
			break;
		}

		// If termination conditions were not met, prepare density distribution for next iteration by removing 
		// elements below minimum stress threshold.
		int no_cells_to_remove = max(1, (int)round(greediness * (float)fe_mesh.surfaces.size()));
		vector<int> removed_cells;
		int no_cells_removed = physics::remove_low_stress_cells(
			&fe_results.data, densities, &fe_case, grid, no_cells_to_remove, removed_cells, max_stress_threshold
		);
		int total_no_cells = fe_mesh.surfaces.size() - no_cells_removed;

		// Check if the resulting shape consists of exactly one piece. If not, the shape has become invalid and a repair operation is necessary.
		handle_floating_pieces(densities, grid, &fe_case, total_no_cells, &removed_cells, &fe_results, no_cells_removed, &fe_mesh, no_cells_to_remove);

		// Check if at least one cell was actually removed. If not, there must be insufficient cells left to continue optimization.
		last_iteration_was_valid = false;
		if (no_cells_removed == 0 && !whitelist_flushed) {
			cout << "FESS: Unable to remove any more cells.\n";
			flush_whitelist(&fe_case, whitelist_flushed);
		}
		else if (no_cells_removed == 0 && whitelist_flushed) {
			cout << "FESS: Termination condition reached: Unable to remove any more cells." << endl;
			log_termination(final_valid_iteration_folder, final_valid_iteration);
			break;
		}
		else {
			cout << "FESS: Removed " << no_cells_removed << " low - stress cells. Relative volume decreased by " << std::fixed
				<< (float)no_cells_removed / (float)grid.size2d << ", to "
				<< (float)(total_no_cells) / (float)grid.size2d << "\n";
			whitelist_flushed = false;
			last_iteration_was_valid = true;
			final_valid_iteration = i;
		}

		i++;
	}
}
