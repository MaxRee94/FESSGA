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
		phys::FEACase _fea_case, msh::SurfaceMesh _mesh, string _base_folder, double _min_stress_threshold,
		grd::Densities2d _densities, int _max_iterations, float _greediness,
		bool _maintain_boundary_connection, bool _export_msh = false, bool _verbose = true
	) : OptimizerBase(_fea_case, _mesh, _base_folder, _densities, _max_iterations, _export_msh, _verbose)
	{
		min_stress_threshold = _min_stress_threshold;
		greediness = _greediness;
		fea_case = _fea_case;
		densities.fea_case->maintain_boundary_connection = _maintain_boundary_connection;
	}
	double min_stress_threshold = 1.0;
	float greediness;
	void run();
	void log_termination(string final_valid_iteration_folder, int final_valid_iteration);
	int FESS::handle_floating_pieces(
		msh::FEMesh2D* fe_mesh, int no_cells_to_remove, int no_cells_removed, bool recurse = true
	);
	void export_stats(string iteration_name) {
		export_base_stats(iteration_name);
	}
};


void FESS::log_termination(string final_valid_iteration_folder, int final_valid_iteration) {
	cout << "\nTerminating FESS algorithm after " << final_valid_iteration << " iterations. Final results were saved to "
		<< final_valid_iteration_folder << endl;
}

int FESS::handle_floating_pieces(msh::FEMesh2D* fe_mesh, int no_cells_to_remove, int no_cells_removed, bool recurse) {
	vector<int> visited_cells;
	densities.init_pieces();
	if (densities.pieces.size() > 1) {
		// Attempt removal of smaller floating pieces
		cout << "FESS: Element removal resulted in a shape consisting of " << densities.pieces.size() << " pieces. Attempting to remove smaller piece(s)..." << endl;
		cout << "piece sizes: "; for (auto& piece : densities.pieces) cout << piece.cells.size() << "  "; cout << endl;
		
		//for (auto& piece : densities.pieces) {
		//	for (auto& cell : piece.cells) {
		//		densities.set(cell, piece.is_main_piece ? 5 : 4);
		//		//if (!piece.is_main_piece) cout << "cell outside main piece: " << cell / densities.dim_y << ", " << cell % densities.dim_y << endl;
		//	};
		//}
		//densities.print();

		densities.remove_smaller_pieces(
			densities.pieces, &densities.removed_cells);

		// Check if all pieces were removed
		int total_no_cells;
		if (densities.pieces.size() == 1) {
			cout << "FESS: All floating pieces successfully removed." << endl;
		}
		else {
			// If not all the floating pieces could be removed, restore the removed cells and re-try cell removal in so-called 'careful mode'.
			// This means that cells are only removed if they will not result in the splitting of the shape into multiple pieces.
			cout << "FESS: Not all smaller pieces could be removed.\n";
			densities.restore_removed_cells(densities.removed_cells);
			if (!densities.is_single_piece()) {
				cout << "Restoring individual cells did not result in unity of the shape. Also restoring removed pieces..\n";
				densities.restore_removed_pieces(densities.removed_pieces);
			}
			int no_cells_removed = densities.get_no_cells_in_removed_pieces();
			cout << "no cells in removed pieces: " << no_cells_removed << endl;
			no_cells_to_remove -= no_cells_removed;

			// If insufficient cells have been removed, remove them one-by-one (each time checking whether the shape still consists of one piece)
			if (no_cells_to_remove > 0) {
				cout << "Re-trying cell removal in 'careful mode'\n";
				no_cells_removed += densities.remove_low_stress_cells(no_cells_to_remove, no_cells_removed, &densities.pieces.at(0));
			}
			
			//int no_removed_pieces = densities.removed_pieces.size();
			/*if (!densities.is_single_piece()) {
				cout << "Restoring individual cells did not result in unity of the shape. Also restoring removed pieces..\n";
				densities.restore_removed_pieces(densities.removed_pieces);
			}*/
			
			//if (no_removed_pieces > 0) {
			//	// At least one piece was succesfully removed previously. These should again be removed.
			//	densities.remove_smaller_pieces(
			//		densities.fea_case, densities.removed_pieces, &densities.removed_cells, max_stress_threshold, densities.fea_results,
			//		maintain_boundary_connection, false, true
			//	);
			//	cout << "Densities after attempt to re-remove pieces:\n";
			//	densities.print();

			//	if (!densities.is_single_piece()) {
			//		cout << "FESS: Something went wrong with piece removal (shape is composed of multiple pieces again). Undoing removal.." << endl;
			//		densities.restore_removed_pieces(previously_removed_pieces); // Undo piece removal if something went wrong
			//	}
			//	else cout << "FESS: " << densities.removed_pieces.size() << " pieces were succesfully removed." << endl;
			//}
			//no_cells_to_remove -= densities.removed_cells.size() + densities.get_no_cells_in_removed_pieces();
		}
	}
	return no_cells_removed;
}


void FESS::run() {
	cout << "Beginning FESS run. Saving results to " << output_folder << endl;
	double min_stress, max_stress;
	string final_valid_iteration_folder;
	int final_valid_iteration = 1;
	int i = 1;
	bool last_iteration_was_valid = true;
	if (verbose) densities.print();
	while (i - 1 < max_iterations) {
		cout << "\nFESS: Starting iteration " << i << ".\n";

		// Create new subfolder for output of current iteration
		iteration_folder = get_iteration_folder(i, true);
		if (last_iteration_was_valid) {
			final_valid_iteration_folder = iteration_folder;
			export_stats(iteration_name);
		}

		// Reset densities object (keeping only the density values themselves)
		grd::Densities2d _densities = grd::Densities2d(densities.dim_x, mesh.diagonal, iteration_folder);
		_densities.copy_from(&densities);
		densities = _densities;

		// Generate new FE mesh using modified density distribution
		cout << "FESS: Generating new FE mesh...\n";
		msh::FEMesh2D fe_mesh;
		msh::create_FE_mesh(mesh, densities, fe_mesh);
		cout << "FESS: FE mesh generation done.\n";

		// Create and export a new version of the case.sif file by updating the boundary ids to fit the topology of the current FE mesh
		map<string, vector<int>> bound_id_lookup;
		msh::create_bound_id_lookup(&fea_case.bound_cond_lines, &fe_mesh, bound_id_lookup);
		msh::assemble_fea_case(densities.fea_case, &bound_id_lookup);
		IO::write_text_to_file(densities.fea_case->content, iteration_folder + "/case.sif");
		cout << "FESS: Exported updated case.sif file.\n";

		// Export newly generated FE mesh
		msh::export_as_elmer_files(&fe_mesh, iteration_folder);
		if (export_msh) msh::export_as_msh_file(&fe_mesh, iteration_folder);
		if (IO::file_exists(iteration_folder + "/mesh.header")) cout << "FESS: Exported new FE mesh.\n";
		else cout << "FESS: ERROR: Failed to export new FE mesh.\n";

		// Export density distribution
		string densities_file = densities.do_export(iteration_folder + "/distribution2d.dens");
		cout << "FESS: Exported current density distribution.\n";

		// Call Elmer to run FEA on new FE mesh
		string batch_file = msh::create_batch_file(iteration_folder);
		cout << "FESS: Calling Elmer .bat file...\n";
		FILE* pipe;
		fessga::phys::call_elmer(batch_file);
		cout << "FESS: ElmerSolver finished. Attempting to read .vtk file...\n";

		// Obtain vonmises stress distribution from the .vtk file
		string cur_case_output_file = iteration_folder + "/case0001.vtk";
		if (!IO::file_exists(cur_case_output_file)) {
			cout << "\nFESS: ERROR: Elmer did not produce a .vtk file (expected path " << cur_case_output_file << ")\n";
			cout << "FESS: Terminating program." << endl;
			exit(1);
		}
		bool physics_loaded = fessga::phys::load_2d_physics_data(
			cur_case_output_file, &densities.fea_results, densities.dim_x, densities.dim_y, densities.cell_size, mesh.offset, "Vonmises"
		);
		if (physics_loaded) cout << "FESS: Finished reading stress distribution from .vtk file." << endl;
		else {
			cout << "FESS: Error: Unable to read physics data from file " << cur_case_output_file << endl;
		}

		// Get minimum and maximum stress values
		max_stress = densities.fea_results.max;
		min_stress = densities.fea_results.min;
		cout << "FESS: Current maximum stress: " << std::setprecision(3) << std::scientific << max_stress << endl;
		cout << "FESS: Current minimum stress: " << std::setprecision(3) << std::scientific << min_stress << endl;

		// Check if maximum stress exceeds threshold
		if (max_stress > fea_case.max_stress_threshold) {
			cout << "FESS: highest stress in FE result (" << std::setprecision(3) << std::scientific << max_stress
				<< ") EXCEEDS MAXIMUM THRESHOLD (" << std::setprecision(3) << std::scientific << fea_case.max_stress_threshold << ")\n";
			final_valid_iteration_folder = get_iteration_folder(i - 1);
			log_termination(final_valid_iteration_folder, final_valid_iteration);
			break;
		}

		// If termination conditions were not met, prepare density distribution for next iteration by removing 
		// elements below minimum stress threshold.
		int no_cells_to_remove = max(1, (int)round(greediness * (float)densities.count()));
		densities.remove_low_stress_cells(no_cells_to_remove, 0);

		// Check if the resulting shape consists of exactly one piece. If not, the shape has become invalid and a repair operation is necessary.
		int no_cells_removed = handle_floating_pieces(&fe_mesh, no_cells_to_remove, densities.removed_cells.size());
		if (last_iteration_was_valid && verbose) densities.print();

		// Check if at least one cell was actually removed. If not, there must be insufficient cells left to continue optimization.
		last_iteration_was_valid = false;
		if (no_cells_removed == 0) {
			cout << "FESS: Termination condition reached: Unable to remove any more cells." << endl;
			log_termination(final_valid_iteration_folder, final_valid_iteration);
			break;
		}
		else {
			cout << "FESS: Removed " << no_cells_removed << " low - stress cells. Relative volume decreased by " << std::fixed
				<< (float)no_cells_removed / (float)densities.size << ", to "
				<< (float)(densities.count()) / (float)densities.size << "\n";
			last_iteration_was_valid = true;
			final_valid_iteration = i;
			densities.save_snapshot();
		}

		i++;
	}
}
