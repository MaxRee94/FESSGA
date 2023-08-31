#pragma once
#include "fess.h"


void FESS::export_meta_parameters(vector<string>*_) {
	string feasibility_filtering = "use feasibility filtering = ";
	string suffix = do_feasibility_filtering ? "yes" : "no";
	feasibility_filtering += suffix;
	string bound_connection = "maintain boundary connection = ";
	suffix = fea_casemanager.maintain_boundary_connection ? "yes" : "no";
	bound_connection += suffix;
	vector<string> additional_metaparameters = {
		feasibility_filtering,
		"initial greediness = " + to_string(greediness),
		bound_connection
	};
	OptimizerBase::export_meta_parameters(&additional_metaparameters);
}

void FESS::export_stats(string iteration_name, bool initialize) {
	cout << "Exporting statistics to " << statistics_file << endl;
	if (initialize) IO::write_text_to_file(
		"Iteration, Iteration time, Relative area, Greediness, #Cells removed",
		statistics_file
	);
	stats.push_back(to_string(iteration_number));
	stats.push_back(to_string(relative_area));
	stats.push_back(to_string(greediness));
	stats.push_back(to_string(no_cells_removed));
	export_base_stats();
}

void FESS::log_termination(string final_valid_iteration_folder, int final_valid_iteration) {
	cout << "\nTerminating FESS algorithm after " << final_valid_iteration
		<< " iterations. Final results were saved to " << final_valid_iteration_folder << endl;
}

int FESS::handle_floating_pieces(msh::FEMesh2D* fe_mesh, int no_cells_to_remove, int no_cells_removed, bool recurse) {
	vector<int> visited_cells;
	densities.init_pieces();
	if (densities.pieces.size() > 1) {
		// Attempt removal of smaller floating pieces
		cout << "FESS: Element removal resulted in a shape consisting of " << densities.pieces.size()
			<< " pieces. Attempting to remove smaller piece(s)..." << endl;
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
			// If not all the floating pieces could be removed, restore the removed cells and re-try cell removal
			// in so-called 'careful mode'. This means that cells are only removed if they will not result in the
			// splitting of the shape into multiple pieces.
			cout << "FESS: Not all smaller pieces could be removed.\n";
			densities.restore_removed_cells(densities.removed_cells);
			if (!densities.is_single_piece()) {
				cout << "Restoring individual cells did not result in unity of the shape. "
					<< "Also restoring removed pieces..\n";
				densities.restore_removed_pieces(densities.removed_pieces);
			}
			int no_cells_removed = densities.get_no_cells_in_removed_pieces();
			cout << "no cells in removed pieces: " << no_cells_removed << endl;
			no_cells_to_remove -= no_cells_removed;

			// If insufficient cells have been removed, remove them one-by-one (each time checking whether the shape still consists of one piece)
			if (no_cells_to_remove > 0) {
				cout << "Re-trying cell removal in 'careful mode'\n";
				no_cells_removed += densities.remove_low_stress_cells(
					no_cells_to_remove, no_cells_removed, &densities.pieces.at(0)
				);
			}

			//int no_removed_pieces = densities.removed_pieces.size();
			/*if (!densities.is_single_piece()) {
				cout << "Restoring individual cells did not result in unity of the shape. Also restoring removed pieces..\n";
				densities.restore_removed_pieces(densities.removed_pieces);
			}*/

			//if (no_removed_pieces > 0) {
			//	// At least one piece was succesfully removed previously. These should again be removed.
			//	densities.remove_smaller_pieces(
			//		densities.fea_cases, densities.removed_pieces, &densities.removed_cells, max_stress_threshold, densities.fea_results,
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

void FESS::fill_design_domain() {
	densities.fill_all();
	for (int& cutout_cell : fea_casemanager.cutout_cells) {
		densities.del(cutout_cell);
	}
}


void FESS::run() {
	cout << "Beginning FESS run. Saving results to " << output_folder << endl;
	export_meta_parameters();
	double min_stress, max_stress;
	string final_valid_iteration_folder;
	int final_valid_iteration = 1;
	int iteration_number = 1;
	bool last_iteration_was_valid = true;
	if (verbose) densities.print();

	fill_design_domain();

	while (iteration_number - 1 < max_iterations) {
		cout << "\nFESS: Starting iteration " << iteration_number << ".\n";

		// Create new subfolder for output of current iteration
		iteration_folder = get_iteration_folder(iteration_number, true);
		if (last_iteration_was_valid) {
			if (verbose) densities.print();
			final_valid_iteration_folder = iteration_folder;
			relative_area = (float)(densities.count()) / (float)(densities.size);
			export_stats(iteration_name, iteration_number == 1);
		}
		densities.output_folder = iteration_folder;

		// Generate new FE mesh using modified density distribution
		cout << "FESS: Generating new FE mesh...\n";
		msh::FEMesh2D fe_mesh;
		msh::create_FE_mesh(mesh, densities, fe_mesh);
		cout << "FESS: FE mesh generation done.\n";

		create_sif_files(&densities, &fe_mesh);

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
		fessga::phys::call_elmer(iteration_folder, &fea_casemanager);
		cout << "FESS: ElmerSolver finished. Attempting to read .vtk file...\n";

		// Obtain vonmises stress distribution from the .vtk files
		load_physics(&densities, &mesh, verbose);

		// Get minimum and maximum stress values
		max_stress = densities.fea_results.max;
		min_stress = densities.fea_results.min;
		cout << "FESS: Current maximum stress: " << std::setprecision(3) << std::scientific << max_stress << endl;
		cout << "FESS: Current minimum stress: " << std::setprecision(3) << std::scientific << min_stress << endl;
		
		int no_cells_to_remove = max(1, (int)round(greediness * (float)densities.count()));

		// Check if maximum stress exceeds threshold
		if (max_stress > fea_casemanager.max_stress_threshold) {
			cout << std::setprecision(3) << std::scientific;
			cout << "FESS: Highest stress in FE result (" << max_stress
				<< ") exceeds maximum stress threshold (" << fea_casemanager.max_stress_threshold << ")\n";

			// Decrease greediness and load snapshot, so that optimization can be retried in a less agressive manner.
			if (no_cells_to_remove > 1) {
				float old_greediness = greediness;
				greediness /= 2.0;
				cout << "Re-trying optimization with reduced greediness...\n";
				cout << "Reducing greediness from " << old_greediness << " to " << greediness << endl;
				densities.load_snapshot();
				no_cells_to_remove = max(1, (int)round(greediness * (float)densities.count()));
			}
			else {
				// If greediness was already at a minimal level (meaning only single cells were removed in 
				// each iteration), terminate FESS altogether.
				final_valid_iteration_folder = get_iteration_folder(iteration_number - 1);
				log_termination(final_valid_iteration_folder, final_valid_iteration);
				break;
			}
		}

		// If termination conditions were not met, prepare density distribution for next iteration by removing 
		// elements with lowest stress
		densities.save_snapshot();
		densities.removed_cells.clear();
		densities.redo_count();
		int initial_no_cells = densities.count();
		no_cells_removed = densities.remove_low_stress_cells(no_cells_to_remove, 0);

		// Check if the resulting shape consists of exactly one piece. 
		// If not, the shape has become invalid and a repair operation is necessary.
		no_cells_removed = handle_floating_pieces(&fe_mesh, no_cells_to_remove, no_cells_removed);
		densities.redo_count();

		// Do feasibility filtering if applicable
		if (do_feasibility_filtering) densities.do_feasibility_filtering();
		no_cells_removed = initial_no_cells - densities.count();

		// Check if at least one cell was actually removed. If not, there must be insufficient cells left to continue
		// optimization.
		last_iteration_was_valid = false;
		if (no_cells_removed == 0) {
			cout << "FESS: Termination condition reached: Unable to remove any more cells." << endl;
			log_termination(final_valid_iteration_folder, final_valid_iteration);
			break;
		}
		else {
			cout << "FESS: Removed " << no_cells_removed << " low - stress cells. Relative volume decreased by "
				<< std::fixed << (float)no_cells_removed / (float)densities.size << ", to "
				<< (float)(densities.count()) / (float)densities.size << "\n";
			last_iteration_was_valid = true;
			final_valid_iteration = iteration_number;
		}

		densities.redo_count();
		iteration_number++;
	}
}