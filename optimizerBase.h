#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <map>
#include "images.h"
#include <list>

using namespace fessga;


bool load_physics(grd::Densities2d* densities, msh::SurfaceMesh* mesh, vector<int>* times, bool verbose = false);

class OptimizerBase {
public:
	OptimizerBase() = default;

	// Constructor for 2d optimization
	OptimizerBase(
		phys::FEACaseManager _fea_casemanager, msh::SurfaceMesh _mesh, string _base_folder,
		grd::Densities2d _densities, int _max_iterations, bool _export_msh, bool _verbose, string existing_folder = ""
	) {
		mesh = _mesh;
		base_folder = IO::get_fullpath(_base_folder);
		if (existing_folder == "") {
			output_folder = IO::get_unique_path(base_folder + "/run_#");
		}
		else {
			output_folder = existing_folder;
		}
		densities = _densities;
		export_msh = _export_msh;
		domain_2d = true;
		no_cells = densities.size;
		max_iterations = _max_iterations;
		verbose = _verbose;
		initial_count = _densities.count();
		IO::create_folder_if_not_exists(output_folder);
		image_folder = IO::create_folder_if_not_exists(output_folder + "/image_output");
		fea_casemanager = _fea_casemanager;
		densities.fea_casemanager = &fea_casemanager;
		statistics_file = output_folder + "/statistics.csv";
		densities.visualize_keep_cells();

		// Create 'copy_superpositions' batfile in output folder
		IO::write_text_to_file("python ../../../../scripts/copy_superpositions.py\npause", output_folder + "/copy_superpositions.bat");

		// Create 'superpositions' folder in output folder
		IO::create_folder_if_not_exists(output_folder + "/superpositions");
	};
	msh::SurfaceMesh mesh;
	bool domain_2d = false;
	int no_cells = 1;
	grd::Densities2d densities;
	bool verbose = true;
	int max_iterations = 0;
	int iteration_number = 0;
	int initial_count = 0;
	bool export_msh = false;
	phys::FEACaseManager fea_casemanager;
	string statistics_file;
	list<string> stats;
	double min_stress, max_stress;
	time_t start_time;
	bool initialize = true;
	string base_folder, image_folder, iteration_name, iteration_folder, output_folder;

	// Function to get the folder corresponding to the given iteration number. If the folder does not exist yet, it will be created.
	string get_iteration_folder(int iteration, bool verbose = false, bool load_existing_population = false, string existing_population = "") {
		/*if (load_existing_population) {
			vector<string> _existing_pop;
			help::split(existing_population, "/", _existing_pop);
			iteration_name = _existing_pop[_existing_pop.size() - 1];
			iteration_folder = existing_population;
		}
		else {
			iteration_name = fessga::help::add_padding("iteration_", iteration) + to_string(iteration);
			iteration_folder = output_folder + "/" + iteration_name;
		}*/
		iteration_name = fessga::help::add_padding("iteration_", iteration) + to_string(iteration);
		iteration_folder = output_folder + "/" + iteration_name;
		if (verbose) {
			cout << "OptimizerBase: Creating folder " << iteration_folder <<
				" for current iteration if it does not exist yet...\n";
		}
		IO::create_folder_if_not_exists(iteration_folder);

		return iteration_folder;
	}

	// Write densities to image file
	virtual void write_densities_to_image(bool verbose = false) {
		img::write_distribution_to_image(densities, image_folder + "/" + iteration_name + ".jpg", true);
	}

	virtual void export_meta_parameters(vector<string>* additional_metaparameters) {
		vector<string> _content = {
			"mechanical threshold = " + to_string(fea_casemanager.mechanical_threshold),
			"max compressive strength = " + to_string(fea_casemanager.max_compressive_strength),
			"max tensile strength = " + to_string(fea_casemanager.max_tensile_strength)
		};
		help::append_vector(_content, additional_metaparameters);
		string content = help::join(&_content, "\n");
		string outpath = output_folder + "/metaparameters.txt";
		if (IO::file_exists(outpath)) outpath = IO::get_unique_file_path(outpath);
		IO::write_text_to_file(content, outpath);
	}

	// Function to write statistics, images, and other generated data to files (called once per iteration)
	void export_base_stats() {
		write_densities_to_image();

		// Collect memory stats
		vector<float> memory_status = help::get_free_memory();
		string memory_string = help::join_as_string(memory_status, ", ");
		stats.push_back(memory_string);

		// Collect other universal stats
		float iteration_time = 0;
		if (!initialize) iteration_time = difftime(time(0), start_time);
		stats.push_front(to_string(iteration_time));
		stats.push_front(to_string(iteration_number));

		// Append the stats to the stats file
		IO::append_to_file(statistics_file, help::join(&stats, ", "));
		stats.clear();
		start_time = time(0);
	}

	// Copy all files in a solution folder (iteration folder for FESS, individual folder for emma) to the specified target folder
	void copy_solution_files(string source_dir, string target_dir, bool verbose = false) {
		if (verbose) cout << "Copying solution files from directory " << source_dir << " to " << target_dir << endl;
		vector<string> filepaths;
		IO::get_files_in_directory(filepaths, source_dir);
		if (verbose) cout << "Starting copy process...\n";
		for (auto& fpath : filepaths) {
			if (verbose) cout << "Copying file " << fpath << endl;
			if (!IO::file_exists(fpath)) cout << "ERROR: File " << fpath << " does not exist.\n";
			string target_path = help::replace_occurrences(fpath, source_dir, target_dir);
			IO::copy_file(fpath, target_path);
		}
	}

	// Create and export new versions of the case.sif files by updating the boundary ids to fit the topology of the current FE mesh
	void create_sif_files(grd::Densities2d* densities, msh::FEMesh2D* fe_mesh, bool verbose = false) {
		densities->fea_casemanager->update_casepaths(densities->output_folder);
		for (auto& fea_case : densities->fea_casemanager->active_cases) {
			map<string, vector<int>> bound_id_lookup;
			msh::create_bound_id_lookup(&fea_case.bound_cond_lines, fe_mesh, bound_id_lookup);
			msh::assemble_fea_case(densities->fea_casemanager, &fea_case, &bound_id_lookup);
			IO::write_text_to_file(fea_case.content, fea_case.path);
			if (verbose) cout << "-- Exported case file to path " << fea_case.path << endl;
		}
		if (verbose) cout << "OptimizerBase: Exported updated case.sif files.\n";
	}

	// Run FEA on all vtk cases of a single solution
	bool run_FEA_on_single_solution(string individual_folder, bool is_2nd_attempt = false) {
		string elmer_bat_file = individual_folder + "/run_elmer.bat";
		Timer timer; timer.start();
		while (!IO::file_exists(elmer_bat_file)) {} // Wait for elmer batfile to appear on disk
		timer.stop();
		if (timer.elapsedSeconds() > 1) cout << "WARNING: Waited for a long time for elmer batfile '" << elmer_bat_file << "' to appear on disk (" << timer.elapsedMilliseconds() << " ms)\n";
		fessga::phys::call_elmer(individual_folder, &fea_casemanager);

		// Get vtk paths
		vector<string> vtk_paths;
		msh::get_vtk_paths(&fea_casemanager, individual_folder, vtk_paths);
		/*if (!is_2nd_attempt) IO::remove_file(vtk_paths[0]);
		cout << "Removed file " << vtk_paths[0] << endl;*/

		// Wait for all case files to appear on disk
		time_t start = time(0);
		bool fea_failed = false;
		for (auto& vtk_path : vtk_paths) {
			while (!IO::file_exists(vtk_path)) {
				float seconds_since_start = difftime(time(0), start);
				if (seconds_since_start > 10) {
					// If the vtk file has not appeared after 10 seconds, retry running the elmer batchfile (once).
					fea_failed = true;
					if (!is_2nd_attempt) {
						cout << "- WARNING:    First attempt to produce a vtk file failed on case '" << vtk_path << "'\n";
					}
					else {
						cout << "- ERROR:    Second attempt to produce a vtk file failed on case '" << vtk_path << "'\n";
					}
					break;
				}
			}
			if (fea_failed) break;
		}
		if (fea_failed) {
			if (is_2nd_attempt) { // After two attempts to run the FEA, consider FEA to have failed
				return false;
			}
			cout << "Re-running Elmer bat file in individual folder '" << individual_folder << "'\n";
			return run_FEA_on_single_solution(individual_folder, true);
		}

		if (is_2nd_attempt) {
			cout << "- INFO:    Second attempt to produce a vtk file was successful (folder '" << individual_folder << "'\n";
		}
		return true;
	}
};


