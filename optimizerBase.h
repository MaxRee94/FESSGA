#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <map>
#include "images.h"

using namespace fessga;


class OptimizerBase {
public:
	OptimizerBase() = default;

	// Constructor for 2d optimization
	OptimizerBase(
		string _msh_file, string _fea_case, msh::SurfaceMesh _mesh, string _base_folder,
		double _max_stress_threshold, grd::Densities2d _densities, int _max_iterations, bool _export_msh, bool _verbose
	) {
		mesh = _mesh;
		msh_file = _msh_file; base_folder = IO::get_fullpath(_base_folder);
		output_folder = IO::get_unique_path(base_folder + "/run_#");
		densities = _densities;
		export_msh = _export_msh;
		max_stress_threshold = _max_stress_threshold;
		domain_2d = true;
		no_cells = densities.size;
		max_iterations = _max_iterations;
		verbose = _verbose;
		IO::create_folder_if_not_exists(output_folder);
		image_folder = IO::create_folder_if_not_exists(output_folder + "/image_output");
		msh::derive_boundary_conditions(densities, bound_conds, mesh);
	};
	msh::SurfaceMesh mesh;
	bool domain_2d = false;
	int no_cells = 1;
	double max_stress_threshold = 0.0;
	grd::Densities2d densities;
	bool verbose = true;
	int max_iterations = 0;
	int iteration_number = 0;
	bool export_msh = false;
	double min_stress, max_stress;
	map<string, vector<pair<int, int>>> bound_conds;
	string msh_file, base_folder, image_folder, iteration_name, iteration_folder, output_folder;

	// Function to get the folder corresponding to the given iteration number. If the folder does not exist yet, it will be created.
	string get_iteration_folder(int iteration, bool verbose = false) {
		iteration_name = fessga::help::add_padding("iteration_", iteration) + to_string(iteration);
		iteration_folder = output_folder + "/" + iteration_name;
		if (verbose) cout << "OptimizerBase: Creating folder " << iteration_folder << " for current iteration if it does not exist yet...\n";
		IO::create_folder_if_not_exists(iteration_folder);

		return iteration_folder;
	}

	// Write densities to image file
	virtual void write_densities_to_image() {
		img::write_distribution_to_image(densities, image_folder + "/" + iteration_name + ".jpg", true);
	}

	virtual void export_meta_parameters(vector<string>* additional_metaparameters) {
		vector<string> _content = {
			"max stress threshold = " + to_string(max_stress_threshold),
		};
		help::append_vector(_content, additional_metaparameters);
		string content = help::join(&_content, "\n");
		IO::write_text_to_file(content, output_folder + "/metaparameters.txt");
	}

	// Function to write statistics, images, and other generated data to files
	void export_base_stats(string iteration_name) {
		write_densities_to_image();
	}

	void load_physics(grd::Densities2d* densities, bool verbose = false) {
		string cur_case_output_file = densities->output_folder + "/case0001.vtk";
		if (!IO::file_exists(cur_case_output_file)) {
			cout << "\nOptimizerBase: ERROR: Elmer did not produce a .vtk file (expected path " << cur_case_output_file << ")\n";
			cout << "OptimizerBase: Terminating program." << endl;
			exit(1);
		}
		bool physics_loaded = fessga::phys::load_2d_physics_data(
			cur_case_output_file, &densities->fea_results, densities->dim_x, densities->dim_y, densities->cell_size, mesh.offset, "Vonmises"
		);
		if (physics_loaded && verbose) cout << "OptimizerBase: Finished reading stress distribution from .vtk file." << endl;
		else if (!physics_loaded) {
			cout << "OptimizerBase: Error: Unable to read physics data from file " << cur_case_output_file << endl;
		}
	}

	// Copy all files in a solution folder (iteration folder for FESS, individual folder for EVOMA) to the specified target folder
	void copy_solution_files(string source_dir, string target_dir, bool verbose = true) {
		if (verbose) cout << "Copying solution files from directory " << source_dir << " to " << target_dir << endl;
		vector<string> filepaths;
		IO::get_files_in_directory(filepaths, source_dir);
		for (auto& fpath : filepaths) {
			string target_path = help::replace_occurrences(fpath, source_dir, target_dir);
			IO::copy_file(fpath, target_path);
		}
	}

};


