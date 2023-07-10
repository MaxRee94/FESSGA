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
		string _msh_file, string _fea_case, msh::SurfaceMesh _mesh, string _output_folder,
		double _max_stress_threshold, grd::Densities2d _densities, int _max_iterations, bool _export_msh, bool _verbose
	) {
		mesh = _mesh;
		msh_file = _msh_file; output_folder = IO::get_fullpath(_output_folder);
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
	int iteration_number = 1;
	bool export_msh = false;
	double min_stress, max_stress;
	map<string, vector<pair<int, int>>> bound_conds;
	string msh_file, output_folder, image_folder, iteration_name, iteration_folder;

	// Function to get the folder corresponding to the given iteration number. If the folder does not exist yet, it will be created.
	string get_iteration_folder(int iteration, bool verbose = false) {
		iteration_name = fessga::help::add_padding("iteration_", iteration) + to_string(iteration);
		string cur_output_folder = output_folder + "/" + iteration_name;
		if (verbose) cout << "FESS: Creating folder " << cur_output_folder << " for current iteration if it does not exist yet...\n";
		IO::create_folder_if_not_exists(cur_output_folder);

		return cur_output_folder;
	}

	// Write densities to image file
	virtual void write_densities_to_image() {
		img::write_distribution_to_image(densities, image_folder + "/" + iteration_name + ".jpg", true);
	}

	// Function to write statistics, images, and other generated data to files
	void export_base_stats(string iteration_name) {
		write_densities_to_image();
	}

	void load_physics(grd::Densities2d* densities) {
		string cur_case_output_file = densities->output_folder + "/case0001.vtk";
		if (!IO::file_exists(cur_case_output_file)) {
			cout << "\nOptimizerBase: ERROR: Elmer did not produce a .vtk file (expected path " << cur_case_output_file << ")\n";
			cout << "OptimizerBase: Terminating program." << endl;
			exit(1);
		}
		bool physics_loaded = fessga::phys::load_2d_physics_data(
			cur_case_output_file, &densities->fea_results, densities->dim_x, densities->dim_y, densities->cell_size, mesh.offset, "Vonmises"
		);
		if (physics_loaded) cout << "OptimizerBase: Finished reading stress distribution from .vtk file." << endl;
		else {
			cout << "OptimizerBase: Error: Unable to read physics data from file " << cur_case_output_file << endl;
		}
	}
};


