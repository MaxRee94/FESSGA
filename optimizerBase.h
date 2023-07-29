#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <map>
#include "images.h"

using namespace fessga;


void load_physics(grd::Densities2d* densities, msh::SurfaceMesh* mesh, bool verbose = false);

class OptimizerBase {
public:
	OptimizerBase() = default;

	// Constructor for 2d optimization
	OptimizerBase(
		phys::FEACaseManager _fea_casemanager, msh::SurfaceMesh _mesh, string _base_folder,
		grd::Densities2d _densities, int _max_iterations, bool _export_msh, bool _verbose
	) {
		mesh = _mesh;
		base_folder = IO::get_fullpath(_base_folder);
		output_folder = IO::get_unique_path(base_folder + "/run_#");
		densities = _densities;
		export_msh = _export_msh;
		domain_2d = true;
		no_cells = densities.size;
		max_iterations = _max_iterations;
		verbose = _verbose;
		IO::create_folder_if_not_exists(output_folder);
		image_folder = IO::create_folder_if_not_exists(output_folder + "/image_output");
		fea_casemanager = _fea_casemanager;
		densities.fea_casemanager = &fea_casemanager;
	};
	msh::SurfaceMesh mesh;
	bool domain_2d = false;
	int no_cells = 1;
	grd::Densities2d densities;
	bool verbose = true;
	int max_iterations = 0;
	int iteration_number = 0;
	bool export_msh = false;
	phys::FEACaseManager fea_casemanager;
	double min_stress, max_stress;
	string base_folder, image_folder, iteration_name, iteration_folder, output_folder;

	// Function to get the folder corresponding to the given iteration number. If the folder does not exist yet, it will be created.
	string get_iteration_folder(int iteration, bool verbose = false) {
		iteration_name = fessga::help::add_padding("iteration_", iteration) + to_string(iteration);
		iteration_folder = output_folder + "/" + iteration_name;
		if (verbose) cout << "OptimizerBase: Creating folder " << iteration_folder <<
			" for current iteration if it does not exist yet...\n";
		IO::create_folder_if_not_exists(iteration_folder);

		return iteration_folder;
	}

	// Write densities to image file
	virtual void write_densities_to_image(bool verbose = false) {
		img::write_distribution_to_image(densities, image_folder + "/" + iteration_name + ".jpg", true);
	}

	virtual void export_meta_parameters(vector<string>* additional_metaparameters) {
		vector<string> _content = {
			"max stress threshold = " + to_string(fea_casemanager.max_stress_threshold),
		};
		help::append_vector(_content, additional_metaparameters);
		string content = help::join(&_content, "\n");
		IO::write_text_to_file(content, output_folder + "/metaparameters.txt");
	}

	// Function to write statistics, images, and other generated data to files
	void export_base_stats() {
		write_densities_to_image();
	}

	// Copy all files in a solution folder (iteration folder for FESS, individual folder for emma) to the specified target folder
	void copy_solution_files(string source_dir, string target_dir, bool verbose = false) {
		if (verbose) cout << "Copying solution files from directory " << source_dir << " to " << target_dir << endl;
		vector<string> filepaths;
		IO::get_files_in_directory(filepaths, source_dir);
		for (auto& fpath : filepaths) {
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

};


