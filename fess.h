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
		phys::FEACaseManager fea_casemanager, msh::SurfaceMesh mesh, string base_folder, double _min_stress_threshold,
		grd::Densities2d _densities, int max_iterations, float _greediness, bool _do_feasibility_filtering,
		bool export_msh = false, bool verbose = true
	) : OptimizerBase(fea_casemanager, mesh, base_folder, _densities, max_iterations, export_msh, verbose)
	{
		min_stress_threshold = _min_stress_threshold;
		greediness = _greediness;
		do_feasibility_filtering = _do_feasibility_filtering;
	}
	double min_stress_threshold = 1.0;
	float greediness;
	double relative_area = INFINITY;
	bool do_feasibility_filtering = false;
	int no_cells_removed = 0;

	void run();
	void log_termination(string final_valid_iteration_folder, int final_valid_iteration);
	int handle_floating_pieces(
		msh::FEMesh2D* fe_mesh, int no_cells_to_remove, int no_cells_removed, bool recurse = true
	);
	void fill_design_domain();
	void export_stats(string iteration_name);
	void export_meta_parameters(vector<string>* _ = 0);
};

