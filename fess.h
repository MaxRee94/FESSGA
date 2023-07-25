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
		phys::FEACaseManager _fea_manager, msh::SurfaceMesh _mesh, string _base_folder, double _min_stress_threshold,
		grd::Densities2d _densities, int _max_iterations, float _greediness,
		bool _export_msh = false, bool _verbose = true
	) : OptimizerBase(_fea_manager, _mesh, _base_folder, _densities, _max_iterations, _export_msh, _verbose)
	{
		min_stress_threshold = _min_stress_threshold;
		greediness = _greediness;
	}
	double min_stress_threshold = 1.0;
	float greediness;

	void run();
	void log_termination(string final_valid_iteration_folder, int final_valid_iteration);
	int handle_floating_pieces(
		msh::FEMesh2D* fe_mesh, int no_cells_to_remove, int no_cells_removed, bool recurse = true
	);
	void export_stats(string iteration_name);
};

