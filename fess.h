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
		string _msh_file, string _case_file, string _output_folder, double _min_stress_threshold, double _max_stress_threshold, 
		uint* _starting_densities, int _dim_x, int _dim_y, int _dim_z = 0
	) : OptimizerBase(_msh_file, _case_file, _output_folder, _max_stress_threshold, _starting_densities, _dim_x, _dim_y, _dim_z)
	{
		min_stress_threshold = _min_stress_threshold;
	}
	double min_stress_threshold = 1.0;
	void run();
};


void FESS::run() {
	cout << "Beginning FESS run. Saving results to " << output_folder << endl;

	double min_stress = 1.0;
	string msh = msh_file;
	string cur_iteration_name = "iteration_0001";
	string cur_output_folder = output_folder + "/" + cur_iteration_name;
	int i = 0;
	while (min_stress > min_stress_threshold) {
		cout << "--- Starting iteration " << i << ". Previous lowest stress value: " << min_stress << " N/m^2" << endl;

		// Call Elmer to run FEA on .msh file in cur_output_folder, using .sif file
		//call_elmer(output_folder + "/" + cur_iteration_name);

		// Wait for Elmer's analysis to complete. This is the case when a new .vtk file has appeared

		// Obtain stress distribution from the .vtk file
		min_stress = 0.0;

		// Set all cells that have corresponding stress values lower than <min_stress_threshold> to density=0.

		// Update current output folder and copy Elmer's project files to it
		string cur_iteration_name = fessga::help::add_padding("iteration_", i + 1) + to_string(i + 1);
		string cur_output_folder = output_folder + "/" + cur_iteration_name;
		
		// Generate new .msh file using modified density distribution

		i++;
	}

	cout << "Finished FESS run. Results were saved to " << output_folder << endl;
}
