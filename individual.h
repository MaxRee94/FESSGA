#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <map>
#include "meshing.h"


namespace fessga {
	class evo {
	public:

		class Individual2d : public grd::Densities2d {
		public:
			Individual2d() = default;
			Individual2d(
				int _dim_x, Vector3d _diagonal, string _output_folder) : grd::Densities2d(_dim_x, _diagonal, _output_folder)
			{};
			Individual2d(grd::Densities2d* densities) : grd::Densities2d(densities) {
				copy_from(densities);
			};
			Individual2d(Individual2d* individual) {
				construct_grid();
				cell_size = individual->cell_size;
				fea_results = individual->fea_results;
				fea_casemanager = individual->fea_casemanager;
				copy_from_individual(individual);
			}
			int phenotype_count() {
				if (_phenotype_count == -1) redo_count();
				return _phenotype_count;
			}
			void update_phenotype();
			void do_ground_element_filtering();
			void copy_from_individual(Individual2d* source);
			double fitness = 0;
			msh::FEMesh2D fe_mesh;
		protected:
			uint* phenotype = 0;
			int _phenotype_count = -1;
		};
	};
}

