#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <map>
#include "helpers.h"
#include "densities.h"


namespace fessga {
	class evo {
	public:

		class Individual2d : public grd::Densities2d {
		public:
			Individual2d() = default;
			Individual2d(
				int _dim_x, int _dim_y, Vector3d diagonal, phys::FEAResults2D* _fea_results,
				phys::FEACase* _fea_case) : grd::Densities2d(_dim_x, _dim_y, diagonal, _fea_results, _fea_case
			) {};
			Individual2d(grd::Densities2d densities, Vector3d diagonal) : grd::Densities2d(
				densities.dim_x, densities.dim_y, diagonal, densities.fea_results, densities.fea_case
			) {
				copy_from(&densities);
			};
			Individual2d(Individual2d* individual) {
				construct_grid(individual->dim_x, individual->dim_y);
				cell_size = individual->cell_size;
				fea_results = individual->fea_results;
				fea_case = individual->fea_case;
				copy_from_individual(individual);
			}
			int phenotype_count() {
				if (_phenotype_count == -1) redo_count();
				return _phenotype_count;
			}
			void update_phenotype();
			bool repair();
			bool remove_isolated_material();
			void do_single_feasibility_filtering_pass();
			void do_feasibility_filtering(bool verbose = false);
			void do_ground_element_filtering();
			void fill_voids(int no_true_neighbors = 4);
			void copy_from_individual(Individual2d* source);
		protected:
			uint* phenotype = 0;
			int _phenotype_count = -1;
		};
	};
}

