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
			Individual2d(int _dim_x, int _dim_y, Vector3d diagonal) : grd::Densities2d(_dim_x, _dim_y, diagonal) {
			};
			int phenotype_count() {
				if (_phenotype_count == -1) redo_count();
				return _phenotype_count;
			}
			void update_phenotype();
			void remove_isolated_material();
			void do_feasibility_filtering();
			void do_ground_element_filtering();
		protected:
			uint* phenotype = 0;
			int _phenotype_count = -1;
		};
	};
}

