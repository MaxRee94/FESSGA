#pragma once


namespace fessga {
	class evo {
	public:

		class Individual2d : public grd::Densities2d {
		public:
			Individual2d() = default;
			Individual2d(int _dim_x, int _dim_y, Vector3d diagonal) : grd::Densities2d(_dim_x, _dim_y, diagonal) {};
			void update_phenotype();
			void remove_isolated_material();
			void do_feasibility_filtering();
			void do_ground_element_filtering();
		};
	};
}

