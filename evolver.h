#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <map>
#include "helpers.h"
#include <functional>


class Evolver {
public:
	Evolver() = default;
	Evolver(
		int _pop_size, float _mutation_rate, function<bool(uint*, int)> _terminate,
		int _tournament_size, int _dim_x, int _dim_y, int _dim_z = 0
	) {
		dim_x = _dim_x;
		dim_y = _dim_y;
		dim_z = _dim_z;
		if (dim_z == 0) {
			dim_z = 1; domain_2d = true;
		}
		pop_size = _pop_size;
		mutation_rate = _mutation_rate;
		tournament_size = _tournament_size;
		terminate = _terminate;
	}
	void do_2d_crossover(uint* parent1, uint* parent2, uint* child1, uint* child2, int dim_x, int dim_y);
	void init_population();
	void do_evolution();
private:
	int dim_x = 1;
	int dim_y = 1;
	int dim_z = 1;
	bool domain_2d = false;
	int pop_size = 1;
	float mutation_rate = 0;
	int tournament_size = 1;
	uint* population = 0;
	function<bool(uint*, int)> terminate = 0;
};

