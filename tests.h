#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <map>
#include "io.h"
#include "gui.h"
#include "meshing.h"
#include "helpers.h"
#include "evolver.h"


class Tester {
public:
    Tester()
    {
        // Initialize RNG
        srand(time(0));

        // Initialize gui
        gui = GUI(V_list, F_list);

        // Load example object
        gui.load_example();

        // Generate a 3d grid-based binary density distribution of the chosen mesh
        domain_size = 2.0;
        dim_x = 30;
        dim_y = 30;
        dim_z = 30;
        cell_size = domain_size / (float)dim_x;
        Vector3d offset = -cell_size * 0.5 * Vector3d((double)dim_x, (double)dim_y, (double)dim_z);
        densities = new uint[dim_x * dim_y * dim_z];
        generate_3d_density_distribution(dim_x, dim_y, dim_z, offset, cell_size, &gui.V_list[0], &gui.F_list[0], densities);
    };
    void create_parents(uint* parent1, uint* parent2);
    bool test_2d_crossover();
    bool test_full_evolution();
private:
    vector<MatrixXd> V_list;
    vector<MatrixXi> F_list;
    GUI gui;
    map<uint, uint> line_bounds;
    string mesh_description = "";
    float domain_size;
    int dim_x;
    int dim_y;
    int dim_z;
    float cell_size = 0;
    Vector3d offset;
    uint* densities = 0;
};

/*
Create 2 parent slices from the 3d binary density distribution for 2d test
*/
void Tester::create_parents(uint* parent1, uint* parent2) {
    
    int z1 = dim_x / 2;
    int z2 = dim_x / 2 + (dim_x / 5);
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            parent1[x * dim_x + y] = densities[z1 * dim_x * dim_y + x * dim_x + y];
            parent2[x * dim_x + y] = densities[z2 * dim_x * dim_y + x * dim_x + y];
        }
    }
}

/*
Test 2-point crossover of two 2d parent solutions. Print parents and children to console
*/
bool Tester::test_2d_crossover() {
    uint* parent1 = new uint[dim_x * dim_y];
    uint* parent2 = new uint[dim_x * dim_y];
    create_parents(parent1, parent2);
    cout << "\nParent 1: \n";
    print_2d_density_distrib(parent1, dim_x);
    cout << "\nParent 2: \n";
    print_2d_density_distrib(parent2, dim_x);

    // Do crossover
    uint* child1 = new uint[dim_x * dim_y];
    uint* child2 = new uint[dim_x * dim_y];
    Evolver evolver = Evolver();
    evolver.do_2d_crossover(parent1, parent2, child1, child2, dim_x, dim_y);

    cout << "\Child 1: \n";
    print_2d_density_distrib(child1, dim_x);
    cout << "\Child 2: \n";
    print_2d_density_distrib(child2, dim_x);

    return true;
}

bool Tester::test_full_evolution() {
    return true;
}
