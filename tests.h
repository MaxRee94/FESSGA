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
#include "fess.h"

using namespace fessga;

class Tester {
public:
    Tester()
    {
        // Initialize RNG
        srand(time(0));

        // Grid parameters
        mesher::Grid3D grid3d;
        grid3d.x = 5;
        grid3d.y = 5;
        grid3d.z = 5;
        float domain_size = 2.0; // TODO: replace

        // Initialize mesh lists
        vector<MatrixXd> V_list;
        vector<MatrixXi> F_list;

        // Initialize gui
        GUI gui = GUI(V_list, F_list);

        // Load example object
        MatrixXd V;
        MatrixXi F;
        fessga::IO::ReadMesh("../data/test_objects/teapot.obj", V, F);
        gui.load_example(&V, &F);
        mesher::SurfaceMesh surface_mesh = mesher::create_surface_mesh(&V, &F);

        // -- Normalize mesh --
        // Align barycenter to world origin
        // Get bounding box (min and max for x,y,z)
        // Get cell size along each dimension
        double cell_size = domain_size / (double)grid3d.x;
        // Get offset along each dimension
        Vector3d offset = -cell_size * 0.5 * Vector3d((double)grid3d.x, (double)grid3d.y, (double)grid3d.z);

        string output_folder = "E:/Development/FESSGA/data/msh_output/test/7_element_project";
        
        // Change domain to 2d
        dim_z = 0;
        no_cells = dim_x * dim_y;
    };
    void create_parents(uint* parent1, uint* parent2);
    bool test_2d_crossover();
    bool test_full_evolution();
    bool test_fess();
private:
    vector<MatrixXd> V_list;
    vector<MatrixXi> F_list;
    GUI gui;
    map<uint, uint> line_bounds;
    string mesh_description = "";
    uint no_cells;
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
    double max_stress = 1e9; // arbitrary maximum stress

    create_parents(parent1, parent2);
    cout << "\nParent 1: \n";
    mesher::print_2d_density_distrib(parent1, dim_x, dim_y);
    cout << "\nParent 2: \n";
    mesher::print_2d_density_distrib(parent2, dim_x, dim_y);

    string msh_file = "../data/msh_output/test.msh";
    string case_file = "../data/msh_output/case.sif";
    string output_folder = "../data/msh_output/FESSGA_test_output";

    // Do crossover
    uint* child1 = new uint[dim_x * dim_y];
    uint* child2 = new uint[dim_x * dim_y];
    Evolver evolver = Evolver(
        msh_file, case_file, output_folder, 4, (float)0.01, &variation_minimum_passed, 2, max_stress, parent1, dim_x, dim_y);
    evolver.do_2d_crossover(parent1, parent2, child1, child2);

    cout << "\Child 1: \n";
    mesher::print_2d_density_distrib(child1, dim_x, dim_y);
    cout << "\Child 2: \n";
    mesher::print_2d_density_distrib(child2, dim_x, dim_y);

    return true;
}

bool Tester::test_full_evolution() {
    Evolver evolver = Evolver();
    
    return true;
}

bool Tester::test_fess() {
    // Parameters
    double max_stress = 1e9;
    double min_stress = 1e-3;
    string msh_file = "../data/msh_output/test.msh";
    string case_file = "../data/msh_output/case.sif";
    string output_folder = "../data/msh_output/FESSGA_test_output";
    
    // Run optimization
    FESS fess = FESS(msh_file, case_file, output_folder, min_stress, max_stress, densities, dim_x, dim_y);
    fess.run();

    return true;
}
