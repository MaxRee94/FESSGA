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

class Controller {
public:
    Controller()
    {
        // Initialize RNG
        srand(time(0));

        // Initialize mesh lists
        vector<MatrixXd> V_list;
        vector<MatrixXi> F_list;

        // Initialize gui
        GUI gui = GUI(V_list, F_list);

        // Load example object
        MatrixXd V;
        MatrixXi F;
        fessga::IO::read_mesh("../data/test_objects/teapot.obj", V, F);
        gui.load_example(&V, &F);
        mesher::SurfaceMesh surface_mesh = mesher::create_surface_mesh(&V, &F);

        // Create 3d Grid
        grid = mesher::create_grid3d(40, 40, 40, surface_mesh.diagonal);

        // Set output folder
        string output_folder = "E:/Development/FESSGA/data/msh_output/test";
        
        // Compute no of cells
        no_cells = dim_x * dim_y;

        // Generate grid-based binary density distribution based on the given (unstructured) mesh file
        uint32_t* densities3d = new uint32_t[grid.x * grid.y * grid.z];
        mesher::generate_3d_density_distribution(grid, surface_mesh, &gui.V_list[0], &gui.F_list[0], densities3d);

        // Create slice from 3d binary density distribution for 2d test
        int z = grid.x / 2;
        densities = new uint[grid.x * grid.y];
        mesher::create_2d_slice(densities3d, densities, grid, z);
        mesher::filter_2d_density_distrib(densities, grid.x, grid.y);
        mesher::print_density_distrib(densities, grid.x, grid.y);
    };
    void create_parents(uint* parent1, uint* parent2);
    bool test_2d_crossover();
    bool test_full_evolution();
    bool run_fess();
private:
    vector<MatrixXd> V_list;
    vector<MatrixXi> F_list;
    mesher::SurfaceMesh mesh;
    mesher::Grid3D grid;
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
void Controller::create_parents(uint* parent1, uint* parent2) {
    
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
bool Controller::test_2d_crossover() {
    uint* parent1 = new uint[dim_x * dim_y];
    uint* parent2 = new uint[dim_x * dim_y];
    double max_stress = 1e9; // arbitrary maximum stress
    int max_iterations = 100;

    create_parents(parent1, parent2);
    cout << "\nParent 1: \n";
    mesher::print_density_distrib(parent1, dim_x, dim_y);
    cout << "\nParent 2: \n";
    mesher::print_density_distrib(parent2, dim_x, dim_y);

    string msh_file = "../data/msh_output/test.msh";
    string case_file = "../data/msh_output/case.sif";
    string output_folder = "../data/msh_output/FESSGA_test_output";

    // Do crossover
    uint* child1 = new uint[dim_x * dim_y];
    uint* child2 = new uint[dim_x * dim_y];
    Evolver evolver = Evolver(
        msh_file, case_file, mesh, output_folder, 4, (float)0.01, &variation_minimum_passed, 2, max_stress, parent1, grid, max_iterations
    );
    evolver.do_2d_crossover(parent1, parent2, child1, child2);

    cout << "\Child 1: \n";
    mesher::print_density_distrib(child1, dim_x, dim_y);
    cout << "\Child 2: \n";
    mesher::print_density_distrib(child2, dim_x, dim_y);

    return true;
}

bool Controller::test_full_evolution() {
    Evolver evolver = Evolver();
    
    return true;
}

bool Controller::run_fess() {
    // Parameters
    double max_stress = 1e10;
    double min_stress = 7e3;
    string msh_file = "../data/msh_output/test.msh";
    string output_folder = "../data/msh_output/FESSGA_test_output_40elements";
    string case_file = output_folder + "/case.sif";
    int max_iterations = 100;
    float greediness = 0.02;
    
    // Run optimization
    FESS fess = FESS(msh_file, case_file, mesh, output_folder, min_stress, max_stress, densities, grid, max_iterations, greediness);
    fess.run();

    return true;
}