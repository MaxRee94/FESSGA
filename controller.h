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
        GUI gui(V_list, F_list);

        // Load example object
        MatrixXd V;
        MatrixXi F;
        fessga::IO::read_mesh("../data/objects/teapot.obj", V, F);
        gui.load_example(&V, &F);
        mesh = mesher::SurfaceMesh(&V, &F);

        // Create 3d Grid
        grid = mesher::create_grid3d(200, 200, 200, mesh.diagonal);

        // Set output folder
        output_folder = "E:/Development/FESSGA/data/FESS/FESS_teapot_40_elements";
        
        // Compute no of cells
        no_cells = grid.x * grid.y;
        densities = new uint[no_cells];
#if 0:
        // Generate grid-based binary density distribution based on the given (unstructured) mesh file
        uint* densities3d = new uint[grid.x * grid.y * grid.z];
        mesher::generate_3d_density_distribution(grid, surface_mesh, &gui.V_list[0], &gui.F_list[0], densities3d);

        // Create slice from 3d binary density distribution for 2d test
        int z = grid.x / 2;
        mesher::create_2d_slice(densities3d, densities, grid, z);
        mesher::filter_2d_density_distrib(densities, grid.x, grid.y);
        mesher::print_density_distrib(densities, grid.x, grid.y);
#else
        string densities_path = output_folder + "/distribution.dens";
        mesher::import_densities(densities_path, densities);
#endif
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
    string output_folder;
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
    bool export_msh = false;

    create_parents(parent1, parent2);
    cout << "\nParent 1: \n";
    mesher::print_density_distrib(parent1, dim_x, dim_y);
    cout << "\nParent 2: \n";
    mesher::print_density_distrib(parent2, dim_x, dim_y);

    string msh_file = "../data/msh_output/test.msh";
    string fe_case = "../data/msh_output/case.sif";
    string output_folder = "../data/msh_output/FESSGA_test_output";

    // Do crossover
    uint* child1 = new uint[dim_x * dim_y];
    uint* child2 = new uint[dim_x * dim_y];
    Evolver evolver = Evolver(
        msh_file, fe_case, mesh, output_folder, 4, (float)0.01, &variation_minimum_passed, 2, max_stress, parent1, grid, max_iterations
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
    double max_stress = 6e4;
    double min_stress = 7e3;
    string msh_file = output_folder + "/mesh.msh";
    string fe_case = output_folder + "/case.sif";
    int max_iterations = 100;
    bool export_msh = true;
    float greediness = 0.2;
    bool verbose = true;
    bool maintain_boundary_cells = false;
    
    // Run optimization
    FESS fess = FESS(
        msh_file, fe_case, mesh, output_folder, min_stress, max_stress, densities, grid, max_iterations, greediness,
        maintain_boundary_cells, export_msh, verbose
    );
    fess.run();

    return true;
}
