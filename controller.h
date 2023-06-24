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

struct Input {
    string type = "object";
    string path;
    string name;
};

class Controller {
public:
    Controller(Input _input, string _output_folder, string action, int dim_x, int dim_y, int dim_z, int size)
    {
        input = _input; output_folder = _output_folder;

        // Initialize mesh lists
        vector<MatrixXd> V_list;
        vector<MatrixXi> F_list;

        // Initialize gui
        GUI gui = GUI(V_list, F_list);

        // Load mesh
        MatrixXd V;
        MatrixXi F;
        if (input.type == "image" || input.type == "distribution") {
            Vector3d object_size(size, size, size);
            mesh = mesher::SurfaceMesh(object_size);
            cout << "Created dummy mesh" << endl;
        }
        else {
            fessga::IO::read_mesh(input.path, V, F);
            gui.load_example(&V, &F);
            mesh = mesher::SurfaceMesh(&V, &F);
        }

        // Create 3d grid
        grid = mesher::create_grid3d(dim_x, dim_y, dim_z, mesh.diagonal);

        // Create output folder
        IO::create_folder_if_not_exists(output_folder);

        // Load density distribution or generate it on the fly
        load_density_distribution();

        // Execute action
        if (action == "export_mesh") {
            mesher::print_density_distrib(densities, grid.x, grid.y);

            // Obtain a grid-based FE representation based on the chosen mesh
            mesher::FEMesh2D fe_mesh;
            mesher::generate_FE_mesh(grid, mesh, densities, fe_mesh);

            // Export the FE mesh in .msh format
            mesher::export_as_msh_file(&fe_mesh, output_folder);

            // Export the FE mesh in Elmer format (in .header, .nodes, .elments, .boundaries files)
            mesher::export_as_elmer_files(&fe_mesh, output_folder);

            cout << "Finished exporting FE mesh. Continuing to post actions.." << endl;

            // -- Post actions
            // Write distribution to image
            img::write_distribution_to_image(densities, grid, output_folder + input.name + ".jpg", 1000, 1000);
        }
        else if (action == "evolve") {
            // Do crossover test
            test_2d_crossover();
        }
        else if (action == "fess") {
            // Do FESS
            run_fess();
        }
        else if (action == "export_distribution") {
            string densities_file = mesher::export_density_distrib(output_folder, densities, grid.x, grid.y);
            cout << "Exported density distribution to " << densities_file << endl;
        }
        else if (action == "test") {
            run_tests();
        }
    };
    void load_density_distribution();
    void create_parents(uint* parent1, uint* parent2);
    void run_tests();
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
    Input input;
};


void Controller::run_tests() {

}


void Controller::load_density_distribution() {
    // Load distribution from file or generate it on the fly
    densities = new uint[grid.x * grid.y];
    if (input.type == "distribution") {
        densities = mesher::import_densities(input.path, grid, mesh.diagonal);
        mesher::print_density_distrib(densities, grid.x, grid.y);
    }
    else if (input.type == "image") {
        cout << "Loading densities from image...\n";
        char img_path[200];
        strcpy(img_path, input.path.c_str());
        img::load_distribution_from_image(img_path, densities, grid);
    }
    else {
        // Generate grid-based binary density distribution based on the given (unstructured) mesh file
        uint* densities3d = new uint[grid.x * grid.y * grid.z];
        mesher::generate_3d_density_distribution(grid, mesh, &gui.V_list[0], &gui.F_list[0], densities3d);

        // Create slice from 3d binary density distribution for 2d surface generation
        int z = grid.x / 2;
        mesher::create_2d_slice(densities3d, densities, grid, z);
        mesher::filter_2d_density_distrib(densities, grid.x, grid.y);

        // Export density distribution
        string densities_file = mesher::export_density_distrib(output_folder, densities, grid.x, grid.y);
        cout << "FESS: Exported current density distribution to file: " << densities_file << endl;
    }
}

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
    double max_stress = 1.5e9;
    double min_stress = 7e3;
    string msh_file = output_folder + "/mesh.msh";
    string fe_case = output_folder + "/case.sif";
    int max_iterations = 100;
    bool export_msh = true;
    float greediness = 0.1;
    bool verbose = true;
    bool maintain_boundary_cells = true;
    
    // Run optimization
    FESS fess = FESS(
        msh_file, fe_case, mesh, output_folder, min_stress, max_stress, densities, grid, max_iterations, greediness,
        maintain_boundary_cells, export_msh, verbose
    );
    fess.run();

    return true;
}
