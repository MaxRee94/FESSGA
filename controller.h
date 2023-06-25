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
#include "evolver.h"
#include "fess.h"


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
        grid = grd::Grid3D(dim_x, dim_y, dim_z, mesh.diagonal);

        // Create output folder
        IO::create_folder_if_not_exists(output_folder);

        // Load density distribution or generate it on the fly
        load_density_distribution();

        // Execute action
        if (action == "export_mesh") {
            densities.print();

            // Obtain a grid-based FE representation based on the chosen mesh
            mesher::FEMesh2D fe_mesh;
            mesher::generate_FE_mesh(mesh, densities, fe_mesh);

            // Export the FE mesh in .msh format
            mesher::export_as_msh_file(&fe_mesh, output_folder);

            // Export the FE mesh in Elmer format (in .header, .nodes, .elments, .boundaries files)
            mesher::export_as_elmer_files(&fe_mesh, output_folder);

            cout << "Finished exporting FE mesh. Continuing to post actions.." << endl;

            // -- Post actions
            // Write distribution to image
            img::write_distribution_to_image(densities, output_folder + input.name + ".jpg", 1000, 1000);
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
            string densities_file = densities.do_export(output_folder);
            cout << "Exported density distribution to " << densities_file << endl;
        }
        else if (action == "test") {
            run_tests();
        }
    };
    void load_density_distribution();
    void create_parents(grd::Densities2d parent1, grd::Densities2d parent2);
    void run_tests();
    bool test_2d_crossover();
    bool test_full_evolution();
    bool run_fess();
private:
    vector<MatrixXd> V_list;
    vector<MatrixXi> F_list;
    mesher::SurfaceMesh mesh;
    grd::Grid3D grid;
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
    grd::Densities2d densities;
    Input input;
};


void Controller::run_tests() {
    cout << "Running tests...\n";
}


void Controller::load_density_distribution() {
    // Load distribution from file or generate it on the fly
    densities = grd::Densities2d(grid.x, grid.y, mesh.diagonal);
    if (input.type == "distribution") {
        densities.do_import(input.path, mesh.diagonal);
    }
    else if (input.type == "image") {
        cout << "Loading densities from image...\n";
        char img_path[200];
        strcpy(img_path, input.path.c_str());
        img::load_distribution_from_image(densities, img_path);
    }
    else {
        // Generate grid-based binary density distribution based on the given (unstructured) mesh file
        uint* densities3d = new uint[grid.x * grid.y * grid.z];
        mesher::generate_3d_density_distribution(grid, mesh, &gui.V_list[0], &gui.F_list[0], densities3d);

        // Create slice from 3d binary density distribution for 2d surface generation
        int z = grid.x / 2;
        mesher::create_2d_slice(densities3d, densities, grid, z);
        densities.filter();

        // Export density distribution
        string densities_file = densities.do_export(output_folder);
        cout << "FESS: Exported 2d slice density distribution to file: " << densities_file << endl;
    }
}

/*
Create 2 parent slices from the 3d binary density distribution for 2d test
*/
void Controller::create_parents(grd::Densities2d parent1, grd::Densities2d parent2) {
    int z1 = dim_x / 2;
    int z2 = dim_x / 2 + (dim_x / 5);
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            parent1.set(x * dim_x + y, densities[z1 * dim_x * dim_y + x * dim_x + y]);
            parent2.set(x * dim_x + y, densities[z1 * dim_x * dim_y + x * dim_x + y]);
        }
    }
    parent1.update_count();
    parent2.update_count();
}

/*
Test 2-point crossover of two 2d parent solutions. Print parents and children to console
*/
bool Controller::test_2d_crossover() {
    grd::Densities2d parent1(dim_x, dim_y, mesh.diagonal);
    grd::Densities2d parent2(dim_x, dim_y, mesh.diagonal);
    double max_stress = 1e9; // arbitrary maximum stress
    int max_iterations = 100;
    bool export_msh = false;

    create_parents(parent1, parent2);
    cout << "\nParent 1: \n";
    parent1.print();
    cout << "\nParent 2: \n";
    parent2.print();

    string msh_file = "../data/msh_output/test.msh";
    string fe_case = "../data/msh_output/case.sif";
    string output_folder = "../data/msh_output/FESSGA_test_output";

    // Do crossover
    grd::Densities2d child1(dim_x, dim_y, mesh.diagonal);
    grd::Densities2d child2(dim_x, dim_y, mesh.diagonal);
    Evolver evolver = Evolver(
        msh_file, fe_case, mesh, output_folder, 4, (float)0.01, &variation_minimum_passed, 2,
        max_stress, parent1, max_iterations
    );
    evolver.do_2d_crossover(parent1, parent2, child1, child2);

    cout << "\Child 1: \n";
    child1.print();
    cout << "\Child 2: \n";
    child2.print();

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
        msh_file, fe_case, mesh, output_folder, min_stress, max_stress, densities, max_iterations, greediness,
        maintain_boundary_cells, export_msh, verbose
    );
    fess.run();

    return true;
}
