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
    Controller(Input _input, string _output_folder, string action, int _dim_x, int _dim_y, int _dim_z, int size)
    {
        input = _input; output_folder = _output_folder;
        dim_x = _dim_x, dim_y = _dim_y, dim_z = _dim_z;

        // Initialize mesh lists
        vector<MatrixXd> V_list;
        vector<MatrixXi> F_list;

        // Initialize gui
        GUI gui = GUI(V_list, F_list);

        // Load mesh
        MatrixXd V;
        MatrixXi F;
        if (input.type == "image" || input.type == "distribution2d" || input.type == "distribution3d") {
            Vector3d object_size(size, size, size);
            mesh = msh::SurfaceMesh(object_size);
            cout << "Created dummy mesh" << endl;
        }
        else {
            fessga::IO::read_mesh(input.path, V, F);
            gui.load_example(&V, &F);
            mesh = msh::SurfaceMesh(&V, &F);
        }

        // Create output folder
        IO::create_folder_if_not_exists(output_folder);

        // Load density distribution from a file or generate it on the fly
        init_densities();

        // Execute action
        if (action == "export_mesh") {
            densities2d.print();

            // Obtain a grid-based FE representation based on the chosen mesh
            msh::FEMesh2D fe_mesh;
            msh::generate_FE_mesh(mesh, densities2d, fe_mesh);

            // Export the FE mesh in .msh format
            msh::export_as_msh_file(&fe_mesh, output_folder);

            // Export the FE mesh in Elmer format (in .header, .nodes, .elments, .boundaries files)
            msh::export_as_elmer_files(&fe_mesh, output_folder);

            cout << "Finished exporting FE mesh. Continuing to post actions.." << endl;           
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
            string densities_file = densities2d.do_export(output_folder);
            cout << "Exported density distribution to " << densities_file << endl;
        }
        else if (action == "test") {
            run_tests();
        }
    };
    void init_densities();
    void create_parents(grd::Densities2d parent1, grd::Densities2d parent2);
    void run_tests();
    bool test_2d_crossover();
    bool test_full_evolution();
    bool run_fess();
private:
    vector<MatrixXd> V_list;
    vector<MatrixXi> F_list;
    msh::SurfaceMesh mesh;
    grd::Densities2d densities2d;
    grd::Densities3d densities3d;
    grd::Grid3D grid;
    GUI gui;
    string mesh_description = "";
    uint no_cells;
    float domain_size;
    int dim_x;
    int dim_y;
    int dim_z;
    float cell_size = 0;
    string output_folder;
    Vector3d offset;
    Input input;
};


void Controller::run_tests() {
    cout << "Running tests...\n";
    bool success = test_2d_crossover();
    cout << "Crossover test success? " << to_string(success) << endl;
}


void Controller::init_densities() {
    // Load distribution from file or generate it on the fly
    densities2d = grd::Densities2d(dim_x, dim_y, mesh.diagonal);
    if (input.type == "distribution2d") {
        cout << "Importing 2d density distribution from location " << input.path << "\n";
        densities2d.do_import(input.path, mesh.diagonal);
    }
    else if (input.type == "distribution3d") {
        cout << "Importing 3d density distribution from location " << input.path << "\n";
        densities3d.do_import(input.path, mesh.diagonal);
    }
    else if (input.type == "image") {
        cout << "Loading densities from image...\n";
        char img_path[200];
        strcpy(img_path, input.path.c_str());
        img::load_distribution_from_image(densities2d, img_path);
    }
    else {
        // Init 3d density distribution
        densities3d = grd::Densities3d(dim_x, dim_y, dim_z, mesh.diagonal);

        // Generate grid-based binary density distribution based on the given (unstructured) mesh file
        densities3d.generate(mesh.offset, mesh.V, mesh.F);

        // Create slice from 3d binary density distribution for 2d surface generation
        int z = dim_z / 2;
        densities3d.create_slice(densities2d, 2, z);
        densities2d.filter();

        // Export 3d density distribution
        string densities3d_file = densities3d.do_export(output_folder);
        cout << "FESS: Exported 3d density distribution to file: " << densities3d_file << endl;

        // Export 2d density distribution
        string densities_file = densities2d.do_export(output_folder);
        cout << "FESS: Exported 2d density distribution to file: " << densities_file << endl;

        // Write 2d distribution to image
        img::write_distribution_to_image(densities2d, output_folder + input.name + ".jpg", 1000, 1000);
    }
}

/*
Create 2 parent slices from the 3d binary density distribution for 2d test
*/
void Controller::create_parents(grd::Densities2d parent1, grd::Densities2d parent2) {
    int z1 = dim_x / 2;
    int z2 = dim_x / 2 + (dim_x / 5);
    densities3d.create_slice(parent1, 2, z1);
    densities3d.create_slice(parent2, 2, z2);
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

    string fe_case = "../data/msh_output/case.sif";
    string msh_file = "../data/msh_output/test.msh";
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
        msh_file, fe_case, mesh, output_folder, min_stress, max_stress, densities2d, max_iterations, greediness,
        maintain_boundary_cells, export_msh, verbose
    );
    fess.run();

    return true;
}
