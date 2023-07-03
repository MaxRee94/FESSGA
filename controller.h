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
        // Initialize member variables
        input = _input; output_folder = _output_folder;
        dim_x = _dim_x, dim_y = _dim_y, dim_z = _dim_z;
        vector<MatrixXd> V_list;
        vector<MatrixXi> F_list;
        GUI gui = GUI(V_list, F_list);
        
        // Initialize RNG
        help::init_RNG();

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
            // Run evolutionary algorithm
        }
        else if (action == "fess") {
            // Run FESS
            run_fess();
        }
        else if (action == "export_distribution") {
            string densities_file = densities2d.do_export(output_folder);
            cout << "Exported density distribution to " << densities_file << endl;
        }
        else if (action == "test") cout << "Test mode; controller remains passive." << endl;
        else {
            cerr << "Error: Action '" << action << "' not recognized.\n" << endl;
        }
    };
    void init_densities();
    bool run_fess();

    vector<MatrixXd> V_list;
    vector<MatrixXi> F_list;
    msh::SurfaceMesh mesh;
    grd::Densities2d densities2d;
    grd::Densities3d densities3d;
    phys::FEAResults2D fea_results;
    phys::FEACase fea_case;
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
    vector<grd::Piece> pieces;
};

// Load distribution from file or generate it on the fly
void Controller::init_densities() {
    fea_results = phys::FEAResults2D(dim_x, dim_y);
    fea_case = phys::FEACase(output_folder + "/case.sif");
    densities2d = grd::Densities2d(dim_x, dim_y, mesh.diagonal, &fea_results, &fea_case);
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
        densities3d = grd::Densities3d(dim_x, dim_y, dim_z, mesh.diagonal, &fea_results, &fea_case);

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
        string densities2d_file = densities2d.do_export(output_folder);
        cout << "FESS: Exported 2d density distribution to file: " << densities2d_file << endl;

        // Write 2d distribution to image
        img::write_distribution_to_image(densities2d, output_folder + "/" + input.name + ".jpg", 1000, 1000);
    }
}


bool Controller::run_fess() {
    // Parameters
    double max_stress = 1.5e9;
    double min_stress = 7e3;
    string msh_file = output_folder + "/mesh.msh";
    string fea_case = output_folder + "/case.sif";
    int max_iterations = 100;
    bool export_msh = true;
    float greediness = 0.05;
    bool verbose = true;
    bool maintain_boundary_connection = true;
    
    
    // Run optimization
    FESS fess = FESS(
        msh_file, fea_case, mesh, output_folder, min_stress, max_stress, densities2d, max_iterations, greediness,
        maintain_boundary_connection, export_msh, verbose
    );
    fess.run();

    return true;
}


