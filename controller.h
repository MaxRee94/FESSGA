#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <map>
#include "gui.h"
#include "evolver.h"
#include "fess.h"


struct Input {
    string type = "object";
    string path;
    string name;
};

class Controller {
public:
    Controller(Input _input, string _base_folder, string action, int _dim_x, int size)
    {
        // Initialize member variables
        input = _input; base_folder = _base_folder;
        dim_x = _dim_x;
        vector<MatrixXd> V_list;
        vector<MatrixXi> F_list;
        GUI gui = GUI(V_list, F_list);
        fea_casemanager = phys::FEACaseManager();
        
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
        IO::create_folder_if_not_exists(base_folder);

        // Execute action
        if (action == "export_mesh") {
            init_densities();
            densities2d.print();
            densities2d.do_export(base_folder + "/distribution2d.dens");

            // Obtain a grid-based FE representation based on the chosen mesh
            msh::FEMesh2D fe_mesh;
            msh::create_FE_mesh(mesh, densities2d, fe_mesh);

            // Export the FE mesh in .msh format
            msh::export_as_msh_file(&fe_mesh, base_folder);

            // Export the FE mesh in Elmer format (in .header, .nodes, .elments, .boundaries files)
            msh::export_as_elmer_files(&fe_mesh, base_folder);

            // TODO: Update the case.sif file to fill in missing target boundary edges

            cout << "Finished exporting FE mesh. Continuing to post actions.." << endl;
        }
        else if (action == "evolve_dynamic") {
            run_emma_dynamic();
        }
        else if (action == "evolve_static") {
            run_emma_static();
        }
        else if (action == "fess") {
            run_fess();
        }
        else if (action == "export_distribution") {
            init_densities();
            string densities_file = densities2d.do_export(base_folder + "/distribution2d.dens");
            cout << "Exported density distribution to " << densities_file << endl;
        }
        else if (action == "test") cout << "Test mode; controller remains passive." << endl;
        else {
            cerr << "Error: Action '" << action << "' not recognized.\n" << endl;
        }
    };
    void init_densities(phys::FEACaseManager* fea_casemanager = 0);
    void run_fess(FESS& _fess);
    void run_fess();
    void run_emma_dynamic(Evolver& _evolver);
    void run_emma_dynamic();
    void run_emma_static(Evolver& _evolver);
    void run_emma_static();
    void run_emma(Evolver& _evolver, phys::FEACaseManager* fea_manager);

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
    int dim_x = 1;
    int dim_y = 1;
    int dim_z = 1;
    float cell_size = 0;
    string base_folder;
    Vector3d offset;
    Input input;
    vector<grd::Piece> pieces;
    phys::FEACaseManager fea_casemanager;
};
