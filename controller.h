#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <map>
#include "io.h"
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
        fea_case = phys::FEACase(base_folder + "/case.sif", dim_x + 1, dim_y + 1, INFINITY);
        
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

        // Load density distribution from a file or generate it on the fly
        if (action != "evolve") init_densities();

        // Execute action
        if (action == "export_mesh") {
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
        else if (action == "evolve") {
            run_evoma();
        }
        else if (action == "fess") {
            // Run FESS
            run_fess();
        }
        else if (action == "export_distribution") {
            string densities_file = densities2d.do_export(base_folder + "/distribution2d.dens");
            cout << "Exported density distribution to " << densities_file << endl;
        }
        else if (action == "test") cout << "Test mode; controller remains passive." << endl;
        else {
            cerr << "Error: Action '" << action << "' not recognized.\n" << endl;
        }
    };
    void init_densities(phys::FEACase* fea_case = 0);
    void run_fess(FESS& _fess);
    void run_fess();
    void run_evoma(Evolver& _evolver);
    void run_evoma();

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
    phys::FEACase fea_case;
};

// Load distribution from file or generate it on the fly
void Controller::init_densities(phys::FEACase* _fea_case) {
    if (!_fea_case) _fea_case = &fea_case;
    densities2d = grd::Densities2d(dim_x, mesh.diagonal, base_folder);
    densities2d.fea_case = _fea_case;
    if (input.type == "distribution2d") {
        cout << "Importing 2d density distribution from location " << input.path << "\n";
        densities2d.do_import(input.path, mesh.diagonal(0));
        mesh = msh::SurfaceMesh(densities2d.diagonal);
    }
    else if (input.type == "distribution3d") {
        cout << "Importing 3d density distribution from location " << input.path << "\n";
        densities3d.do_import(input.path, mesh.diagonal(0));
        mesh = msh::SurfaceMesh(densities3d.diagonal);
    }
    else if (input.type == "image") {
        cout << "Loading densities from image...\n";
        char img_path[200];
        strcpy(img_path, input.path.c_str());
        img::load_distribution_from_image(densities2d, mesh, img_path);
    }
    else {
        // Init 3d density distribution
        densities3d = grd::Densities3d(dim_x, mesh.diagonal, base_folder);

        // Generate grid-based binary density distribution based on the given (unstructured) mesh file
        densities3d.generate(mesh.offset, mesh.V, mesh.F);

        // Create slice from 3d binary density distribution for 2d surface generation
        int z = dim_z / 2;
        densities3d.create_slice(densities2d, 2, z);
        densities2d.filter();

        // Export 3d density distribution
        string densities3d_file = densities3d.do_export(base_folder);
        cout << "FESS: Exported 3d density distribution to file: " << densities3d_file << endl;

        // Export 2d density distribution
        string densities2d_file = densities2d.do_export(base_folder);
        cout << "FESS: Exported 2d density distribution to file: " << densities2d_file << endl;

        // Write 2d distribution to image
        img::write_distribution_to_image(densities2d, base_folder + "/" + input.name + ".jpg");
    }
    _fea_case->dim_x = densities2d.dim_x + 1;
    _fea_case->dim_y = densities2d.dim_y + 1;
}


void Controller::run_fess(FESS& _fess) {
    // Parameters
    fea_case.max_stress_threshold = 1.5e9;
    fea_case.maintain_boundary_connection = true;
    double min_stress = 7e3;
    string msh_file = base_folder + "/mesh.msh";
    int max_iterations = 100;
    bool export_msh = true;
    float greediness = 0.05;
    bool verbose = true;
    phys::FEACaseInterpolator fea_interpolator;
    fea_interpolator.source = fea_case;
    
    // Run optimization
    FESS fess = FESS(
        fea_interpolator, mesh, base_folder, min_stress, densities2d, max_iterations, greediness,
        export_msh, verbose
    );
    _fess = fess;
    fess.run();
}

void Controller::run_evoma(Evolver& _evolver) {
    // Generic optimizer parameters
    string msh_file = base_folder + "/mesh.msh";
    int max_iterations = 100000;
    bool export_msh = true;
    bool verbose = true;

    // FEA Parameters
    string source_case_folder = base_folder + "/cases/coelophysis_bauri";
    string target_case_folder = base_folder + "/cases/trex_v02";

    // Load source- and target shape's densities and corresponding FEACase
    // -- Target
    input.path = target_case_folder + "/img.jpg";
    input.type = "image";
    phys::FEACase fea_case_target = phys::FEACase(target_case_folder + "/case.sif", densities2d.dim_x + 1, densities2d.dim_y + 1, INFINITY);
    init_densities(&fea_case_target);
    grd::Densities2d target_densities = densities2d;
    
    // -- Source
    input.path = source_case_folder + "/img.jpg";
    input.type = "image";
    phys::FEACase fea_case_source = phys::FEACase(source_case_folder + "/case.sif", densities2d.dim_x + 1, densities2d.dim_y + 1, INFINITY);
    init_densities(&fea_case_source);
    grd::Densities2d source_densities = densities2d;

    // Create interpolator from source and target
    phys::FEACaseInterpolator fea_interpolator(*source_densities.fea_case, *target_densities.fea_case);
    source_densities.fea_case = &fea_interpolator.source;
    target_densities.fea_case = &fea_interpolator.target;
    fea_interpolator.source.max_stress_threshold = 1e4;
    fea_interpolator.target.max_stress_threshold = 1.5e6;

    // Derive bound conditions for source and target
    msh::derive_boundary_conditions(*source_densities.fea_case, source_densities, mesh);
    msh::derive_boundary_conditions(*target_densities.fea_case, target_densities, mesh);

    // Initialize fea case interpolator
    fea_interpolator.initialize();
    
    // temp
    cout << "\nsource keep:\n";
    source_densities.visualize_keep_cells();
    cout << "\nsource cutouts:\n";
    source_densities.visualize_cutout_cells();
    
    cout << "\ntarget keep:\n";
    target_densities.visualize_keep_cells();
    cout << "\ntarget cutouts:\n";
    target_densities.visualize_cutout_cells();

    densities2d.fea_case = &fea_interpolator.interpolated;
    cout << "\interpolated keep:\n";
    densities2d.visualize_keep_cells();
    cout << "\interpolated cutouts:\n";
    densities2d.visualize_cutout_cells();


    // Evolver-specific parameters
    string crossover_method = "2x";
    float initial_perturb_level0 = 0.1;
    float initial_perturb_level1 = 0.2;
    int pop_size = 12; // NOTE: must be divisible by 4
    float mutation_rate_level0 = 0.0002;
    float mutation_rate_level1 = 0.001;
    int max_iterations_without_change = 150;
    float variation_trigger = 1.5;
    int no_static_iterations_trigger = 6;
    
    // Initialize and run evolver
    Evolver evolver(
        fea_interpolator, mesh, base_folder, pop_size, no_static_iterations_trigger, mutation_rate_level0, mutation_rate_level1, densities2d,
        variation_trigger, max_iterations, max_iterations_without_change, export_msh, verbose, initial_perturb_level0, initial_perturb_level1,
        crossover_method
    );
    _evolver = _evolver;
    evolver.evolve();
}

void Controller::run_evoma() {
    Evolver evolver;
    run_evoma(evolver);
}

void Controller::run_fess() {
    FESS fess;
    run_fess(fess);
}

