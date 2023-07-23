#pragma once
#include "controller.h"


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
    init_densities();

    // Parameters
    fea_case.max_stress_threshold = 1.5e9;
    fea_case.maintain_boundary_connection = true;
    double min_stress = 7e3;
    string msh_file = base_folder + "/mesh.msh";
    int max_iterations = 100;
    bool export_msh = true;
    float greediness = 0.05;
    bool verbose = true;
    phys::FEACaseManager fea_manager;
    fea_manager.current = fea_case;
    densities2d.fea_case = &fea_manager.current;
    msh::derive_boundary_conditions(fea_case, densities2d, mesh);

    // Run optimization
    FESS fess = FESS(
        fea_manager, mesh, base_folder, min_stress, densities2d, max_iterations, greediness,
        export_msh, verbose
    );
    _fess = fess;
    fess.run();
}

void Controller::run_emma_static(Evolver& _evolver) {
    init_densities();
    phys::FEACaseManager fea_manager;
    fea_manager.current = fea_case;
    densities2d.fea_case = &fea_manager.current;
    msh::derive_boundary_conditions(*densities2d.fea_case, densities2d, mesh);
    run_emma(_evolver, &fea_manager);
}

void Controller::run_emma_dynamic(Evolver& _evolver) {
    // FEA Parameters
    string source_case_folder = base_folder + "/cases/coelophysis_bauri";
    string target_case_folder = base_folder + "/cases/trex_v02";

    // Load source- and target shape's densities and corresponding FEACase
    // -- Target
    input.path = target_case_folder + "/img.jpg";
    input.type = "image";
    phys::FEACase fea_case_target = phys::FEACase(
        target_case_folder + "/case.sif", densities2d.dim_x + 1, densities2d.dim_y + 1, INFINITY
    );
    init_densities(&fea_case_target);
    grd::Densities2d target_densities = densities2d;

    // -- Source
    input.path = source_case_folder + "/img.jpg";
    input.type = "image";
    phys::FEACase fea_case_source = phys::FEACase(
        source_case_folder + "/case.sif", densities2d.dim_x + 1, densities2d.dim_y + 1, INFINITY
    );
    init_densities(&fea_case_source);
    grd::Densities2d source_densities = densities2d;

    // Create interpolator from source and target
    phys::FEACaseManager fea_manager(*source_densities.fea_case, *target_densities.fea_case);
    source_densities.fea_case = &fea_manager.source;
    target_densities.fea_case = &fea_manager.target;
    fea_manager.source.max_stress_threshold = 1e4;
    fea_manager.target.max_stress_threshold = 1.5e6;

    // Derive bound conditions for source and target
    msh::derive_boundary_conditions(*source_densities.fea_case, source_densities, mesh);
    msh::derive_boundary_conditions(*target_densities.fea_case, target_densities, mesh);

    // Initialize fea case interpolator
    fea_manager.initialize();
    fea_manager.interpolate(0.5);

    // temp
    /*cout << "\nsource keep:\n";
    source_densities.visualize_keep_cells();
    cout << "\nsource cutouts:\n";
    source_densities.visualize_cutout_cells();

    cout << "\ntarget keep:\n";
    target_densities.visualize_keep_cells();
    cout << "\ntarget cutouts:\n";
    target_densities.visualize_cutout_cells();*/

    /*densities2d.fea_case = &fea_manager.current;
    cout << "\interpolated keep:\n";
    densities2d.visualize_keep_cells();
    cout << "\interpolated cutouts:\n";
    densities2d.visualize_cutout_cells();*/

    run_emma(_evolver, &fea_manager);
}

void Controller::run_emma(Evolver& _evolver, phys::FEACaseManager* fea_manager) {
    // Generic optimizer parameters
    string msh_file = base_folder + "/mesh.msh";
    int max_iterations = 100000;
    bool export_msh = true;
    bool verbose = true;

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
        *fea_manager, mesh, base_folder, pop_size, no_static_iterations_trigger, mutation_rate_level0,
        mutation_rate_level1, densities2d, variation_trigger, max_iterations, max_iterations_without_change,
        export_msh, verbose, initial_perturb_level0, initial_perturb_level1, crossover_method
    );
    _evolver = evolver;
    evolver.evolve();
}

void Controller::run_emma_dynamic() {
    Evolver evolver;
    run_emma_dynamic(evolver);
}

void Controller::run_emma_static() {
    Evolver evolver;
    run_emma_static(evolver);
}

void Controller::run_fess() {
    FESS fess;
    run_fess(fess);
}