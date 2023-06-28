#pragma once
#include "controller.h"


class Tester {
public:
    Tester() = default;
    Tester(Controller* _ctrl) : ctrl(_ctrl) {};
    void run_tests();
    void init_dummy_optimizer(OptimizerBase& optimizer);
    bool test_is_single_piece();
    bool do_individual_is_single_piece_test(string type, string path);
    void create_parents(grd::Densities2d parent1, grd::Densities2d parent2);
    bool test_2d_crossover();
    bool test_evolution();

    Controller* ctrl = 0;
};

void Tester::init_dummy_optimizer(OptimizerBase& optimizer) {
    // Parameters
    double max_stress = 1.5e9;
    double min_stress = 7e3;
    string msh_file = ctrl->output_folder + "/mesh.msh";
    string fea_case = ctrl->output_folder + "/case.sif";
    int max_iterations = 100;
    bool export_msh = true;
    float greediness = 0.1;
    bool verbose = true;
    bool maintain_boundary_cells = true;
    optimizer = OptimizerBase(msh_file, fea_case, ctrl->mesh, ctrl->output_folder, max_stress, ctrl->densities2d, max_iterations, export_msh, verbose);
}

bool Tester::do_individual_is_single_piece_test(string type, string path) {
    // Setup
    Input _input = ctrl->input;
    ctrl->input.path = path;
    ctrl->input.type = type;
    ctrl->init_densities();
    OptimizerBase optimizer;
    init_dummy_optimizer(optimizer);
    optimizer.densities.print();

    // Test
    bool output = optimizer.densities.is_single_piece(ctrl->densities2d.fea_case, &ctrl->densities2d.removed_cells);

    // Teardown
    ctrl->input = _input;

    return output;
}

bool Tester::test_is_single_piece() {
    // Do multi-piece and single-piece tests
    bool success = true;
    success = success && !do_individual_is_single_piece_test("distribution2d", "../data/unit_tests/distribution2d_multi_piece_1.dens");
    success = success && !do_individual_is_single_piece_test("distribution2d", "../data/unit_tests/distribution2d_multi_piece_2.dens");
    success = success && do_individual_is_single_piece_test("distribution2d", "../data/unit_tests/distribution2d_single_piece.dens");

    cout << "TESTING: is_single_piece(). Test " << (success ? "succeeded." : "failed.") << endl;

    return success;
}

void Tester::run_tests() {
    cout << "Running tests...\n";
    int successes = 0;
    int failures = 0;

    bool is_single_piece_success = test_is_single_piece();
    successes += is_single_piece_success;
    failures += !is_single_piece_success;

    cout << "ALL TESTS FINISHED. " << successes << " / " << (failures + successes) << " tests successful.\n";
}


/*
Create 2 parent slices from the 3d binary density distribution for 2d test
*/
void Tester::create_parents(grd::Densities2d parent1, grd::Densities2d parent2) {
    int z1 = ctrl->dim_x / 2;
    int z2 = ctrl->dim_x / 2 + (ctrl->dim_x / 5);
    ctrl->densities3d.create_slice(parent1, 2, z1);
    ctrl->densities3d.create_slice(parent2, 2, z2);
}

/*
Test 2-point crossover of two 2d parent solutions. Print parents and children to console
*/
bool Tester::test_2d_crossover() {
    evo::Individual2d parent1(ctrl->dim_x, ctrl->dim_y, ctrl->mesh.diagonal);
    evo::Individual2d parent2(ctrl->dim_x, ctrl->dim_y, ctrl->mesh.diagonal);
    double max_stress = 1e9; // arbitrary maximum stress
    int max_iterations = 100;
    bool export_msh = false;

    create_parents(parent1, parent2);
    cout << "\nParent 1: \n";
    parent1.print();
    cout << "\nParent 2: \n";
    parent2.print();

    string fea_case = "../data/msh_output/case.sif";
    string msh_file = "../data/msh_output/test.msh";
    string output_folder = "../data/msh_output/FESSGA_test_output";

    // Do crossover
    evo::Individual2d child1(ctrl->dim_x, ctrl->dim_y, ctrl->mesh.diagonal);
    evo::Individual2d child2(ctrl->dim_x, ctrl->dim_y, ctrl->mesh.diagonal);
    Evolver evolver = Evolver(
        msh_file, fea_case, ctrl->mesh, output_folder, 4, (float)0.01, &variation_minimum_passed, 2,
        max_stress, parent1, max_iterations
    );
    evolver.do_2d_crossover(parent1, parent2, child1, child2);

    cout << "\Child 1: \n";
    child1.print();
    cout << "\Child 2: \n";
    child2.print();

    return true;
}

bool Tester::test_evolution() {
    Evolver evolver = Evolver();

    return true;
}

