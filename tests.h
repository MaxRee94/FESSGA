#pragma once
#include "controller.h"


class Tester {
public:
    Tester() = default;
    Tester(Controller* _ctrl) : ctrl(_ctrl) {};
    void run_tests();
    void init_dummy_fess(FESS& optimizer);
    void init_dummy_evolver(Evolver& evolver);
    bool test_is_single_piece();
    bool do_individual_is_single_piece_test(string type, string path);
    bool do_individual_init_pieces_test(string type, string path, int expected_result, int dim_x = -1, int dim_y = -1, bool verbose = false);
    bool do_individual_remove_smaller_pieces_test(string type, string path, bool expected_result, bool verbose = false);
    bool do_individual_restore_pieces_test(string type, string path, bool verbose = false);
    bool do_individual_remove_isolated_material_test(string type, string path, bool expected_validity, bool verbose = false);
    bool do_individual_fill_voids_test(string type, string path, bool verbose = false);
    bool do_individual_feasibility_filtering_test(string type, string path, bool verbose = false);
    bool do_individual_repair_test(string type, string path, bool verbose = false);
    bool do_individual_init_population_test(string type, string path, bool verbose = false);
    bool do_individual_image_loader_test(string type, string path, bool verbose = true);
    void create_parents(grd::Densities2d parent1, grd::Densities2d parent2);
    bool test_2x_crossover();
    bool test_evolution();
    bool test_init_pieces();
    bool test_remove_smaller_pieces();
    bool test_restore_removed_pieces();
    bool test_remove_isolated_material();
    bool test_fill_voids();
    bool test_feasibility_filtering();
    bool test_repair();
    bool test_init_population();
    bool test_image_loader();
    void do_teardown();
    OptimizerBase do_setup(
        string type, string path, bool verbose = false, int dim_x = -1, int dim_y = -1, string base_folder = ""
    );
    Evolver do_evolver_setup(
        string type, string path, bool verbose = false, int dim_x = -1, int dim_y = -1, string base_folder = ""
    );

    Controller* ctrl = 0;
    Input _input;
    int _dim_x, _dim_y;
    string _base_folder;
};



