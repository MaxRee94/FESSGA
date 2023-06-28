#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <map>
#include "physics.h"

using namespace Eigen;
using namespace std;
using namespace fessga;

namespace fessga {
    class grd {
    public:

        struct Piece {
            vector<int> cells;
            bool is_removable = true;
        };

        struct Grid3D {
            Grid3D() = default;
            Grid3D(int dim_x, int dim_y, int dim_z, Vector3d diagonal) {
                x = dim_x;
                y = dim_y;
                z = dim_z;
                cell_size = diagonal.cwiseProduct(
                    Vector3d(1.0 / (double)x, 1.0 / (double)y, 1.0 / (double)z)
                );
                size3d = dim_x * dim_y * dim_z;
                size2d = dim_x * dim_y;
            }
            int x, y, z, size2d, size3d;
            Vector3d cell_size;
        };

        class Densities2d {
        public:
            Densities2d() = default;
            Densities2d(int _dim_x, int _dim_y, Vector2d diagonal) {
                construct_grid(_dim_x, _dim_y);
                cell_size = diagonal.cwiseProduct(
                    Vector2d(1.0 / (double)dim_x, 1.0 / (double)dim_y)
                );
            }
            Densities2d(int _dim_x, int _dim_y, Vector3d diagonal) {
                construct_grid(_dim_x, _dim_y);
                compute_cellsize(diagonal);
            }
            int get_idx(int x, int y) {
                return x * dim_y + y;
            }
            pair<int, int> get_coords(int idx) {
                pair<int, int> coords(idx / dim_y, idx % dim_y);
                return coords;
            }
            uint operator[](int idx) {
                return values[idx];
            }
            uint operator()(int x, int y) {
                return values[get_idx(x, y)];
            }
            void redo_count() {
                _count = 0;
                for (int i = 0; i < size; i++) _count += values[i];
            }
            void update_count() {
                if (_count == -2) redo_count(); // Values have been set using the set() function and can be safely recounted.
                else if (_count == -1) cerr << "Cannot count density values because the array has not been completely filled.\n" << endl;
            }
            void fill(int idx) {
                update_count();
                if (values[idx] != 1) {
                    values[idx] = 1;
                    _count++;
                }
            }
            void fill(vector<int> indices) {
                update_count();
                for (auto& idx : indices) fill(idx);
            }
            // Set a value without updating the _count. Count will be reset to -1.
            void set(int idx, uint value) {
                _count = -2;
                values[idx] = value;
            }
            void replace_values(uint* values_ptr) {
                values = values_ptr;
                redo_count();
            }
            void del(int idx) {
                update_count();
                if (values[idx] != 0) {
                    values[idx] = 0;
                    _count--;
                }
            }
            void del(vector<int> indices) {
                update_count();
                for (auto& idx : indices) del(idx);
            }
            void delete_all() {
                for (int i = 0; i < size; i++) values[i] = 0;
                _count = 0;
            }
            void fill_all() {
                for (int i = 0; i < size; i++) values[i] = 1;
                _count = size;
            }
            void print() {
                for (int y = dim_y - 1; y > -1; y--) {
                    for (int x = 0; x < dim_x; x++) {
                        cout << values[get_idx(x, y)];
                    }
                    cout << endl;
                }
            }
            int count() {
                if (_count == -1) redo_count();
                return _count;
            }
            uint at(int idx) {
                return values[idx];
            }

            bool cell_is_safe_to_delete(int cell_coord, vector<int>* removed_cells, int& no_deleted_neighbors, phys::FEACase* fea_case);
            bool is_single_piece(phys::FEACase* fea_case, vector<int>* removed_cells, int _start_cell = -1, bool verbose = false);
            int get_cell_not_in_vector(vector<int>* cells_vector);
            int get_no_connected_cells(int cell_coord, Piece& piece, phys::FEACase* fea_case = 0, bool verbose = false);
            int get_unvisited_cell(vector<int>* visited_cells);
            int get_unvisited_neighbor_of_removed_cell(vector<int>* visited_cells);
            string do_export(string output_folder);
            vector<int> get_true_neighbors(int idx);
            vector<int> get_true_neighbors(int x, int y);
            vector<int> remove_smaller_pieces(
                phys::FEACase* fea_case, vector<grd::Piece>* pieces, vector<int>* removed_cells, double max_stress_threshold, phys::FEAResults2D* fea_results,
                bool maintain_boundary_cells,bool _remove_largest_piece = true, bool check_if_single_piece = false
            );
            void copy_from(Densities2d* source);
            void do_import(string path, Vector3d diagonal);
            void filter();
            void get_pieces(vector<int>* visited_cells, int cells_left = -1, int _start_cell = -1);
            void load_snapshot();
            bool remove_floating_piece(
                grd::Piece* piece, double max_stress_threshold, phys::FEACase* fea_case, phys::FEAResults2D* fea_results, int& total_no_cells,
                bool maintain_boundary_cells
            );
            void remove_largest_piece(vector<Piece>* pieces, int& max_size);
            int remove_low_stress_cells(int no_cells_to_remove, grd::Piece* smaller_piece = 0);
            void restore_removed_cells(vector<int>* removed_cells);
            void restore_removed_pieces(vector<Piece>* removed_pieces);
            void save_snapshot();

            int dim_x = 0;
            int dim_y = 0;
            int size = 0;
            int no_dimensions = 2;
            Vector2d cell_size;
            phys::FEACase* fea_case = 0;
            phys::FEAResults2D* fea_results = 0;
            vector<int> removed_cells = {};
            vector<Piece>* pieces;

        protected:
            uint* values = 0;
            uint* snapshot = 0;
            int _count = 0;
            int _snapshot_count = 0;

            virtual void construct_grid(int _dim_x, int _dim_y, int _ = 0) {
                dim_x = _dim_x;
                dim_y = _dim_y;
                size = dim_x * dim_y;
                values = new uint[size];
                snapshot = new uint[size];
                delete_all(); // Initialize all values to zero
            }
            virtual void compute_cellsize(Vector3d diagonal) {
                Vector2d _diagonal2d = Vector2d(diagonal(0), diagonal(1));
                cell_size = _diagonal2d.cwiseProduct(
                    Vector2d(1.0 / (double)dim_x, 1.0 / (double)dim_y)
                );
            }
            void _copy(uint* source, uint* target, int source_count, int target_count);
        };

        class Densities3d : public Densities2d {
        public:
            Densities3d(){ no_dimensions = 3; };
            Densities3d(int _dim_x, int _dim_y, int _dim_z, Vector3d diagonal) {
                construct_grid(_dim_x, _dim_y, _dim_z);
                compute_cellsize(diagonal);
            }
            string do_export(string output_folder);
            void generate(Vector3d offset, MatrixXd* V, MatrixXi* F);
            void filter();
            void create_slice(Densities2d& densities2d, int dimension, int offset);

            Vector3d cell_size;
            int dim_z = 0;

        protected:
            void construct_grid(int _dim_x, int _dim_y, int _dim_z) override {
                dim_x = _dim_x;
                dim_y = _dim_y;
                dim_z = _dim_z;
                size = dim_x * dim_y * dim_z;
                values = new uint[size];
                snapshot = new uint[size];
                delete_all(); // Initialize all values to zero
            }
            void compute_cellsize(Vector3d diagonal) override {
                cell_size = diagonal.cwiseProduct(
                    Vector3d(1.0 / (double)dim_x, 1.0 / (double)dim_y, 1.0 / (double)dim_z)
                );
            }
            void fill_cells_inside_mesh(Vector3d offset, MatrixXd* V, MatrixXi* F);
            void create_x_slice(grd::Densities2d& densities2d, int x);
            void create_y_slice(grd::Densities2d& densities2d, int y);
            void create_z_slice(grd::Densities2d& densities2d, int z);
        };
    };
}

