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

        class Piece {
        public:
            Piece() { id = help::get_rand_uint(0, 1e9); };
            vector<int> cells;
            bool is_removable = true;
            int id;
            bool is_main_piece = false;
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
            Densities2d(int _dim_x, int _dim_y, Vector3d diagonal, phys::FEAResults2D* _fea_results, phys::FEACase* _fea_case) {
                construct_grid(_dim_x, _dim_y);
                compute_cellsize(diagonal);
                fea_results = _fea_results;
                fea_case = _fea_case;
            }
            int get_idx(int x, int y) {
                return x * dim_y + y;
            }
            pair<int, int> get_coords(int cell) {
                pair<int, int> coords(cell / dim_y, cell % dim_y);
                return coords;
            }
            uint operator[](int cell) {
                return values[cell];
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
            void fill(int cell) {
                update_count();
                if (values[cell] != 1) {
                    values[cell] = 1;
                    _count++;
                }
            }
            void fill(vector<int> cells) {
                update_count();
                for (auto& cell : cells) fill(cell);
            }
            // Set a value without updating the _count. Count will be reset to -1.
            void set(int cell, uint value) {
                _count = -2;
                values[cell] = value;
            }
            void replace_values(uint* values_ptr) {
                values = values_ptr;
                redo_count();
            }
            // Delete cell, update count, but do not record the removal
            void del(int cell) {
                update_count();
                if (values[cell] != 0) {
                    values[cell] = 0;
                    _count--;
                }
            }
            // Delete cells, update count, but do not record the removal
            void del(vector<int> cells) {
                update_count();
                for (auto& cell : cells) del(cell);
            }
            // Delete cell, update count, AND record the removal
            void remove_and_remember(int cell) {
                update_count();
                if (values[cell] != 0) {
                    values[cell] = 0;
                    _count--;
                }
                removed_cells.push_back(cell);
            }
            // Delete cells, update count, AND record the removal
            void remove_and_remember(vector<int> cells) {
                update_count();
                for (auto& cell : cells) {
                    del(cell);
                    removed_cells.push_back(cell);
                }
            }
            void set_main_piece(Piece* piece) {
                main_piece = *piece;
                for (auto& _piece : pieces) _piece.is_main_piece = false;
                piece->is_main_piece = true;
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
            uint at(int cell) {
                return values[cell];
            }
            void restore(int cell) {
                fill(cell);
                help::remove(&removed_cells, cell);
            }
            int get_piece_index(vector<Piece>* _pieces, Piece* piece) {
                for (int i = 0; i < _pieces->size(); i++) {
                    if (piece->id == _pieces->at(i).id) return i;
                }
                throw("ERROR: piece not found\n");
            }

            void move_piece_to_trash(Piece* piece);
            void move_piece_from_trash(Piece* piece);
            bool cell_is_safe_to_delete(int cell_coord, vector<int>* removed_cells, int& no_deleted_neighbors, phys::FEACase* fea_case);
            bool is_single_piece(int _start_cell = -1, bool verbose = false);
            int get_cell_not_in_vector(vector<int>* cells_vector);
            int get_no_connected_cells(int cell_coord, Piece& piece, phys::FEACase* fea_case = 0, bool verbose = false);
            int get_unvisited_cell(vector<int>* visited_cells);
            bool is_in(vector<Piece>* pieces, Piece* piece);
            int get_unvisited_neighbor_of_removed_cell(vector<int>* visited_cells);
            string do_export(string output_path);
            vector<int> get_true_neighbors(int idx);
            vector<int> get_true_neighbors(int x, int y);
            void remove_smaller_pieces(
                vector<grd::Piece> pieces, vector<int>* removed_cells, bool _remove_largest_piece_from_vector = true, bool check_if_single_piece = false
            );
            void remove_smaller_pieces();
            void copy_from(Densities2d* source);
            void do_import(string path, Vector3d diagonal);
            void filter();
            void init_pieces(int _start_cell = -1);
            void load_snapshot();
            bool remove_floating_piece(grd::Piece* piece);
            void remove_largest_piece_from_vector(vector<Piece>* pieces, int& max_size);
            int remove_low_stress_cells(int no_cells_to_remove, int no_cells_removed, grd::Piece* smaller_piece = 0);
            void restore_removed_cells(vector<int> removed_cells);
            void restore_removed_pieces(vector<Piece> pieces_to_restore);
            void save_snapshot();
            void flush_edit_memory();
            int get_no_cells_in_removed_pieces();
            void visualize_bound_cells();

            int dim_x = 0;
            int dim_y = 0;
            int size = 0;
            int no_dimensions = 2;
            Vector2d cell_size;
            phys::FEACase* fea_case = 0;
            phys::FEAResults2D* fea_results = 0;
            vector<int> removed_cells = {};
            vector<Piece> pieces;
            vector<Piece> removed_pieces;
            Piece main_piece = Piece();

        protected:
            uint* values = 0;
            uint* snapshot = 0;
            int _count = 0;
            int _snapshot_count = 0;

            void init_pieces(vector<int>* visited_cells, int cells_left, int _start_cell);
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
            Densities3d(int _dim_x, int _dim_y, int _dim_z, Vector3d diagonal, phys::FEAResults2D* _fea_results, phys::FEACase* _fea_case) {
                construct_grid(_dim_x, _dim_y, _dim_z);
                compute_cellsize(diagonal);
                fea_results = _fea_results;
                fea_case = _fea_case;
            }
            string do_export(string output_path);
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

