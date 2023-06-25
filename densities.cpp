#include "densities.h"



// Remove floating cells (i.e. cells that have no direct neighbors)
void fessga::grd::Densities2d::filter() {
    cout << "Filtering 2d floating cells..." << endl;
    update_count();
#pragma omp parallel for
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            int filled = values[x * dim_y + y];
            if (!filled) continue;
            int neighbor = 0;
            for (int _x = -1; _x <= 1; _x++) {
                for (int _y = -1; _y <= 1; _y++) {
                    if (_y == 0 && _x == 0) continue;
                    if (x + _x == dim_x || y + _y == dim_y) continue;
                    if (x + _x <= 0 || y + _y <= 0) continue;
                    neighbor = values[(x + _x) * dim_y + (y + _y)];
                    if (neighbor) break;
                }
                if (neighbor) break;
            }
            if (!neighbor) {
                cout << "floating cell detected. Setting to 0" << endl;
                del(x * dim_y + y); // Remove the cell if it has no neighbors
            }
        }
    }
    cout << "Finished filtering floating cells." << endl;
}

/*
* Export 2d density distribution to file
*/
string fessga::grd::Densities2d::do_export(string output_folder) {
    string content = "2\n";
    content += to_string(dim_x) + "\n" + to_string(dim_y) + "\n";
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            content += to_string(values[x * dim_y + y]);
        }
    }
    content += "\n";
    IO::write_text_to_file(content, output_folder + "/distribution.dens");
    return output_folder + "/distribution.dens";
}


/*
* Import 2d density distribution from file
*/
void fessga::grd::Densities2d::do_import(string path, Vector3d diagonal) {
    // Get a vector of strings representing the lines in the file
    vector<string> lines;
    IO::read_file_content(path, lines);

    // Get the number of dimensions of the distribution
    int no_dimensions = stoi(lines[0]);

    // Get the sizes of the grid along each axis
    int dim_x = stoi(lines[1]), dim_y = stoi(lines[2]);
    if (no_dimensions == 3) {
        cerr << "Error: Densities file " << path << " contains an (incompatible) 3d density distribution.\n";
        return;
    }

    // Construct a grid
    int grid_size = dim_x * dim_y;
    construct_grid(dim_x, dim_y);

    // Fill the densities array with the binary values stored in the last line of the file
    string densities_line = lines[no_dimensions + 1];
    for (int i = 0; i < grid_size; i++) {
        set(i, densities_line[i] - '0');
    }
    update_count();

    // Update cell size
    init_cell_size(diagonal);
}

// Return the indices of the 'true neighbors' of the cell at the given coordinates.
// True neighbors are here defined as filled neighbor cells that share a line with the given cell
vector<int> fessga::grd::Densities2d::get_true_neighbors(int x, int y) {
    vector<pair<int, int>> offsets = { pair(0,1), pair(1,0), pair(-1, 0), pair(0, -1) };
    vector<int> true_neighbors;
    for (auto& offset : offsets) {
        int _x = x + offset.first;
        int _y = y + offset.second;
        if (_x == dim_x || _y == dim_y || _x < 0 || _y < 0) continue;
        int neighbor_coord = _x * dim_y + _y;
        if (values[neighbor_coord]) true_neighbors.push_back(neighbor_coord);
    }
    return true_neighbors;
}

// Return the indices of the 'true neighbors' of the cell at the given coordinates.
// True neighbors are here defined as filled neighbor cells that share a line with the given cell
vector<int> fessga::grd::Densities2d::get_true_neighbors(int idx) {
    pair<int, int> coords = get_coords(idx);
    return get_true_neighbors(coords.first, coords.second);
}

// Get number of connected cells of the given cell using a version of floodfill
int fessga::grd::Densities2d::get_no_connected_cells(int cell_coord, Piece& piece, grd::Case* fe_case, bool verbose) {
    piece.cells = { cell_coord };
    int i = 0;
    while (i < piece.cells.size()) {
        vector<int> neighbors = get_true_neighbors(piece.cells[i]);
        for (int j = 0; j < neighbors.size(); j++) {
            if (!help::is_in(&piece.cells, neighbors[j])) {
                piece.cells.push_back(neighbors[j]);
                if (fe_case != 0 && piece.is_removable) {
                    // If the piece was flagged as removable but one of its cells is a boundary cell, flag it as non-removable.
                    if (help::is_in(&fe_case->boundary_cells, neighbors[j])) {
                        piece.is_removable = false;
                    }
                }
            }
        }
        i++;
    }
    return piece.cells.size();
}

// Return whether the cell at the given coordinates is safe to remove. Also remove neighbor cells that become invalid as a result of
// deleting the given cell.
bool fessga::grd::Densities2d::cell_is_safe_to_delete(
    int cell_coord, vector<int>* removed_cells, int& no_deleted_neighbors, grd::Case* fe_case
) {
    vector<int> neighbors = get_true_neighbors(cell_coord);
    for (auto& neighbor : neighbors) {
        vector<int> sub_neighbors = get_true_neighbors(neighbor);

        // If the neighboring cell has only one true neighbor itself, deleting the current cell would make it float in mid-air and thus invalid.
        // Therefore we either delete the neighboring cell too, or - in case the neighbor is a bound condition cell - skip deletion alltogether.
        if (sub_neighbors.size() <= 1) {
            // If the cell has a line on which a boundary condition was applied, skip deletion
            if (help::is_in(&fe_case->boundary_cells, neighbor) || help::is_in(&fe_case->whitelisted_cells, neighbor)) {
                return false;
            }
            no_deleted_neighbors++;
            del(neighbor); // Delete the neighboring cell, since deleting the cell at <cell_coord> would make it invalid.
            removed_cells->push_back(neighbor);
        }
    }
    return true;
}

// Return cell that was not yet visited
int fessga::grd::Densities2d::get_unvisited_cell(
    vector<int>* visited_cells, vector<int>* removed_cells
) {
    int unvisited_cell = -1;
    for (auto& removed_cell : (*removed_cells)) {
        // At least one of the neighbors of the last-removed cells must belong to the other piece
        vector<int> neighbors = get_true_neighbors(removed_cell);
        for (auto& neighbor : neighbors) {
            if (!help::is_in(visited_cells, neighbor)) {
                unvisited_cell = neighbor;
                return unvisited_cell;
            }
        }
    }
    return unvisited_cell;
}

// Get the pieces inside the density distribution
void fessga::grd::Densities2d::get_pieces(
    grd::Case* fe_case, vector<Piece>* pieces, vector<int>* visited_cells, int& cells_left,
    vector<int>* removed_cells, int& no_pieces, int _start_cell
) {
    Piece piece;
    int start_cell = -1;
    if (_start_cell > -1) start_cell = _start_cell;
    else {
        // If no start cell was provided as an argument, pick one of the neighbors of one of the removed cells
        for (int i = 0; i < removed_cells->size(); i++) {
            vector<int> neighbors = get_true_neighbors(removed_cells->at(i));
            if (neighbors.size() > 0) { start_cell = neighbors[0]; break; }
        }
    }
    if (start_cell == -1) start_cell = fe_case->boundary_cells[0];
    int piece_size = get_no_connected_cells(start_cell, piece, fe_case);
    pieces->push_back(piece);

    // Check if the shape consists of one piece or several
    if (piece_size < cells_left) {
        no_pieces++;
        cells_left -= piece_size;
        for (int i = 0; i < piece_size; i++) visited_cells->push_back(piece.cells[i]);

        // Recurse if there are still unvisited cells left
        if (cells_left > 0) {
            int unvisited_cell = get_unvisited_cell(visited_cells, removed_cells);
            if (unvisited_cell == -1) return;
            get_pieces(fe_case, pieces, visited_cells, cells_left, removed_cells, no_pieces, unvisited_cell);
        }
    }
}

// Restore all cells that were removed
void fessga::grd::Densities2d::restore_removed_cells(vector<int>* removed_cells) {
    for (auto& cell : (*removed_cells)) {
        fill(cell);
    }
}

// Restore all cells in the provided pieces
void fessga::grd::Densities2d::restore_removed_pieces(vector<grd::Piece>* removed_pieces) {
    for (auto& piece : (*removed_pieces)) {
        for (auto& cell : piece.cells) {
            fill(cell);
        }
    }
}

// Return whether or not the shape consists of a single piece
bool fessga::grd::Densities2d::is_single_piece(
    grd::Case* fe_case, int& total_no_cells, vector<int>* removed_cells, int _start_cell, bool verbose
) {
    grd::Piece piece;
    int start_cell;
    if (_start_cell > -1) start_cell = _start_cell;
    else start_cell = fe_case->boundary_cells[0];
    int piece_size = get_no_connected_cells(start_cell, piece);
    if (verbose) {
        cout << "\npiece size: " << piece_size << endl;
        cout << "total no cells: " << total_no_cells << endl;
    }

    // Check if the shape consists of one piece or multiple
    if (piece_size < total_no_cells) {
        int unvisited_cell = get_unvisited_cell(&piece.cells, removed_cells);
        if (unvisited_cell < 0 && (total_no_cells - piece_size == 1 || piece_size == 1)) return true;

        return false;
    }
    else return true;
}

// Remove the largest piece from the given pieces vector
void fessga::grd::Densities2d::remove_largest_piece(vector<grd::Piece>* pieces, int& max_size) {
    max_size = 0;
    int largest_piece_idx = -1;
    for (int i = 0; i < pieces->size(); i++) {
        if (pieces->at(i).cells.size() > max_size) {
            max_size = pieces->at(i).cells.size();
            largest_piece_idx = i;
        }
    }
    pieces->erase(pieces->begin() + largest_piece_idx);
}

int fessga::grd::Densities2d::remove_low_stress_cells(
    PairSet* fe_results, grd::Case* fe_case, int no_cells_to_remove,
    vector<int>& removed_cells, double max_stress_threshold, grd::Piece* smaller_piece, int total_no_cells
) {
    int count = 0;
    int cell_from_smaller_piece;
    int no_iterations_without_removal = -1;
    for (auto& item : (*fe_results)) {
        no_iterations_without_removal++;
        if (no_iterations_without_removal > 100) {
            cout << "WARNING: 100 cells in a row could not be removed. There may be no more cells to remove.\n";
        }
        int cell_coord = item.first;
        double cell_stress = item.second;

        // If the cell's stress exceeds the maximum, break the loop (all subsequent cells will also exceed the maximum since the list is ordered).
        if (cell_stress > max_stress_threshold) continue;

        // If the cell has a line on which a boundary condition was applied, skip deletion
        if (help::is_in(&fe_case->boundary_cells, cell_coord)) {
            continue;
        }

        // If the cell was previously whitelisted, skip deletion
        if (help::is_in(&fe_case->whitelisted_cells, cell_coord)) {
            continue;
        }

        // If the cell was already empty, skip deletion (more importantly: don't count this as a deletion)
        if (!values[cell_coord]) continue;

        // If a 'smaller piece' vector was provided, perform cell removal in 'careful mode'. This means: check whether cell deletion results in 
        // multiple pieces. If so, undo the deletion and whitelist the cell.
        if (smaller_piece != 0) {
            vector<int> _removed_cell = { cell_coord };
            //cout << "\ntotal no cells: " << total_no_cells - count << endl;
            //cout << "is single piece before element removal: " << mesher::is_single_piece(
            //    densities, grid, fe_case, total_no_cells - count, &_removed_cell, smaller_piece->cells[0]) << endl;
            del(cell_coord);
            cell_from_smaller_piece = smaller_piece->cells[0];
            bool verbose = no_iterations_without_removal > 100;
            bool _is_single_piece = is_single_piece(
                fe_case, total_no_cells, &_removed_cell, cell_from_smaller_piece, verbose
            );
            //if (verbose) cout << "is single piece after deleting cell " << cell_coord / grid.y << ", " << cell_coord % grid.y << " ? " << _is_single_piece << endl;
            if (!_is_single_piece) {
                if (no_iterations_without_removal > 100) {
                    set(cell_coord, 5); // TEMP
                    print(); // TEMP
                }
                // Cell removal resulted in multiple pieces. Restore cell and whitelist it.
                fill(cell_coord);
                //cout << "restored cell." << endl;
                fe_case->whitelisted_cells.push_back(cell_coord);
                continue;
            }
            if (count % (no_cells_to_remove / 10) == 0) cout << "no removed cells: " << count << " / " << no_cells_to_remove << endl;
            count++;
            total_no_cells--;
        }

        // If cell deletion leads to infeasibility, skip deletion
        int no_deleted_neighbors = 0;
        if (!cell_is_safe_to_delete(cell_coord, &removed_cells, no_deleted_neighbors, fe_case)) {
            fill(cell_coord);
            //cout << "cell is not safe to delete" << endl;
            continue;
        }
        count += no_deleted_neighbors;
        total_no_cells -= no_deleted_neighbors;

        if (smaller_piece != 0) {
            vector<int> _removed_cell = { cell_coord };
            //cout << "total no cells: " << total_no_cells - count << endl;
            //cout << "is single piece after neighbor checking / deletion: " << mesher::is_single_piece(
            //    densities, grid, fe_case, total_no_cells - count, &_removed_cell, smaller_piece->cells[0]) << endl;
        }

        // Set cell to zero, making it empty
        del(cell_coord);
        if (!smaller_piece) count++;
        removed_cells.push_back(cell_coord);
        no_iterations_without_removal = 0;
        if (count >= no_cells_to_remove) break;
    }
    fe_case->whitelisted_cells.clear();
    return count;
}

// Remove a floating piece (if possible)
bool fessga::grd::Densities2d::remove_floating_piece(
    grd::Piece* piece, double max_stress_threshold, grd::Case* fe_case, FEResults2D* fe_results,
    int& total_no_cells, bool maintain_boundary_cells
) {
    vector<int> removed_cells;
    for (auto& cell : piece->cells) {
        // TODO: The last two conditions of the following if-statement should not be necessary, but they are. Figure out why.
        bool cell_cannot_be_removed = fe_results->data_map[cell] > max_stress_threshold ||
            (help::is_in(&fe_case->boundary_cells, cell) || help::is_in(&fe_case->whitelisted_cells, cell)
                );
        if (cell_cannot_be_removed) {
            if (maintain_boundary_cells) {
                restore_removed_cells(&removed_cells);
                total_no_cells += removed_cells.size();
                piece->is_removable = false;
                return false;
            }
        }
        else {
            removed_cells.push_back(cell);
            total_no_cells -= 1;
            del(cell);
        }
    }
    return true;
}

// Remove all pieces smaller than the largest piece from the mesh, if possible
vector<int> fessga::grd::Densities2d::remove_smaller_pieces(
    grd::Case* fe_case, int total_no_cells, vector<grd::Piece>* pieces, vector<int>* removed_cells,
    double max_stress_threshold, FEResults2D* fe_results, bool maintain_boundary_cells,
    bool _remove_largest_piece, bool check_if_single_piece
) {
    bool verbose = false;
    vector<int> removed_piece_indices;
    int size_largest_piece;
    if (_remove_largest_piece) {
        cout << "DENSITIES BEFORE ANY PIECES ARE REMOVED " << endl;
        remove_largest_piece(pieces, size_largest_piece);
    }
    if (verbose) print();

    for (int i = 0; i < pieces->size(); i++) {
        grd::Piece piece = pieces->at(i);
        bool piece_removed = false;
        if (piece.is_removable || !maintain_boundary_cells) {
            piece_removed = remove_floating_piece(
                &piece, max_stress_threshold, fe_case, fe_results, total_no_cells, maintain_boundary_cells
            );
        }

        if (piece_removed) {
            cout << "Removed piece " << i + 1 << " / " << pieces->size() << endl;
            if (verbose) {
                print();
            }
            removed_piece_indices.push_back(i);
        }
        else {
            cout << "Piece " << i + 1 << " was not removed." << endl;
        }

        // If this flag is set, we want to prevent that the removal of a piece results in the splitting of the shape into multiple pieces.
        // This can be the case if we've just restored the naively removed cells, and are trying to re-remove the pieces that could be successfully
        // removed.
        if (check_if_single_piece) {
            bool _is_single_piece = is_single_piece(fe_case, total_no_cells, removed_cells, -1, true);
            vector<grd::Piece> removed_piece = { piece };
            cout << "CHECKING WHETHER PIECE REMOVAL RESULTED IN MULTIPLE PIECES....\n";

            // If the shape is broken into multiple pieces, undo the removal of the current piece
            if (!_is_single_piece) {
                cout << "UNDOING THE REMOVAL OF PIECE " << i + 1 << " BECAUSE ITS REMOVAL SPLIT THE SHAPE INTO MULTIPLE PIECES" << endl;
                restore_removed_pieces(&removed_piece);
                total_no_cells += piece.cells.size();
                removed_piece_indices.erase(removed_piece_indices.begin() + removed_piece_indices.size() - 1);
            }
        }
    }

    return removed_piece_indices;
}

// Save a copy of the current state of the density distribution
void fessga::grd::Densities2d::save_snapshot() {
    for (int i = 0; i < size; i++) snapshot[i] = values[i];
    _snapshot_count = _count;
}

// Load the previously saved state (i.e. the snapshot)
void fessga::grd::Densities2d::load_snapshot() {
    for (int i = 0; i < size; i++) values[i] = snapshot[i];
    _count = _snapshot_count;
}

// Copy the density values from the given Densities2d-object to the current object
void fessga::grd::Densities2d::copy_from(Densities2d* source) {
    for (int i = 0; i < size; i++) values[i] = source->at(i);
    _count = source->count();
}
