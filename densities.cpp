#include "densities.h"



// Remove floating cells (i.e. cells that have no direct neighbors)
void fessga::grd::Densities2d::filter(int no_neighbors, bool restore_bound_cells) {
    save_internal_snapshot();
    vector<int> x_bounds = { 0, dim_x - 1 };
    vector<int> y_bounds = { 0, dim_y - 1 };
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            int cell = x * dim_y + y;

            // If the cell is already empty, continue to the next cell
            if (!values[cell]) continue;

            // Delete the cell if the number of neighbors is equal to the given number
            vector<int> neighbors = get_neighbors(x, y, snapshot_internal);
            if (neighbors.size() == no_neighbors) {
                del(cell);
            }
        }
    }
    if (restore_bound_cells) {
        // Restore cells marked as 'to keep'
        for (auto& cell : fea_casemanager->cells_to_keep) {
            fill(cell);
        }
    }
}

/*
* Export 2d density distribution to given path
*/
string fessga::grd::Densities2d::do_export(string output_path) {
    string content = "2\n";
    content += to_string(dim_x) + "\n" + to_string(dim_y) + "\n";
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            content += to_string(values[x * dim_y + y]);
        }
    }
    content += "\n";
    IO::write_text_to_file(content, output_path);
    return output_path;
}

void fessga::grd::Densities2d::init_vtk_paths() {
    vtk_paths.clear();
    for (auto& fea_case : fea_casemanager->active_cases) {
        vtk_paths.push_back(output_folder + "/" + fea_case.name + "_0001.vtk");
    }
}


/*
* Export 3d density distribution to given path
*/
string fessga::grd::Densities3d::do_export(string output_path) {
    string content = "3\n";
    content += to_string(dim_x) + "\n" + to_string(dim_y) + "\n" + to_string(dim_z) + "\n";
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            for (int z = 0; z < dim_z; z++) {
                content += to_string(values[x * dim_y * dim_z + y * dim_z + z]);
            }
        }
    }
    content += "\n";
    IO::write_text_to_file(content, output_path);
    return output_path;
}

/*
* Import density distribution from file
*/
void fessga::grd::Densities2d::do_import(string path, float width) {
    // Get a vector of strings representing the lines in the file
    vector<string> lines;
    IO::read_file_content(path, lines);

    // Get the number of dimensions of the distribution
    int _no_dimensions = stoi(lines[0]);
    int _dim_x = stoi(lines[1]); int _dim_y = stoi(lines[2]);
    if (_no_dimensions == 3) {
        if (no_dimensions == 2) throw std::runtime_error("Error: Cannot load a 3d density distribution using a Densities2d object. Use Densities3d instead.\n");
        else {
            int _dim_z = stoi(lines[3]);
            compute_diagonal_and_cellsize(width, _dim_x, _dim_y, _dim_z);
        }
    }
    else compute_diagonal_and_cellsize(width, _dim_x, _dim_y);

    // Re-initialize grid
    construct_grid();

    // Fill the densities array with the binary values stored in the last line of the file
    string densities_line = lines[no_dimensions + 1];
    for (int i = 0; i < size; i++) {
        set(i, densities_line[i] - '0');
    }
    update_count();
}


// Return the indices of the 'true neighbors' of the cell at the given coordinates.
// True neighbors are here defined as filled neighbor cells that share a line with the given cell
vector<int> fessga::grd::Densities2d::get_neighbors(int x, int y, uint* _values) {
    vector<pair<int, int>> offsets = { pair(0,1), pair(1,0), pair(-1, 0), pair(0, -1) };
    vector<int> true_neighbors;
    if (_values == 0) _values = values;
    for (auto& offset : offsets) {
        int _x = x + offset.first;
        int _y = y + offset.second;
        if (_x == dim_x || _y == dim_y || _x < 0 || _y < 0) continue;
        int neighbor_coord = _x * dim_y + _y;
        if (_values[neighbor_coord] == 1) true_neighbors.push_back(neighbor_coord);
    }
    return true_neighbors;
}

// Return the indices of the 'true neighbors' of the cell at the given coordinates.
// True neighbors are here defined as filled neighbor cells that share a line with the given cell
vector<int> fessga::grd::Densities2d::get_neighbors(int idx) {
    pair<int, int> coords = get_coords(idx);
    return get_neighbors(coords.first, coords.second);
}

// Return the indices of the cells neighboring the given cell that are currently void
vector<int> fessga::grd::Densities2d::get_empty_neighbors(int x, int y, bool get_diagonal_neighbors) {
    vector<pair<int, int>> offsets = { pair(0,1), pair(1,0), pair(-1, 0), pair(0, -1) };
    if (get_diagonal_neighbors) {
        help::append_vector(offsets, { pair(-1,-1), pair(-1, 1), pair(1,1), pair(1,-1) });
    }
    vector<int> void_neighbors;
    for (auto& offset : offsets) {
        int _x = x + offset.first;
        int _y = y + offset.second;
        if (_x == dim_x || _y == dim_y || _x < 0 || _y < 0) continue;
        int neighbor_coord = _x * dim_y + _y;
        if (values[neighbor_coord] == 0) void_neighbors.push_back(neighbor_coord);
    }
    return void_neighbors;
}

// Return the indices of the cells neighboring the given cell that are currently void
vector<int> fessga::grd::Densities2d::get_empty_neighbors(int idx, bool get_diagonal_neighbors) {
    pair<int, int> coords = get_coords(idx);
    return get_empty_neighbors(coords.first, coords.second, get_diagonal_neighbors);
}

// Get number of connected cells of the given cell using a version of floodfill
int fessga::grd::Densities2d::get_no_connected_cells(int cell_coord, Piece& piece, bool verbose) {
    piece.cells = { cell_coord };
    int i = 0;
    while (i < piece.cells.size()) {
        vector<int> neighbors = get_neighbors(piece.cells[i]);
        for (int j = 0; j < neighbors.size(); j++) {
            if (!help::is_in(&piece.cells, neighbors[j])) {
                piece.cells.push_back(neighbors[j]);
                if (piece.is_removable) {
                    // If the piece was flagged as removable but one of its cells is a boundary cell, flag it as non-removable.
                    if (fea_casemanager && help::is_in(&fea_casemanager->cells_to_keep, neighbors[j])) {
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
bool fessga::grd::Densities2d::cell_is_safe_to_delete(int cell_coord, int& no_deleted_neighbors) {
    vector<int> neighbors = get_neighbors(cell_coord);
    for (auto& neighbor : neighbors) {
        vector<int> sub_neighbors = get_neighbors(neighbor);

        // If the neighboring cell has only one true neighbor itself, deleting the current cell would make it float in mid-air and thus invalid.
        // Therefore we either delete the neighboring cell too, or - in case the neighbor is a bound condition cell - skip deletion alltogether.
        if (sub_neighbors.size() <= 1) {
            // If the cell has a line on which a boundary condition was applied, skip deletion
            if (help::is_in(&fea_casemanager->cells_to_keep, neighbor)) {
                return false;
            }
            no_deleted_neighbors++;
            remove_and_remember(neighbor); // Delete the neighboring cell, since deleting the cell at <cell_coord> would make it invalid.
        }
    }
    return true;
}

int fessga::grd::Densities2d::get_unvisited_neighbor_of_removed_cell(vector<int>* visited_cells) {
    int unvisited_cell = -1;
    for (auto& removed_cell : removed_cells) {
        // At least one of the neighbors of the last-removed cells must belong to the other piece
        vector<int> neighbors = get_neighbors(removed_cell);
        for (auto& neighbor : neighbors) {
            if (!help::is_in(visited_cells, neighbor)) {
                if (VERBOSE) cout << "neighbor " << neighbor << " of removed cell " << removed_cell << " was not visited yet.\n";
                unvisited_cell = neighbor;
                return unvisited_cell;
            }
        }
    }
    return unvisited_cell;
}

// Return a filled cell that is not in the given vector of cells
int fessga::grd::Densities2d::get_cell_not_in_vector(vector<int>* cells_vector) {
    for (int cell = 0; cell < size; cell++) {
        if (values[cell] != 1) continue;
        bool cell_is_member = help::is_in(cells_vector, cell);
        if (!cell_is_member) return cell;
    }
    return -1;
}


// Return a filled cell that was not yet visited
int fessga::grd::Densities2d::get_unvisited_cell(vector<int>* visited_cells) {
    int unvisited_cell = -1;
    //if (removed_cells.size() > 0) unvisited_cell = get_unvisited_neighbor_of_removed_cell(visited_cells);
    if (unvisited_cell == -1) return get_cell_not_in_vector(visited_cells);
}


// Get the pieces inside the density distribution (public function)
void fessga::grd::Densities2d::init_pieces(int _start_cell) {
    vector<int> visited_cells;
    int cells_left = count();
    pieces.clear();
    int start_cell = -1;
    if (_start_cell != -1) start_cell = _start_cell;
    else if (removed_cells.size() > 0) {
        // If no start cell was provided as an argument, pick one of the neighbors of a removed cell (if a removed cell is available)
        for (int i = 0; i < removed_cells.size(); i++) {
            vector<int> neighbors = get_neighbors(removed_cells.at(i));
            if (neighbors.size() > 0) { start_cell = neighbors[0]; break; }
        }
    }
    else start_cell = fea_casemanager->cells_to_keep[0]; // If no removed cells were stored, pick a cell on which a boundary condition was applied
    init_pieces(&visited_cells, cells_left, start_cell);
}

// Get the pieces inside the density distribution (protected internal function)
void fessga::grd::Densities2d::init_pieces(vector<int>* visited_cells, int cells_left, int start_cell) {
    Piece piece = Piece();
    int piece_size = get_no_connected_cells(start_cell, piece);
    pieces.push_back(piece);

    // Check if the shape consists of one piece or several
    if (piece_size < cells_left) {
        // If the piece has boundary cells and it is the largest piece, it is considered the main piece
        if (
            fea_casemanager && piece.cells.size() > main_piece.cells.size() &&
            help::have_overlap(&piece.cells, &fea_casemanager->cells_to_keep)
            ) {
            set_main_piece(&piece);
        }
        cells_left -= piece_size;
        for (int i = 0; i < piece_size; i++) visited_cells->push_back(piece.cells[i]);

        // Recurse if there are still unvisited cells left
        if (cells_left > 0 && visited_cells->size() < count()) {
            int unvisited_cell = get_unvisited_cell(visited_cells);
            //cout << "unvisited cell (" << unvisited_cell << ")\n";
            /*set(unvisited_cell, 5);
            print();*/
            if (unvisited_cell == -1) {
                redo_count();
                if (visited_cells->size() >= count()) return;
                for (auto& cell : *visited_cells) set(cell, 5);
                cout << "Error: No unvisited cell found, but there are still cells left that haven't been visited (" + to_string(cells_left) + ")\n";
                print();
                for (auto& cell : *visited_cells) set(cell, 1);
                //throw("Error: No unvisited cell found, but there are still cells left that haven't been visited (" + to_string(cells_left) + ")\n");
                return;
            }
            init_pieces(visited_cells, cells_left, unvisited_cell);
        }
    }
}

// Visualize distribution and highlight keep cells
void fessga::grd::Densities2d::visualize_keep_cells() {
    //for (auto& cell : fea_case->cells_to_keep) values[cell] ? set(cell, 5) : throw("Error: Keep cell not filled.\n");
    for (auto& cell : fea_casemanager->cells_to_keep) set(cell, 5);
    print();
    for (auto& cell : fea_casemanager->cells_to_keep) set(cell, 1);
}

// Visualize distribution and highlight cutout cells
void fessga::grd::Densities2d::visualize_cutout_cells() {
    //for (auto& cell : fea_casemanager->cutout_cells) !values[cell] ? set(cell, 5) : throw("Error: Cutout cell is filled.\n");
    for (auto& cell : fea_casemanager->cutout_cells) {
        if (!values[cell]) set(cell, 5);
        else set(cell, 8);
    }
    print();
    for (auto& cell : fea_casemanager->cutout_cells) set(cell, 0);
}

// Visualize distribution and highlight removed cells
void fessga::grd::Densities2d::visualize_removed_cells() {
    //for (auto& cell : fea_casemanager->removed_cells) !values[cell] ? set(cell, 5) : throw("Error: 'Removed' cell is filled.\n");
    for (auto& cell : removed_cells) {
        if (!values[cell]) set(cell, 5);
        else set(cell, 8);
    }
    print();
    for (auto& cell : removed_cells) set(cell, 0);
}

bool fessga::grd::Densities2d::is_in(vector<grd::Piece>* _pieces, grd::Piece* piece) {
    for (int i = 0; i < _pieces->size(); i++) {
        if (_pieces->at(i).id == piece->id) return true;
    }
    return false;
}

// Get the fraction (#no_filled_cells / #total_no_cells)
double fessga::grd::Densities2d::get_relative_area() {
    return (double)count() / (double)size;
}

// Get the fraction (#total_no_cells / #no_filled_cells)
double fessga::grd::Densities2d::get_inverse_relative_volume() {
    return (double)size / (double)count();
}

// Restore all cells that were removed
void fessga::grd::Densities2d::restore_removed_cells(vector<int> _removed_cells) {
    for (auto& cell : _removed_cells) {
        restore(cell);
    }
}

// Restore all cells in the provided pieces
void fessga::grd::Densities2d::restore_removed_pieces(vector<grd::Piece> pieces_to_restore) {
    for (int i = 0; i < pieces_to_restore.size(); i++) {
        Piece piece = pieces_to_restore.at(i);
        fill(piece.cells);
        move_piece_from_trash(&piece);
    }
}

void fessga::grd::Densities2d::move_piece_to_trash(Piece* piece) {
    removed_pieces.push_back(*piece);
    /*cout << "Moving piece with id " << piece->id << " to trash, existing piece ids: ";
    for (auto& piece : pieces) cout << piece.id << ", ";
    cout << endl;*/
    int idx = get_piece_index(&pieces, piece);
    pieces.erase(pieces.begin() + idx);
}

void fessga::grd::Densities2d::move_piece_from_trash(Piece* piece) {
    pieces.push_back(*piece);
    int idx = get_piece_index(&removed_pieces, piece);
    removed_pieces.erase(removed_pieces.begin() + idx);
}

// Return whether or not the shape consists of a single piece
bool fessga::grd::Densities2d::is_single_piece(int _start_cell, bool verbose) {
    Piece piece = Piece();
    int start_cell;
    if (_start_cell > -1) start_cell = _start_cell;
    else start_cell = fea_casemanager->cells_to_keep[0];
    int piece_size = get_no_connected_cells(start_cell, piece);
    if (verbose) {
        cout << "\npiece size: " << piece_size << endl;
        cout << "total no cells: " << count() << endl;
    }

    // Check if the shape consists of one piece or multiple
    if (piece_size < count()) {
        int unvisited_cell = get_unvisited_cell(&piece.cells);
        if (unvisited_cell < 0 && (count() - piece_size == 1 || piece_size == 1)) return true;

        return false;
    }
    else return true;
}

// Remove the largest piece from the given pieces vector
void fessga::grd::Densities2d::remove_largest_piece_from_vector(vector<grd::Piece>* _pieces, int& max_size) {
    max_size = 0;
    int main_piece_idx = -1;
    for (int i = 0; i < _pieces->size(); i++) {
        if (_pieces->at(i).cells.size() > max_size) {
            max_size = _pieces->at(i).cells.size();
            main_piece_idx = i;
        }
    }
    _pieces->erase(_pieces->begin() + main_piece_idx);
}

void fessga::grd::Densities2d::flush_edit_memory() {
    removed_cells.clear();
    removed_pieces.clear();
    main_piece = grd::Piece();
}

int fessga::grd::Densities2d::remove_low_stress_cells(
    int no_cells_to_remove, int no_cells_removed, grd::Piece* smaller_piece
) {
    int cell_from_smaller_piece;
    int no_iterations_without_removal = -1;
    int initial_count = count();
    for (auto [cell, cell_stress] : fea_results.data) {
        no_iterations_without_removal++;
        if (no_iterations_without_removal > 200) {
            if (true) cout << "WARNING: 200 cells in a row could not be removed. "
                << "There may be no more cells to remove.\n";
            return removed_cells.size();
        }

        // If the cell's stress exceeds the maximum, break the loop (all subsequent cells will also exceed the
        // maximum since the list is ordered).
        if (cell_stress > fea_casemanager->max_stress_threshold) break;

        //if (no_iterations_without_removal > 100) cout << "before boundcells\n";
        
        // If the cell has a line on which a boundary condition was applied, skip deletion
        if (help::is_in(&fea_casemanager->cells_to_keep, cell)) {
            continue;
        }

        //if (no_iterations_without_removal > 100) cout << "before value check\n";

        // If the cell was already empty, skip deletion (more importantly: don't count this as a deletion)
        if (!values[cell]) continue;

        //if (no_iterations_without_removal > 100) cout << "before removal\n";

        // If a 'smaller piece' vector was provided, perform cell removal in 'careful mode'. This means: check
        // whether cell deletion results in multiple pieces. If so, undo the deletion and whitelist the cell.
        if (smaller_piece != 0) {
            vector<int> _removed_cell = { cell };
            remove_and_remember(cell);
            cell_from_smaller_piece = smaller_piece->cells[0];
            int total_no_cells = count();
            bool _is_single_piece = is_single_piece(cell_from_smaller_piece, VERBOSE);
            if (!_is_single_piece) {

                // Try removing the smaller pieces
                // Flush edit memory except for removed_cells
                pieces.clear(); removed_pieces.clear(); main_piece = grd::Piece(); 
                init_pieces();
                vector<int> removed_cell = { cell };
                remove_smaller_pieces(pieces, &removed_cell, true);
                if (pieces.size() == 1) {
                    //cout << "removing pieces restored unity. Continuing..\n";
                    no_cells_removed += get_no_cells_in_removed_pieces();
                }
                else {
                    // If removing the pieces failed to restore unity, undo the removal of the pieces and the cell
                    if (VERBOSE) cout << "Removal of cell " << cell
                        << " resulted in splitting of the shape, and removing pieces failed to restore unity.\n";
                    restore_removed_pieces(removed_pieces);
                    restore(cell);
                    if (VERBOSE) cout << "restored cell." << endl;
                    continue;
                }
            }
            if (no_cells_removed % (no_cells_to_remove / 10) == 0) cout << "no removed cells: "
                << no_cells_removed << " / " << no_cells_to_remove << endl;
            no_cells_removed++;
        }

        // If cell deletion leads to infeasibility, skip deletion
        int no_deleted_neighbors = 0;
        if (!cell_is_safe_to_delete(cell, no_deleted_neighbors)) {
            restore(cell);
            continue;
        }

        // Set cell to zero, making it empty
        if (smaller_piece == 0) remove_and_remember(cell);
        no_iterations_without_removal = 0;
        no_cells_removed = initial_count - count();
        if (no_cells_removed >= no_cells_to_remove) {
            break;
        }
    }
    return no_cells_removed;
}

int fessga::grd::Densities2d::get_no_cells_in_removed_pieces() {
    int no_cells = 0;
    for (auto& piece : removed_pieces) no_cells += piece.cells.size();
    return no_cells;
}

// Remove a floating piece (if possible)
bool fessga::grd::Densities2d::remove_floating_piece(
    grd::Piece* piece
) {
    vector<int> _removed_cells;
    for (auto& cell : piece->cells) {
        // TODO: The last two conditions of the following if-statement should not be necessary, but they are.
        // Figure out why.
        bool cell_cannot_be_removed = (
            fea_results.data_map[cell] > fea_casemanager->max_stress_threshold ||
            help::is_in(&fea_casemanager->cells_to_keep, cell)
        );
        if (cell_cannot_be_removed) {
            if (fea_casemanager->maintain_boundary_connection) {
                fill(_removed_cells);
                _count += removed_cells.size();
                piece->is_removable = false;
                return false;
            }
        }
        else {
            del(cell);
            _removed_cells.push_back(cell);
        }
    }
    move_piece_to_trash(piece);
    return true;
}

void fessga::grd::Densities2d::remove_smaller_pieces() {
    remove_smaller_pieces(pieces, &removed_cells);
}

// Remove all pieces smaller than the largest piece from the mesh, if possible
void fessga::grd::Densities2d::remove_smaller_pieces(
    vector<grd::Piece> _pieces, vector<int>* removed_cells, bool ignore_largest_piece, bool check_if_single_piece
) {
    int no_pieces = _pieces.size();
    int i = 1;
    for (auto& piece : _pieces) {
        // Ignore largest piece if this flag is set
        if (ignore_largest_piece && piece.id == main_piece.id) {
            i++; continue;
        }

        // Attempt to remove the floating piece
        bool piece_removed = false;
        if (piece.is_removable || !fea_casemanager->maintain_boundary_connection) {
            piece_removed = remove_floating_piece(&piece);
        }

        // Check if the piece was removed
        if (piece_removed) {
            if (VERBOSE) cout << "Removed piece " << i << " / " << no_pieces << endl;
        }
        else if (VERBOSE) cout << "Piece " << i << " was not removed." << endl;

        // If this flag is set, we must prevent that the removal of a piece results in the splitting of the shape
        // into multiple pieces.
        if (check_if_single_piece) {
            bool _is_single_piece = is_single_piece(-1, false);
            if (VERBOSE) cout << "CHECKING WHETHER PIECE REMOVAL RESULTED IN MULTIPLE PIECES....\n";

            // If the shape is broken into multiple pieces, undo the removal of the current piece
            if (!_is_single_piece) {
                if (VERBOSE) cout << "UNDOING THE REMOVAL OF PIECE " << i
                    << " BECAUSE ITS REMOVAL SPLIT THE SHAPE INTO MULTIPLE PIECES" << endl;
                vector<grd::Piece> removed_piece = { piece };
                restore_removed_pieces(removed_piece);
            }
        }

        i++;
    }
}

// Filter out floating cells that have no direct neighbors
void fessga::grd::Densities3d::filter() {
#pragma omp parallel for
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            for (int z = 0; z < dim_z; z++) {
                int filled = values[x * dim_z * dim_y + y * dim_z + z];
                if (!filled) continue;
                int neighbor = 0;
                for (int _x = -1; _x <= 1; _x++) {
                    for (int _y = -1; _y <= 1; _y++) {
                        for (int _z = -1; _z <= 1; _z++) {
                            if (x + _x == dim_x || y + _y == dim_y || z + _z == dim_z) continue;
                            if (x + _x <= 0 || y + _y <= 0 || z + _z <= 0) continue;
                            neighbor = values[(x + _x) * dim_z * dim_y + (y + _y) * dim_z + (z + _z)];
                            if (neighbor) break;
                        }
                        if (neighbor) break;
                    }
                    if (neighbor) break;
                }
                if (!neighbor) {
                    del(x * dim_z * dim_y + y * dim_z + z); // Set density to 0 if the cell has no neighbors
                }
            }
        }
    }
    cout << "Finished filtering floating cells." << endl;
}

void fessga::grd::Densities3d::fill_cells_inside_mesh(Vector3d offset, MatrixXd* V, MatrixXi* F) {
    // Compute the barycenter of the mesh
    Vector3d mesh_barycent = V->colwise().mean();

    // Compute vector to center of a grid cell from its corner
    Vector3d to_cell_center = Vector3d(0.5, 0.5, 0.5).cwiseProduct(cell_size);

    // Initialize different ray directions, (this is a temporary fix for a bug whereby some cells
    // are not properly assigned a density of 1)
    vector<Vector3d> ray_directions = {
        Vector3d(0, 1.0, 0), Vector3d(1.0, 1.0, 1.0).normalized(), Vector3d(1.0, 0, 0)
    };

#if 0:
    // Create list of triangles
    std::vector<tracer::Triangle> triangles;
    for (int face_idx = 0; face_idx < F->rows(); face_idx++) {
        tracer::Triangle triangle;
        triangle.v0 = V->row(F->coeff(face_idx, 0));
        triangle.v1 = V->row(F->coeff(face_idx, 1));
        triangle.v2 = V->row(F->coeff(face_idx, 2));
        triangles.push_back(triangle);
    }

    int slices_done = 0;
#pragma omp parallel for
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            for (int z = 0; z < dim_z; z++) {
                tracer::Cell3D cell;
                Vector3d indices; indices << x, y, z;
                cell.position = offset + indices.cwiseProduct(cell_size) + to_cell_center;
                cell.density = 0;
                // Try three ray directions (temporary bug fix, see 'Initialize two different ray directions' above) 
                for (int i = 0; i < 3; i++) {
                    // Cast ray and check for hits
                    tracer::Ray ray;
                    ray.origin = cell.position;
                    ray.direction = ray_directions[i];
                    Vector3d hitPoint;
                    Vector3d hit_normal;
                    bool hit = tracer::trace_ray(ray, triangles, hitPoint, hit_normal);

                    // If there was a hit, check if the hit triangle's normal points in the same direction as the ray
                    // If so, the cell must be inside the mesh
                    bool inside = false;
                    if (hit) {
                        inside = hit_normal.dot(ray.direction) > 0.0;
                    }

                    // If the cell is inside the mesh, assign density 1
                    if (inside) {
                        cell.density = 1;
                        break;
                    }
                }
                set(x * dim_z * dim_y + y * dim_z + z, cell.density);
            }
        }
        cout << "    Processed slice " << slices_done + 1 << " / " << dim_x << endl;
        slices_done++;
    }
#endif
}

/* Generate a binary density distribution on the grid based on the given mesh
*/
void fessga::grd::Densities3d::generate(Vector3d offset, MatrixXd* V, MatrixXi* F) {
    cout << "Generating 3d grid-based density distribution..." << endl;
    fill_cells_inside_mesh(offset, V, F);
    update_count();
    filter();
    cout << "Finished generating density distribution." << endl;
}

void fessga::grd::Densities3d::create_x_slice(grd::Densities2d& densities2d, int x) {
    for (int z = 0; z < dim_z; z++) {
        for (int y = 0; y < dim_y; y++) {
            densities2d.set(x * dim_y + y, values[z * dim_x * dim_y + x * dim_y + y]);
        }
    }
}

void fessga::grd::Densities3d::create_y_slice(grd::Densities2d& densities2d, int y) {
    for (int x = 0; x < dim_x; x++) {
        for (int z = 0; z < dim_z; z++) {
            densities2d.set(x * dim_y + y, values[z * dim_x * dim_y + x * dim_y + y]);
        }
    }
}

void fessga::grd::Densities3d::create_z_slice(grd::Densities2d& densities2d, int z) {
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            densities2d.set(x * dim_y + y, values[z * dim_x * dim_y + x * dim_y + y]);
        }
    }
}

void fessga::grd::Densities3d::create_slice(grd::Densities2d& densities2d, int dimension, int offset) {
    switch (dimension) {
        case 0: create_x_slice(densities2d, offset); return;
        case 1: create_y_slice(densities2d, offset); return;
        case 2: create_z_slice(densities2d, offset); return;
    }
    densities2d.update_count();
}

bool fessga::grd::Densities2d::repair() {
    do_feasibility_filtering();
    bool _is_single_piece = remove_isolated_material();
    bool is_valid = _is_single_piece;
    return is_valid;
}

bool fessga::grd::Densities2d::remove_isolated_material() {
    flush_edit_memory();
    init_pieces();
    if (pieces.size() == 1) return true;
    remove_smaller_pieces();
    return pieces.size() == 1;
}

void fessga::grd::Densities2d::fill_voids(int target_no_neighbors) {
    save_internal_snapshot();
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            int cell = x * dim_y + y;

            // If the cell is already filled, skip it
            if (values[cell]) continue;

            // Fill the cell if the number of neighbors is equal to the required number computed earlier.
            vector<int> neighbors = get_neighbors(x, y, snapshot_internal);
            if (neighbors.size() == target_no_neighbors) {
                fill(cell);
            }
        }
    }
    for (auto& cutout_cell : fea_casemanager->cutout_cells) {
        del(cutout_cell);
    }
}

// Use a 2x2 kernel to check for level1 voids (2x2 pockets of cells that are empty) and fill them.
void fessga::grd::Densities2d::fix_level1_void_or_pinch(vector<pair<int, int>>* offsets, int x, int y) {
    // Check if the kernel contains a void
    bool is_level1_void = true;
    bool is_pinch = false;
    for (auto& offset : *offsets) {
        int cell = get_idx(pair(x + offset.first, y + offset.second));

        // If the cell is not empty or has fewer than 2 neighbors, the kernel apparently does not contain
        // a 2x2 void (voids that are not strictly 2x2 are not considered level1 voids)
        if (values[cell] || get_neighbors(cell).size() < 2) {
            is_level1_void = false;

            // Check for the presence of a pinch. A pinch is defined as a node that forms the only connection between
            // two cells which neighbor each other diagonally. Therefore, check for the presence of two diagonally
            // neighboring filled cells and two empty cells.
            int bottomleft = at(x, y); int topright = at(x + 1, y + 1);
            int topleft = at(x, y + 1); int bottomright = at(x + 1, y);
            bool diagonal_1_full = bottomleft && topright;
            bool diagonal_2_full = bottomright && topleft;
            bool diagonal_1_empty = !bottomleft && !topright;
            bool diagonal_2_empty = !bottomright && !topleft;
            is_pinch = diagonal_1_full && diagonal_2_empty; // 1st diagonal is filled, 2nd empty
            is_pinch = is_pinch || (diagonal_1_empty && diagonal_2_full); // Same but reverse
            break;
        }
    }

    // If the kernel contains a void or a pinch, fill it
    if (is_level1_void || is_pinch) {
        for (auto& offset : *offsets) {
            int cell = get_idx(pair(x + offset.first, y + offset.second));
            fill(cell);
        }
    }
}

void fessga::grd::Densities2d::fill_level1_voids_and_fix_pinches(bool verbose) {
    vector<pair<int, int>> offsets = { pair(0,0), pair(0,1), pair(1,1), pair(1,0) };
    for (int x = 0; x < dim_x - 1; x++) {
        for (int y = 0; y < dim_y - 1; y++) {
            fix_level1_void_or_pinch(&offsets, x, y);
        }
    }
    // If void and/or pinch filling filled some cutout cells, delete these.
    for (auto& cutout_cell : fea_casemanager->cutout_cells) {
        del(cutout_cell);
    }
}

void fessga::grd::Densities2d::do_feasibility_filtering(bool verbose) {
    grd::Densities2d previous_state(this);
    bool filtering_had_effect = true;
    int i = 1;

    // Run level0 (meaning 'acting on individual cells') filtering loop 
    while (filtering_had_effect) {
        do_single_feasibility_filtering_pass();
        filtering_had_effect = !previous_state.is_identical_to(values);
        previous_state.copy_from(this);
        i++;
    }

    // Fill level1 voids (2x2 pockets of cells that are empty)
    fill_level1_voids_and_fix_pinches(verbose);
    
    if (verbose) cout << "Performed feasibility filtering (" << i << " passes).\n";
}

void fessga::grd::Densities2d::do_single_feasibility_filtering_pass() {
    // Step 1: Fill void elements that are surrounded on all sides by solid elements
    fill_voids(4);

    // Step 2: Remove solid elements that have 0 true neighbors
    filter(0, true);

    // Step 3: Fill void elements that are surrounded on all but one side by solid elements
    fill_voids(3);

    // Step 4: Remove solid elements that have exactly 1 true neighbor
    filter(1, true);
}

// Get empty neighbor that shares the cell's given line
int fessga::grd::Densities2d::get_empty_neighbor_cell_of_line(int cell, int local_line_idx) {
    pair<int, int> cell_coords = get_coords(cell);
    pair<int, int> neighbor;
    if (local_line_idx == 0) neighbor = pair(cell_coords.first - 1, cell_coords.second);
    else if (local_line_idx == 1) neighbor = pair(cell_coords.first, cell_coords.second + 1);
    else if (local_line_idx == 2) neighbor = pair(cell_coords.first + 1, cell_coords.second);
    else if (local_line_idx == 3) neighbor = pair(cell_coords.first, cell_coords.second - 1);
    
    // Check if neighbor is inside the design domain
    if (neighbor.first > dim_x || neighbor.first < 0 || neighbor.second > dim_y || neighbor.second < 0) return -1;

    return get_idx(neighbor);
}

// Save a copy of the current state of the density distribution
void fessga::grd::Densities2d::save_snapshot() {
    for (int i = 0; i < size; i++) snapshot[i] = values[i];
    _snapshot_count = _count;
    fea_results_snapshot = fea_results;
}

// Save a copy of the current state of the density distribution, to be used only internally
void fessga::grd::Densities2d::save_internal_snapshot() {
    for (int i = 0; i < size; i++) snapshot_internal[i] = values[i];
    _snapshot_internal_count = _count;
}

// Load the previously saved state (i.e. the snapshot)
void fessga::grd::Densities2d::load_snapshot() {
    for (int i = 0; i < size; i++) values[i] = snapshot[i];
    _count = _snapshot_count;
    fea_results = fea_results_snapshot;
}

// Load the previously saved state (i.e. the snapshot). Only for internal use.
void fessga::grd::Densities2d::load_internal_snapshot() {
    for (int i = 0; i < size; i++) values[i] = snapshot_internal[i];
    _count = _snapshot_internal_count;
}

// Copy the density values from the given Densities2d-object to the current object
void fessga::grd::Densities2d::copy_from(Densities2d* source) {
    assert(source->dim_x == dim_x && source->dim_y == dim_y);
    for (int i = 0; i < size; i++) values[i] = source->at(i);
    _count = source->count();
}

// Copy the density values from one array to another
void fessga::grd::Densities2d::_copy(uint* source, uint* target, int source_count, int target_count) {
    for (int i = 0; i < size; i++) values[i] = source[i];
    target_count = source_count;
}

// Compute Center Of Mass
void fessga::grd::Densities2d::compute_center_of_mass(bool verbose) {
    float com_x = 0; float com_y = 0;
    for (int i = 0; i < size; i++) {
        int x = i / dim_y;
        int y = i % dim_y;
        if (values[i] == 1) {
            com_x += (float)x * cell_size[0] + cell_size[0] * 0.5;
            com_y += (float)y * cell_size[1] + cell_size[1] * 0.5;
        }
    }
    com_x /= (float)count();
    com_y /= (float)count();
    if (verbose) cout << "center of gravity: " << com_x << ", " << com_y << endl;
    int cent_coord = (int)(com_x / cell_size[0]) * dim_y + com_y / cell_size[1];
}

void fessga::grd::Densities2d::compute_area(bool verbose) {
    area = (cell_size[0] * cell_size[1]) * (double)count();
}

void fessga::grd::Densities2d::invert() {
    for (int i = 0; i < count(); i++) values[i] = !values[i];
}
