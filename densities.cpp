#include "densities.h"
#include "raytracing.h"


#define VERBOSE false


// Remove floating cells (i.e. cells that have no direct neighbors)
void fessga::grd::Densities2d::filter(int no_neighbors) {
    save_internal_snapshot();
    vector<int> x_bounds = { 0, dim_x - 1 };
    vector<int> y_bounds = { 0, dim_y - 1 };
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            int cell = x * dim_y + y;

            // If the cell is already filled, continue to the next cell
            if (values[cell]) continue;

            // Delete the cell if the number of neighbors is equal to the given number
            vector<int> neighbors = get_true_neighbors(x, y, snapshot_internal);
            if (neighbors.size() == no_neighbors) {
                del(cell);
            }
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

/*
* Export 3d density distribution to file
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
void fessga::grd::Densities2d::do_import(string path, Vector3d diagonal) {
    // Get a vector of strings representing the lines in the file
    vector<string> lines;
    IO::read_file_content(path, lines);

    // Get the number of dimensions of the distribution
    int _no_dimensions = stoi(lines[0]);
    int _dim_x = stoi(lines[1]), _dim_y = stoi(lines[2]);
    if (_no_dimensions == 3) {
        if (no_dimensions == 2) {
            cerr << "Error: Cannot load a 3d density distribution using a Densities2d object. Use Densities3d instead.\n";
            return;
        }
        int _dim_z = stoi(lines[3]);
        construct_grid(_dim_x, _dim_y, _dim_z);
    }
    else construct_grid(_dim_x, _dim_y);

    // Fill the densities array with the binary values stored in the last line of the file
    string densities_line = lines[_no_dimensions + 1];
    for (int i = 0; i < size; i++) {
        set(i, densities_line[i] - '0');
    }
    update_count();

    // Initialize cell size
    compute_cellsize(diagonal);
}


// Return the indices of the 'true neighbors' of the cell at the given coordinates.
// True neighbors are here defined as filled neighbor cells that share a line with the given cell
vector<int> fessga::grd::Densities2d::get_true_neighbors(int x, int y, uint* _values) {
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
vector<int> fessga::grd::Densities2d::get_true_neighbors(int idx) {
    pair<int, int> coords = get_coords(idx);
    return get_true_neighbors(coords.first, coords.second);
}

// Get number of connected cells of the given cell using a version of floodfill
int fessga::grd::Densities2d::get_no_connected_cells(int cell_coord, Piece& piece, phys::FEACase* fea_case, bool verbose) {
    piece.cells = { cell_coord };
    int i = 0;
    while (i < piece.cells.size()) {
        vector<int> neighbors = get_true_neighbors(piece.cells[i]);
        for (int j = 0; j < neighbors.size(); j++) {
            if (!help::is_in(&piece.cells, neighbors[j])) {
                piece.cells.push_back(neighbors[j]);
                if (fea_case != 0 && piece.is_removable) {
                    // If the piece was flagged as removable but one of its cells is a boundary cell, flag it as non-removable.
                    if (help::is_in(&fea_case->boundary_cells, neighbors[j])) {
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
    int cell_coord, vector<int>* removed_cells, int& no_deleted_neighbors, phys::FEACase* fea_case
) {
    vector<int> neighbors = get_true_neighbors(cell_coord);
    for (auto& neighbor : neighbors) {
        vector<int> sub_neighbors = get_true_neighbors(neighbor);

        // If the neighboring cell has only one true neighbor itself, deleting the current cell would make it float in mid-air and thus invalid.
        // Therefore we either delete the neighboring cell too, or - in case the neighbor is a bound condition cell - skip deletion alltogether.
        if (sub_neighbors.size() <= 1) {
            // If the cell has a line on which a boundary condition was applied, skip deletion
            if (help::is_in(&fea_case->boundary_cells, neighbor) || help::is_in(&fea_case->whitelisted_cells, neighbor)) {
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
        vector<int> neighbors = get_true_neighbors(removed_cell);
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
            vector<int> neighbors = get_true_neighbors(removed_cells.at(i));
            if (neighbors.size() > 0) { start_cell = neighbors[0]; break; }
        }
    }
    else start_cell = fea_case->boundary_cells[0]; // If no removed cells were stored, pick a cell on which a boundary condition was applied
    init_pieces(&visited_cells, cells_left, start_cell);
}

// Get the pieces inside the density distribution (protected internal function)
void fessga::grd::Densities2d::init_pieces(vector<int>* visited_cells, int cells_left, int start_cell) {
    Piece piece = Piece();
    int piece_size = get_no_connected_cells(start_cell, piece, fea_case);
    pieces.push_back(piece);

    // Check if the shape consists of one piece or several
    if (piece_size < cells_left) {
        // If the piece has boundary cells and it is the largest piece, it is considered the main piece
        if (piece.cells.size() > main_piece.cells.size() && help::have_overlap(&piece.cells, &fea_case->boundary_cells)) {
            set_main_piece(&piece);
        }
        cells_left -= piece_size;
        for (int i = 0; i < piece_size; i++) visited_cells->push_back(piece.cells[i]);

        // Recurse if there are still unvisited cells left
        if (cells_left > 0) {
            int unvisited_cell = get_unvisited_cell(visited_cells);
            //cout << "unvisited cell (" << unvisited_cell << ")\n";
            /*set(unvisited_cell, 5);
            print();*/
            if (unvisited_cell == -1) throw(
                "Error: No unvisited cell found, but there are still cells left that haven't been visited (" + to_string(cells_left) + ")\n");
            init_pieces(visited_cells, cells_left, unvisited_cell);
        }
    }
}

// Visualize distribution and highlight bound cells
void fessga::grd::Densities2d::visualize_bound_cells() {
    for (auto& cell : fea_case->boundary_cells) values[cell] ? set(cell, 5) : throw("Error: Boundary cell not filled.\n");
    print();
    for (auto& cell : fea_case->boundary_cells) set(cell, 1);
}

bool fessga::grd::Densities2d::is_in(vector<grd::Piece>* _pieces, grd::Piece* piece) {
    for (int i = 0; i < _pieces->size(); i++) {
        if (_pieces->at(i).id == piece->id) return true;
    }
    return false;
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
    else start_cell = fea_case->boundary_cells[0];
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
    fea_case->whitelisted_cells.clear();
    removed_pieces.clear();
    main_piece = grd::Piece();
}

int fessga::grd::Densities2d::remove_low_stress_cells(int no_cells_to_remove, int no_cells_removed, grd::Piece* smaller_piece) {
    int cell_from_smaller_piece;
    int no_iterations_without_removal = -1;
    for (auto& item : fea_results->data) {
        no_iterations_without_removal++;
        if (smaller_piece && no_iterations_without_removal > 200) {
            if (true) cout << "WARNING: 200 cells in a row could not be removed. There may be no more cells to remove.\n";
            return removed_cells.size();
        }
        int cell = item.first;
        double cell_stress = item.second;

        // If the cell's stress exceeds the maximum, break the loop (all subsequent cells will also exceed the maximum since the list is ordered).
        if (cell_stress > fea_case->max_stress_threshold) continue;

        if (smaller_piece && no_iterations_without_removal > 100) cout << "before boundcells\n";
        
        // If the cell has a line on which a boundary condition was applied, skip deletion
        if (help::is_in(&fea_case->boundary_cells, cell)) {
            continue;
        }

        if (smaller_piece && no_iterations_without_removal > 100) cout << "before value check\n";

        // If the cell was already empty, skip deletion (more importantly: don't count this as a deletion)
        if (!values[cell]) continue;

        if (smaller_piece && no_iterations_without_removal > 100) cout << "before removal\n";

        // If a 'smaller piece' vector was provided, perform cell removal in 'careful mode'. This means: check whether cell deletion results in 
        // multiple pieces. If so, undo the deletion and whitelist the cell.
        if (smaller_piece != 0) {
            vector<int> _removed_cell = { cell };
            remove_and_remember(cell);
            cell_from_smaller_piece = smaller_piece->cells[0];
            int total_no_cells = count();
            bool _is_single_piece = is_single_piece(cell_from_smaller_piece, VERBOSE);
            if (!_is_single_piece) {

                // Try removing the smaller pieces
                pieces.clear(); removed_pieces.clear(); main_piece = grd::Piece(); // Flush edit memory except for removed_cells
                init_pieces();
                vector<int> removed_cell = { cell };
                remove_smaller_pieces(pieces, &removed_cell, true);
                if (pieces.size() == 1) {
                    //cout << "removing pieces restored unity. Continuing..\n";
                    no_cells_removed += get_no_cells_in_removed_pieces();
                }
                else {
                    // If removing the pieces failed to restore unity, undo the removal of the pieces and the cell
                    if (VERBOSE) cout << "Removal of cell " << cell << " resulted in splitting of the shape, and removing pieces failed to restore unity.\n";
                    restore_removed_pieces(removed_pieces);
                    //if (no_iterations_without_removal > 100) {
                    //    set(cell, 5); // TEMP
                    //    print(); // TEMP
                    //}
                    restore(cell);
                    if (VERBOSE) cout << "restored cell." << endl;
                    continue;
                }
            }
            if (no_cells_removed % (no_cells_to_remove / 10) == 0) cout << "no removed cells: " << no_cells_removed << " / " << no_cells_to_remove << endl;
            no_cells_removed++;
        }

        // If cell deletion leads to infeasibility, skip deletion
        int no_deleted_neighbors = 0;
        if (!cell_is_safe_to_delete(cell, &removed_cells, no_deleted_neighbors, fea_case)) {
            restore(cell);
            //cout << "cell is not safe to delete" << endl;
            continue;
        }

        // Set cell to zero, making it empty
        //cout << "no removed cells: " << removed_cells.size() << endl;
        if (smaller_piece == 0) remove_and_remember(cell);
        no_iterations_without_removal = 0;
        if (removed_cells.size() >= no_cells_to_remove) break;
    }
    return removed_cells.size();
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
        // TODO: The last two conditions of the following if-statement should not be necessary, but they are. Figure out why.
        bool cell_cannot_be_removed = fea_results->data_map[cell] > fea_case->max_stress_threshold ||
            (help::is_in(&fea_case->boundary_cells, cell) || help::is_in(&fea_case->whitelisted_cells, cell)
                );
        if (cell_cannot_be_removed) {
            if (fea_case->maintain_boundary_connection) {
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
        if (piece.is_removable || !fea_case->maintain_boundary_connection) {
            piece_removed = remove_floating_piece(&piece);
        }

        // Check if the piece was removed
        if (piece_removed) {
            if (VERBOSE) cout << "Removed piece " << i << " / " << no_pieces << endl;
        }
        else if (VERBOSE) cout << "Piece " << i << " was not removed." << endl;

        // If this flag is set, we must prevent that the removal of a piece results in the splitting of the shape into multiple pieces.
        if (check_if_single_piece) {
            bool _is_single_piece = is_single_piece(-1, false);
            if (VERBOSE) cout << "CHECKING WHETHER PIECE REMOVAL RESULTED IN MULTIPLE PIECES....\n";

            // If the shape is broken into multiple pieces, undo the removal of the current piece
            if (!_is_single_piece) {
                if (VERBOSE) cout << "UNDOING THE REMOVAL OF PIECE " << i << " BECAUSE ITS REMOVAL SPLIT THE SHAPE INTO MULTIPLE PIECES" << endl;
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
    vector<Vector3d> ray_directions = { Vector3d(0, 1.0, 0), Vector3d(1.0, 1.0, 1.0).normalized(), Vector3d(1.0, 0, 0) };

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

// Save a copy of the current state of the density distribution
void fessga::grd::Densities2d::save_snapshot() {
    for (int i = 0; i < size; i++) snapshot[i] = values[i];
    _snapshot_count = _count;
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
}

// Load the previously saved state (i.e. the snapshot). Only for internal use.
void fessga::grd::Densities2d::load_internal_snapshot() {
    for (int i = 0; i < size; i++) values[i] = snapshot_internal[i];
    _count = _snapshot_internal_count;
}

// Copy the density values from the given Densities2d-object to the current object
void fessga::grd::Densities2d::copy_from(Densities2d* source) {
    for (int i = 0; i < size; i++) values[i] = source->at(i);
    _count = source->count();
}

// Copy the density values from one array to another
void fessga::grd::Densities2d::_copy(uint* source, uint* target, int source_count, int target_count) {
    for (int i = 0; i < size; i++) values[i] = source[i];
    target_count = source_count;
}
