#include <iostream>
#include <filesystem>
#include "io.h"
#include "helpers.h"
#include "meshing.h"

using namespace Eigen;
using namespace fessga;
using namespace std;


/* 
Generate a 2d Finite Element mesh from the given binary density distribution
*/
void mesher::generate_FE_mesh(
    Grid3D grid, SurfaceMesh mesh, uint* densities, FEMesh2D& fe_mesh
) {
    // Create nodes and surfaces
    vector<vector<double>> nodes;
    vector<Element> surfaces;
    int node_idx = 1;
    std::map<int, int> node_coords;
    for (int x = 0; x < grid.x; x++) {
        for (int y = 0; y < grid.y; y++) {
            int filled = densities[x * grid.y + y];
            if (filled) {
                // Create 4 nodes that will enclose the surface element to be created (if they do not yet exist)
                int node4_idx = mesher::add_node_if_not_exists(x, y, mesh.offset, grid.y, node_coords, grid.cell_size, nodes, node_idx);
                int node1_idx = mesher::add_node_if_not_exists(x + 1, y, mesh.offset, grid.y, node_coords, grid.cell_size, nodes, node_idx);
                int node2_idx = mesher::add_node_if_not_exists(x + 1, y + 1, mesh.offset, grid.y, node_coords, grid.cell_size, nodes, node_idx);
                int node3_idx = mesher::add_node_if_not_exists(x, y + 1, mesh.offset, grid.y, node_coords, grid.cell_size, nodes, node_idx);

                // Compute surface index (equals flattened coordinates + 1)
                int surface_idx = x * grid.y + y + 1;

                // For the 2d case, surfaces are never boundary elements. Therefore, tag is given a default value of 1.
                int tag = 1;

                // Other default indices
                int type = 3; // quad
                int physical_entity = 0;
                int no_tags = 2;

                // Create the surface element
                Element surface;
                surface.id = surface_idx;
                surface.type = type;
                surface.no_tags = no_tags;
                surface.body = physical_entity;
                surface.tag = tag;
                surface.nodes = { node1_idx, node2_idx, node3_idx, node4_idx };
                surfaces.push_back(surface);
            }
        }
    }

    // Find boundary nodes
    std::map<int, int> boundary_node_coords;
    int no_boundary_nodes = 0;
    for (auto const& [node_coord, node_idx] : node_coords) {
        int x = node_coord / (grid.y + 1);
        int y = node_coord % (grid.y + 1);

        // Nodes at an extreme of at least one axis are always boundary nodes
        if (x == 0 || y == 0 || x == grid.x || y == grid.y) {
            boundary_node_coords[node_coord] = node_idx;
            no_boundary_nodes++;
            continue;
        }

        // Nodes with at least one empty neighboring cell are boundary nodes
        bool empty = false;
        for (int x_offset = 0; x_offset < 2; x_offset++) {
            for (int y_offset = 0; y_offset < 2; y_offset++) {
                int cell_is_filled = densities[(x - 1 + x_offset) * grid.y + (y - 1 + y_offset)];
                if (!cell_is_filled) {
                    boundary_node_coords[node_coord] = node_idx;
                    empty = true;
                    no_boundary_nodes++;
                    break;
                }
            }
            if (empty) break;
        }
    }
    //cout << "unordered boundary node coords: \n";
    for (auto [coord, idx] : boundary_node_coords) {
        int x = coord / (grid.y + 1);
        int y = coord % (grid.y + 1);
        //cout << "(" << x << ", " << y << "), ";
    }
    //cout << endl;

    //cout << "node coord of node 6: " << fessga::help::get_key(&boundary_node_coords, 6) << endl;

    // ----  Order boundary nodes according to occurrence along perimeter of the mesh ---- //
    // Initialize starting node and ordered boundary nodes maps
    int node_coord = boundary_node_coords.begin()->first;
    node_idx = boundary_node_coords.begin()->second;
    std::vector<int> ordered_boundary_node_coords = { node_coord };
    int x = node_coord / (grid.y + 1);
    int y = node_coord % (grid.y + 1);
    int start_x = x;
    int start_y = y;

    // Initialize 'previous' coordinates to an arbitrary neighboring location on the grid 
    int previous_x = x;
    int previous_y = y - 1;
    if (previous_y < 0) previous_y = y + 1;

    // Trace perimeter by stepping from one boundary node to a direct neighbor that
    // was not visited in the previous iteration
    int i = 1;
    bool neighbor_found = false;
    int neighbor_idx, neighbor_coord, neighbor_x, neighbor_y;
    int _neighbor_idx, _neighbor_coord, _neighbor_x, _neighbor_y;
    vector<pair<int, int>> offsets = { pair(0,1), pair(1,0), pair(-1, 0), pair(0, -1) };
    int no_components = 1;
    bool force_new_component_search = false;
    while (true) {
        if (i > (no_boundary_nodes + no_components)) {
            node_coord = mesher::find_unvisited_node(&boundary_node_coords, &ordered_boundary_node_coords);
            if (node_coord == -1) {
                break; // If there's no unvisited nodes, break the while-loop since we're done.
            }
            else {
                force_new_component_search = true;
            }
        }

        // If current node is the same as the starting node of the perimeter walk, re-start walk on a
        // node that has not yet been visited (such a node must be part of another component, for example a hole)
        if (force_new_component_search || (i > 1 && x == start_x && y == start_y)) {
            force_new_component_search = false;
            node_coord = mesher::find_unvisited_node(&boundary_node_coords, &ordered_boundary_node_coords);
            if (node_coord == -1) {
                break; // If there's no unvisited nodes, break the while-loop since we're done.
            }
            else {
                no_components++;
                ordered_boundary_node_coords.push_back(-1); // Add a separator between components, to avoid the generation of lines between them
            }
            int node_idx = boundary_node_coords[node_coord];
            ordered_boundary_node_coords.push_back(node_coord);
            x = node_coord / (grid.y + 1);
            y = node_coord % (grid.y + 1);
            start_x = x;
            start_y = y;

            // Re-initialize 'previous' coordinates to an arbitrary neighboring location on the grid
            previous_x = x;
            previous_y = y - 1;
            if (previous_y < 0) previous_y = y + 1;
            
            i++;
        }

        // Check which of the other three neighboring grid indices contains a boundary node
        vector<pair<int, int>> valid_neighbors;
        for (auto _offset : offsets) {
            // Compute neighbor coordinates
            _neighbor_x = (x + _offset.first);
            _neighbor_y = (y + _offset.second);

            // Check coordinate validity
            if (_neighbor_x < 0 || _neighbor_y < 0 || _neighbor_x > grid.x || _neighbor_y > grid.y)
                continue; // coordinates outside design domain are invalid
            if (_neighbor_x == previous_x && _neighbor_y == previous_y)
                continue; // skip the boundary node we found in the previous iteration

            // Check whether line connecting current node to neighbor is blocked by two adjacent filled cells
            bool blocked = false;
            bool _blocked = false;
            bool infeasible = false;
            if (x == _neighbor_x) // Vertical line
            {
                int y_diff = min(0, _neighbor_y - y);
                int right_cell, left_cell;
                if (x < grid.x) right_cell = densities[x * grid.y + y + y_diff];
                if (x > 0) left_cell = densities[(x - 1) * grid.y + y + y_diff];
                if (x == 0) {
                    infeasible = !right_cell; // Infeasible if the line is on the left x-limit and there's no cell to the right
                }
                else if (x == grid.x) {
                    infeasible = !left_cell; // Infeasible if the line is on the right x-limit and there's no cell to the left
                }
                else {
                    _blocked = right_cell && left_cell;
                    infeasible = !right_cell && !left_cell; // Infeasible if there's no cells at all (since then an edge would be drawn in thin air)
                }
            }
            else if (y == _neighbor_y) // Horizontal line
            {
                int x_diff = min(0, _neighbor_x - x);
                int up_cell, down_cell;
                if (y < grid.y) up_cell = densities[(x + x_diff) * grid.y + y];
                if (y > 0) down_cell = densities[(x + x_diff) * grid.y + y - 1];
                if (y == 0) {
                    infeasible = !up_cell; // Infeasible if the line is on the lower y-limit and there's no cell above it
                    //cout << "coming from (" << x << ", " << y << ") a path to " << grid.y + y << " is infeasible." << endl;
                }
                else if (y == grid.y) {
                    infeasible = !down_cell; // Infeasible if the line is on the upper y-limit and there's no cell below it
                }
                else {
                    _blocked = up_cell && down_cell;
                    infeasible = !up_cell && !down_cell; // Infeasible if there's no cells at all (since then an edge would be drawn in thin air)
                }
            }
            blocked = _blocked || infeasible;
            if (blocked) continue; // Skip neighbor if line is blocked by adjacent filled cells

            // Check whether the neighbor coordinates contain a boundary node
            _neighbor_coord = _neighbor_x * (grid.y + 1) + _neighbor_y;
            _neighbor_idx = fessga::help::get_value(&boundary_node_coords, _neighbor_coord);
            neighbor_found = _neighbor_idx != -1;
            if (neighbor_found) {
                neighbor_x = _neighbor_x;
                neighbor_y = _neighbor_y;
                neighbor_idx = _neighbor_idx;
                neighbor_coord = _neighbor_coord;
                valid_neighbors.push_back(pair(_neighbor_coord, _neighbor_idx));
            }
            else {
            }
        }

        // Check whether more than one valid neighbor was found. If so, choose a neighbor that has not been visited yet
        if (valid_neighbors.size() > 1) {
            int j = 0;
            while (j < valid_neighbors.size()) {
                //cout << "checking whether " << valid_neighbors[j].first << " is in ordered boundary node coords" << endl;
                bool neighbor_visited = fessga::help::is_in(&ordered_boundary_node_coords, valid_neighbors[j].first);
                if (!neighbor_visited) {
                    break; // Found a valid neighbor that hasn't been visited yet
                }
                j++;
            }
            if (j == valid_neighbors.size()) {
                // No unvisited neighbor was found. This must mean that we have reached the last line of the current perimeter component.
                help::print("adding node (" + to_string(start_x) + ", " + to_string(start_y) + ")\n");
                ordered_boundary_node_coords.push_back(start_x * (grid.y + 1) + start_y);
                force_new_component_search = true;
                continue;
            }
            else {
                // Set the neighbor data to the data of the unvisited neighbor
                neighbor_coord = valid_neighbors[j].first;
                neighbor_idx = valid_neighbors[j].second;
                neighbor_x = neighbor_coord / (grid.y + 1);
                neighbor_y = neighbor_coord % (grid.y + 1);
                //cout << "Choosing unvisited node " << neighbor_x << ", " << neighbor_y << endl;
            }

        }

        if (neighbor_idx == -1) {
            cout << "ERROR: Neighbor not found. Previous x, y: " << previous_x << ", " << previous_y <<
                ", next x, y: " << neighbor_x << ", " << neighbor_y << endl;
        }

        // Add neighboring boundary node to vector
        help::print("adding node (" + to_string(neighbor_x) + ", " + to_string(neighbor_y) + ")\n");
        ordered_boundary_node_coords.push_back(neighbor_coord);

        // Set previous x and y coordinates to current ones
        previous_x = x;
        previous_y = y;

        // Set next node coordinates to those of neighbor
        x = neighbor_x;
        y = neighbor_y;

        i++;
    }

    // Hotfix: In case the start node of the last-traced component was somehow not added to the ordered boundary node coords, add it.
    if (ordered_boundary_node_coords.back() != (start_x * (grid.y + 1) + start_y))
        ordered_boundary_node_coords.push_back((start_x * (grid.y + 1) + start_y));

    cout << "   no perimeter components: " << no_components << endl;
    cout << "   no unordered bound nodes: " << boundary_node_coords.size() << endl;
    cout << "   no ordered bound nodes: " << ordered_boundary_node_coords.size() << endl;
    //cout << "ordered bound nodes: "; fessga::help::print_vector(&ordered_boundary_node_coords);

    // ---- Generate boundary lines ---- //
    int tag = 1;
    int no_lines = 0;
    vector<Element> lines;
    for (int i = 1; i < ordered_boundary_node_coords.size(); i++) {
        // Get the 2 nodes that define the line
        int node1_coord = ordered_boundary_node_coords[i - 1];
        int node2_coord = ordered_boundary_node_coords[i];

        // Skip node pairs containing a -1. These indicate a transition from one component to another. No lines should be generated here.
        if (node1_coord == -1 || node2_coord == -1) continue;

        int node1_idx = boundary_node_coords[node1_coord];
        int node2_idx = boundary_node_coords[node2_coord];

        // Get line idx. Convention: line_idx = (<surface_coord> << 2) + <local_line_idx>
        // local_line_idx: Stepping clockwise from left edge of cell (local_line_idx = 0)
        int node1_x = node1_coord / (grid.y + 1);
        int node1_y = node1_coord % (grid.y + 1);
        int node2_x = node2_coord / (grid.y + 1);
        int node2_y = node2_coord % (grid.y + 1);
        int cell_x, cell_y;
        uint cell_coord, local_line_idx;
        if (node1_x == node2_x) {   // Vertical line
            cell_x = node1_x;
            cell_y = min(node1_y, node2_y);
            local_line_idx = 0;
            // Choose cell that has the line as its left side unless line is on the rightmost x-limit or the right cell is empty
            if (node1_x == grid.x || !densities[cell_x * grid.y + cell_y]) {
                cell_x = node1_x - 1;
                local_line_idx = 2;
            }
        }
        else { // Horizontal line
            cell_y = node1_y;
            cell_x = min(node1_x, node2_x);
            local_line_idx = 3;
            // Choose cell that has the line as its bottom side unless line is on the upper y-limit or the upper cell is empty
            if (node1_y == grid.y || !densities[cell_x * grid.y + cell_y]) { 
                cell_y = node1_y - 1;
                local_line_idx = 1;
            }
        }
        // Use the last two bits of the integer to store the line idx. The first 30 are used for the cell coordinates
        cell_coord = cell_x * grid.y + cell_y;
        uint line_idx = (cell_coord << 2) + local_line_idx;

        // Increment the tag
        tag++;

        // Unchanging indices
        int type = 1;
        int no_tags = 2;
        int physical_entity = 0;

        // Create the line element
        Element line;
        line.id = line_idx + 1;
        line.type = type;
        line.no_tags = no_tags;
        line.body = physical_entity;
        line.tag = tag;
        line.nodes = { node1_coord + 1, node2_coord + 1 };
        line.boundary_id = i; // Add i to boundary id, creating a continuous path around the perimeter(s) for Elmer's boundary element parser to walk
        lines.push_back(line);

        no_lines++;
    }

    // Populate FE mesh member variables with generated data
    fe_mesh.lines = lines;
    fe_mesh.surfaces = surfaces;
    fe_mesh.nodes = nodes;
}