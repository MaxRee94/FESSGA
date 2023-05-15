#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <map>
#include "Helpers.h"


using namespace Eigen;
using namespace std;
using namespace mvis;


// Define the struct for a cell in the grid
struct Cell {
    Vector3d position;
    int density;
};

void print_density_distrib(uint32_t* densities, int grid_size) {
    int x = grid_size / 2;
    for (int y = 0; y < grid_size; y++) {
        for (int z = 0; z < grid_size; z++) {
            cout << densities[x * grid_size * grid_size + y * grid_size + z];
        }
        cout << endl;
    }
}


struct Ray {
    Vector3d origin;
    Vector3d direction;
};

struct Triangle {
    Vector3d v0 = Vector3d(0.0, 0.0, 0.0);
    Vector3d v1 = Vector3d(0.0, 0.0, 0.0);
    Vector3d v2 = Vector3d(0.0, 0.0, 0.0);
};

bool trace_ray(const Ray& ray, const std::vector<Triangle>& triangles, Vector3d& hitPoint, Vector3d& hit_normal) {
    double closestDist = INFINITY;
    bool hit = false;

    for (const Triangle& triangle : triangles) {
        Vector3d e1 = triangle.v1 - triangle.v0;
        Vector3d e2 = triangle.v2 - triangle.v0;
        Vector3d h = ray.direction.cross(e2);
        double a = e1.dot(h);

        if (a > -1e-8 && a < 1e-8) {
            continue;
        }

        double f = 1.0 / a;
        Vector3d s = ray.origin - triangle.v0;
        double u = f * s.dot(h);

        if (u < 0 || u > 1) {
            continue;
        }

        Vector3d q = s.cross(e1);
        double v = f * ray.direction.dot(q);

        if (v < 0 || u + v > 1) {
            continue;
        }

        double dist = f * e2.dot(q);

        if (dist > 1e-12 && dist < closestDist) {
            closestDist = dist;
            hit = true;
            hitPoint = ray.origin + ray.direction * dist;
            hit_normal = e1.cross(e2).normalized();
        }
    }

    return hit;
}


/* Generate a binary density distribution on the grid based on the given mesh file
Input:
Returns (by reference):
    densities (uint32_t*): Array which contains a binary density value (0 or 1) for each cell in the grid
*/
void generate_3d_density_distribution(
    int dim_x, int dim_y, int dim_z, Vector3d offset, double cell_size, MatrixXd* V, MatrixXi* F,
    uint32_t* densities
) {
    cout << "Generating 3d grid-based density distribution..." << endl;

    // Compute the barycenter of the mesh
    Vector3d mesh_barycent = V->colwise().mean();

    // Compute the size of each cell in the grid
    Vector3d grid_size = Vector3d(dim_x, dim_y, dim_z);

    // Compute vector to center of a grid cell from its corner
    Vector3d to_cell_center = Vector3d(0.5, 0.5, 0.5) * cell_size;

    // Initialize two different ray directions, (this is a temporary fix for a bug whereby some cells
    // are not properly assigned a density of 1)
    vector<Vector3d> ray_directions = { Vector3d(0, 1.0, 0), Vector3d(1.0, 1.0, 1.0).normalized() };

    // Create list of triangles
    std::vector<Triangle> triangles;
    for (int face_idx = 0; face_idx < F->rows(); face_idx++) {
        Triangle triangle;
        triangle.v0 = V->row(F->coeff(face_idx, 0));
        triangle.v1 = V->row(F->coeff(face_idx, 1));
        triangle.v2 = V->row(F->coeff(face_idx, 2));
        triangles.push_back(triangle);
    }

    // Assign density values to cells in the grid
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            for (int z = 0; z < dim_z; z++) {
                Cell cell;
                Vector3d indices; indices << x, y, z;
                cell.position = offset + indices * cell_size + to_cell_center;
                cell.density = 0;
                // Try two ray directions (temporary bug fix, see 'Initialize two different ray directions' above) 
                for (int i = 0; i < 2; i++) { 
                    // Cast ray and check for hits
                    Ray ray;
                    ray.origin = cell.position;
                    ray.direction = ray_directions[i];
                    Vector3d hitPoint;
                    Vector3d hit_normal;
                    bool hit = trace_ray(ray, triangles, hitPoint, hit_normal);

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
                densities[x * dim_x * dim_y + y * dim_y + z] = cell.density;
            }
        }
    }
}


bool node_exists(std::map<int, int>* node_coords, int coords) {
    bool exists = !(node_coords->find(coords) == node_coords->end());
    //cout << "node with coords " << coords << " exists?  " << to_string(exists) << endl;
    //cout << "node coords: ";
    //help::print_map(node_coords);

    return exists;
}

int add_node_if_not_exists(
    int x, int y, Vector3d offset, int dim_x, std::map<int, int>& node_coords, float cell_size,
    vector<string>& nodes, int& _node_idx)
{
    dim_x += 1; // Number of nodes along an axis = number of cells + 1
    int node_idx = _node_idx;
    if (!node_exists(&node_coords, x * dim_x + y)) {
        string node = 
            to_string(_node_idx) + " " + to_string(cell_size * x + offset(0)) + " "
            + to_string(cell_size * y + offset(1)) + "\n";
        nodes.push_back(node);
        node_coords[x * dim_x + y] = _node_idx;
        _node_idx++;
        //cout << "node for (" << x << ", " << y << " ) does not exist.Created new node " << _node_idx-1 << endl;
    }
    else {
        node_idx = node_coords[x * dim_x + y];
        //cout << "node for (" << x << ", " << y << " ) already existed" << endl;
    }

    return node_idx;
}

/* Return the tag that corresponds to the 
*/
int get_tag(map<uint32_t, uint32_t>* bounds, uint32_t element_idx, int type, int& tag) {
    int _tag = mvis::help::get_value(bounds, element_idx);
    if (_tag == -1) {
        tag++; // If no tag was specified simply give each element a unique tag by incrementing it
        return tag;
    }
    else {
        return _tag; // Use the found tag
    }
}

/* Find a node that has not yet been visited
*/
int find_unvisited_node(
    std::map<int, int>* boundary_nodes, std::vector<int>* ordered_boundary_nodes
) {
    for (auto [node_coord, node_idx] : *boundary_nodes) {
        bool visited = std::find(
            ordered_boundary_nodes->begin(), ordered_boundary_nodes->end(), node_idx) != ordered_boundary_nodes->end();
        if (!visited) {
            return node_coord;
        }
    }
}

/* Generate a 2d Finite Element mesh from the given binary density distribution
*/
void generate_2d_FE_mesh(
    int dim_x, int dim_y, Vector3d offset, double cell_size, uint32_t* densities, vector<string>& nodes,
    vector<vector<int>>& elements, map<uint32_t, uint32_t>* bounds
) {
    int node_idx = 1;
    std::map<int, int> node_coords;
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            int filled = densities[x * dim_x + y];
            if (filled) {
                // Create 4 nodes that will enclose the surface element to be created (if they do not yet exist)
                int node4_idx = add_node_if_not_exists(x, y, offset, dim_x, node_coords, cell_size, nodes, node_idx);
                int node1_idx = add_node_if_not_exists(x + 1, y, offset, dim_x, node_coords, cell_size, nodes, node_idx);
                int node2_idx = add_node_if_not_exists(x + 1, y + 1, offset, dim_x, node_coords, cell_size, nodes, node_idx);
                int node3_idx = add_node_if_not_exists(x, y + 1, offset, dim_x, node_coords, cell_size, nodes, node_idx);

                // Compute surface index (equals coordinates + 1)
                int surface_idx = x * dim_x + y + 1;

                // For the 2d case, surfaces are never boundary elements. Therefore, tag is given a default value of 1.                int tag = 1;
                int tag = 1;
                
                // Other default indices
                int type = 3;
                int physical_entity = 0;
                int no_tags = 2;

                // Create the surface element
                vector<int> surface = {
                    surface_idx, type, no_tags, physical_entity, tag, node1_idx, node2_idx, node3_idx, node4_idx
                };
                elements.push_back(surface);
            }
        }
    }

    // Find boundary nodes
    std::map<int, int> boundary_node_coords;
    int no_boundary_nodes = 0;
    for (auto const& [node_coord, node_idx] : node_coords) {
        int x = node_coord / (dim_x + 1);
        int y = node_coord % (dim_x + 1);

        // Nodes at an extreme of at least one axis are always boundary nodes
        if (x == 0 || y == 0 || x == dim_x || y == dim_y) {
            boundary_node_coords[node_coord] = node_idx;
            no_boundary_nodes++;
            continue;
        }

        // Nodes with at least one empty neighboring cell are boundary nodes
        bool empty = false;
        for (int x_offset = 0; x_offset < 2; x_offset++) {
            for (int y_offset = 0; y_offset < 2; y_offset++) {
                int cell_is_filled = densities[(x - 1 + x_offset) * dim_x + (y - 1 + y_offset)];
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

    //cout << "node coord of node 6: " << mvis::help::get_key(&boundary_node_coords, 6) << endl;

    // ----  Order boundary nodes according to occurrence along perimeter of the mesh ---- //
    // Initialize starting node and ordered boundary nodes maps
    int node_coord = boundary_node_coords.begin()->first;
    node_idx = boundary_node_coords.begin()->second;
    std::vector<int> ordered_boundary_node_coords = { node_coord };
    int x = node_coord / (dim_x + 1); 
    int y = node_coord % (dim_x + 1);
    int start_x = x; 
    int start_y = y;

    // Initialize 'previous' coordinates to an arbitrary neighboring location on the grid 
    int previous_x = x; 
    int previous_y = y - 1;
    if (previous_y < 0) previous_y = y + 1;

    // Trace perimeter by stepping from one boundary node to the direct neighbor that
    // was not visited in the previous iteration
    int i = 1;
    bool neighbor_found = false;
    int neighbor_idx, neighbor_coord, neighbor_x, neighbor_y;
    int _neighbor_idx, _neighbor_coord, _neighbor_x, _neighbor_y;
    vector<pair<int, int>> offsets = { pair(0,1), pair(1,0), pair(-1, 0), pair(0, -1) };
    while (i < (no_boundary_nodes + 1)) {
        // If current node is the same as the starting node of the perimeter walk, re-start walk on a
        // node that has not yet been visited (such a node must be part of another component, for example a hole)
        if (i > 1 && x == start_x && y == start_y) {
            node_coord = find_unvisited_node(&boundary_node_coords, &ordered_boundary_node_coords);
            int node_idx = boundary_node_coords[node_coord];
            ordered_boundary_node_coords.push_back(node_coord);
            x = node_coord / (dim_x + 1);
            y = node_coord % (dim_x + 1);
            start_x = x;
            start_y = y;
            cout << "starting perimeter walk on new component. New Idx " << node_idx << ", New start x, y: " << x << ", " << y << endl;

            // Re-initialize 'previous' coordinates to an arbitrary neighboring location on the grid
            previous_x = x;
            previous_y = y - 1;
            if (previous_y < 0) previous_y = y + 1;
        }

        // Check which of the other three neighboring grid indices contains a boundary node
        vector<pair<int,int>> valid_neighbors;
        for (auto _offset : offsets) {
            // Compute neighbor coordinates
            _neighbor_x = (x + _offset.first);
            _neighbor_y = (y + _offset.second);

            //cout << "checking neighbor coord: " << _neighbor_x << ", " << _neighbor_y << " of cur node " << x << ", " << y << endl;

            // Check coordinate validity
            if (_neighbor_x < 0 || _neighbor_y < 0 || _neighbor_x > dim_x || _neighbor_y > dim_y)
                continue; // coordinates outside design domain are invalid
            if (_neighbor_x == previous_x && _neighbor_y == previous_y)
                continue; // skip the boundary node we found in the previous iteration

            // Check whether line connecting current node to neighbor is blocked by two adjacent filled cells
            bool blocked = false;
            if (x == _neighbor_x && x != 0 && x != dim_x) // Vertical line
            {
                int y_diff = _neighbor_y - y;
                int left_cell   = densities[(x - 1) * dim_x + y + y_diff];
                int right_cell  = densities[x * dim_x + y + y_diff];
                blocked = right_cell && left_cell;
                //cout << "Vertical line blocked. x,y: " << x << ", " << y << "  and  neighbor x,y: " << _neighbor_x << ", " << _neighbor_y << endl;
            }
            else if (y == _neighbor_y && y != 0 && y != dim_y) // Horizontal line
            {
                int x_diff = _neighbor_x - x;
                int up_cell     = densities[(x + x_diff) * dim_x + y];
                int down_cell   = densities[(x + x_diff) * dim_x + y - 1];
                blocked = up_cell && down_cell;
                //cout << "Horizontal line blocked. x,y: " << x << ", " << y << "  and  neighbor x,y: " << _neighbor_x << ", " << _neighbor_y << endl;
            }
            if (blocked) continue; // Skip neighbor if line is blocked by adjacent filled cells

            //cout << "coordinates valid. " << endl;

            // Check whether the neighbor coordinates contain a boundary node
            _neighbor_coord = (x + _offset.first) * (dim_x + 1) + (y + _offset.second);
            _neighbor_idx = mvis::help::get_value(&boundary_node_coords, _neighbor_coord);
            neighbor_found = _neighbor_idx != -1;
            if (neighbor_found) {
                neighbor_x = _neighbor_x;
                neighbor_y = _neighbor_y;
                neighbor_idx = _neighbor_idx;
                neighbor_coord = _neighbor_coord;
                valid_neighbors.push_back(pair(_neighbor_coord, _neighbor_idx));
            }
        }

        // Check whether more than one valid neighbor was found. If so, choose a neighbor that has not been visited yet
        if (valid_neighbors.size() > 1) {
            int j = 0;
            while (j < valid_neighbors.size()) {
                bool neighbor_visited = mvis::help::is_in(&ordered_boundary_node_coords, valid_neighbors[j].second);
                if (!neighbor_visited) {
                    break; // Found a valid neighbor that hasn't been visited yet
                }
                j++;
            }
            if (j == valid_neighbors.size()) {
                // No unvisited neighbor was found. 
                j = 0;
            }

            // Set the neighbor data to the data of the unvisited neighbor
            neighbor_coord = valid_neighbors[j].first;
            neighbor_idx = valid_neighbors[j].second;
            neighbor_x = node_coord / (dim_x + 1);
            neighbor_y = node_coord % (dim_x + 1);
            cout << "Choosing unvisited node " << neighbor_x << ", " << neighbor_y << endl;
        }

        if (neighbor_idx == -1) {
            cout << "ERROR: Neighbor not found. Previous x, y: " << previous_x << ", " << previous_y <<
                ", next x, y: " << neighbor_x << ", " << neighbor_y << endl;
        }

        // Add neighboring boundary node to vector
        ordered_boundary_node_coords.push_back(neighbor_coord);

        // Set previous x and y coordinates to current ones
        previous_x = x;
        previous_y = y;

        // Set next node coordinates to those of neighbor
        x = neighbor_x;
        y = neighbor_y;

        i++;
    }
    cout << "no of unordered bound nodes: " << boundary_node_coords.size() << endl;
    cout << "no of ordered bound nodes: " << ordered_boundary_node_coords.size() << endl;
    cout << "ordered bound nodes: "; mvis::help::print_vector(&ordered_boundary_node_coords);

    // ---- Generate boundary lines ---- //
    int tag = 1;
    for (int i = 1; i < ordered_boundary_node_coords.size(); i++) {
        // Get the 2 nodes that define the line
        int node1_coord = ordered_boundary_node_coords[i - 1];
        int node2_coord = ordered_boundary_node_coords[i];
        int node1_idx = boundary_node_coords[node1_coord];
        int node2_idx = boundary_node_coords[node2_coord];

        // Get line idx. Convention: line_idx = (<surface_coord> << 2) + <local_line_idx>
        // local_line_idx: Stepping clockwise from left edge of cell (local_line_idx = 0)
        int node1_x = node1_coord / (dim_x + 1);
        int node1_y = node1_coord % (dim_x + 1);
        int node2_x = node2_coord / (dim_x + 1);
        int node2_y = node2_coord % (dim_x + 1);
        int cell_x, cell_y;
        uint32_t cell_coord, local_line_idx;
        if (node1_x == node2_x) {   // Vertical line
            cell_x = node1_x;
            cell_y = min(node1_y, node2_y);
            local_line_idx = 0;
            if (node1_x == dim_x) { // Choose cell that has the line as its left side unless line is on the upper x-limit
                cell_x = node1_x - 1;
                local_line_idx = 2;
            }
        }
        else { // Horizontal line
            cell_y = node1_y;
            cell_x = min(node1_x, node2_x);
            local_line_idx = 3;
            if (node1_y == dim_y) { // Choose cell that has the line as its bottom side unless line is on the upper y-limit
                cell_y = node1_y - 1;
                local_line_idx = 1;
            }
        }
        // Use the last two bits of the integer to store the line idx. The first 30 are used for the cell coordinates
        cell_coord = cell_x * dim_x + cell_y;
        uint32_t line_idx = (cell_coord << 2) + local_line_idx;

        // Get the tag belonging to the surface element (if it has any)
        int _tag = get_tag(bounds, line_idx, 1, tag);

        // Unchanging indices
        int type = 1;
        int no_tags = 2;
        int physical_entity = 0;

        // Create the line element
        vector<int> line = { (int)line_idx + 1, type, no_tags, physical_entity, _tag, node1_idx, node2_idx };
        elements.push_back(line);
    }
}
