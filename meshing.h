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

    // Compute the barycenter of the mesh
    Vector3d mesh_barycent = V->colwise().mean();

    // Compute the size of each cell in the grid
    Vector3d grid_size = Vector3d(dim_x, dim_y, dim_z);

    // Compute vector to center of a grid cell from its corner
    Vector3d to_cell_center = Vector3d(0.5, 0.5, 0.5) * cell_size;

    // Create list of triangles
    std::vector<Triangle> triangles;
    for (int face_idx = 0; face_idx < F->rows(); face_idx++) {
        Triangle triangle;
        triangle.v0 = V->row(F->coeff(face_idx, 0));
        triangle.v1 = V->row(F->coeff(face_idx, 1));
        triangle.v2 = V->row(F->coeff(face_idx, 2));
        triangles.push_back(triangle);
    }

    // Set ray direction
    Vector3d ray_dir = Vector3d(1.0, 1.0, 1.0).normalized();

    // Assign density values to cells in the grid
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            for (int z = 0; z < dim_z; z++) {
                Cell cell;
                Vector3d indices; indices << x, y, z;
                cell.position = offset + indices * cell_size + to_cell_center;
                cell.density = 0;

                // Cast ray and check for hits
                Ray ray;
                ray.origin = cell.position;
                ray.direction = ray_dir;
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

int get_2d_tag(vector<uint64_t>* bounds, int element_idx, int type) {
    int tag = 1;
    if (type == 3) {
        if (bounds->size() > 0) {
            // TODO: Check if element index corresponds to an item in bounds, and if so, which tag should be applied
        }
    }
    return tag;
}

/* Generate a 2d Finite Element mesh from the given binary density distribution
*/
void generate_2d_FE_mesh(
    int dim_x, int dim_y, Vector3d offset, double cell_size, uint32_t* densities, vector<string>& nodes,
    vector<vector<int>>& elements, vector<uint64_t>* bounds
) {
    int node_idx = 1;
    int surface_idx = 1;
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

                // Get the tag belonging to the surface element (if it has any)
                // This information is stored in the bounds vector
                int tag = get_2d_tag(bounds, surface_idx, 3);

                // Create the surface element
                vector<int> surface = { surface_idx, 3, 2, 1, tag, node1_idx, node2_idx, node3_idx, node4_idx };
                elements.push_back(surface);
                surface_idx++;
            }
        }
    }

    // Find boundary nodes
    std::map<int, int> boundary_node_coords;
    int no_boundary_nodes = 0;
    for (auto const& [node_idx, node_coord] : node_coords) {
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

    // Order boundary nodes according to occurrence along perimeter of the mesh
    int node_coord = boundary_node_coords.begin()->first;
    node_idx = boundary_node_coords.begin()->second;
    std::vector<int> ordered_boundary_node_coords = {node_idx};

    // Initialize previous coordinates to an arbitrary neighboring location on the grid
    int x = node_coord / (dim_x + 1);
    int y = node_coord % (dim_x + 1);
    int previous_x = x;
    int previous_y = y - 1;
    if (previous_y < 0) previous_y = y + 1;

    // Trace perimeter by stepping from one boundary node to the direct neighbor that
    // was not visited in the previous iteration
    int i = 0;
    while (i < no_boundary_nodes) {
        // Check which of the other three neighboring grid indices contains a boundary node
        bool neighbor_found = false;
        int neighbor_idx, neighbor_coord, neighbor_x, neighbor_y;
        for (int x_offset = 0; x_offset < 3; x_offset++) {
            for (int y_offset = 0; y_offset < 3; y_offset++) {
                // Skip coordinates on opposite corner of neighboring cell
                if (x_offset == y_offset == 1 || x_offset == y_offset == -1) continue;

                // Compute neighbor coordinates
                neighbor_x = (x + x_offset);
                neighbor_y = (y + y_offset);

                // Check coordinate validity
                if (neighbor_x < 0 || neighbor_y < 0 || neighbor_x > dim_x || neighbor_y > dim_y)
                    continue; // coordinates outside design domain are invalid
                if (neighbor_x == previous_x && neighbor_y == previous_y)
                    continue; // skip the boundary node we found in the previous iteration

                // Check whether the neighbor coordinates contain a boundary node
                neighbor_coord = (x + x_offset) * dim_x + (y + y_offset);
                neighbor_idx = mvis::help::get_value(&boundary_node_coords, neighbor_coord);
                neighbor_found = neighbor_idx != -1;
                if (neighbor_found) {
                    break;
                }
            }
            if (neighbor_found) break;
        }

        // Get node coordinates of neighbor
        x = node_coord / (dim_x + 1);
        y = node_coord % (dim_x + 1);

        // Add neighboring boundary node to vector
        ordered_boundary_node_coords.push_back(neighbor_idx);

        // Set previous x and y coordinates to current ones
        previous_x = x;
        previous_y = y;

        // Set next node coordinates to those of neighbor
        x = neighbor_x;
        y = neighbor_y;

        i++;
    }
    cout << "no of ordered bound elements: " << i << endl;
    cout << "ordered bound elements: ";
    mvis::help::print_vector(&ordered_boundary_node_coords);
}
