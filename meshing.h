#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <map>
#include "helpers.h"

#define TEST2D true;
#define TEST3D false;
#define TESTCROSSOVER true;

using namespace Eigen;
using namespace std;


namespace fessga {
    class mesher
    {
    public:
        // Define the struct for a 3d Grid
        struct Grid3D {
            int dim_x, dim_y, dim_z;
        };

        // Define the struct for a 2d Grid
        struct Grid2D {
            int dim_x, dim_y;
        };

        // Define the struct for a Surface mesh
        struct SurfaceMesh {
            Vector3d offset;
            MatrixXd bounding_box = MatrixXd(2, 3);
            MatrixXd* V = 0;
            MatrixXi* F = 0;
        };

        // Define the struct for a 2D Finite Element mesh
        struct FEMesh2D {
            map<uint32_t, uint32_t> line_boundaries;
            vector<vector<int>> elements;
            vector<string> nodes;
        };

        // Define the struct for a cell in the grid
        struct Cell {
            Vector3d position;
            int density;
        };

        static void print_density_distrib(uint* densities, int grid_size, int z = -1) {
            if (z == -1) z = grid_size / 2;
            for (int x = 0; x < grid_size; x++) {
                for (int y = 0; y < grid_size; y++) {
                    cout << densities[z * grid_size * grid_size + x * grid_size + y];
                }
                cout << endl;
            }
        }

        static void print_2d_density_distrib(uint* densities, int dim_x, int dim_y) {
            for (int y = dim_y -1; y > -1; y--) {
                for (int x = 0; x < dim_x; x++) {
                    cout << densities[x * dim_x + y];
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

        static bool trace_ray(const Ray& ray, const std::vector<Triangle>& triangles, Vector3d& hitPoint, Vector3d& hit_normal) {
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
            densities (uint*): Array which contains a binary density value (0 or 1) for each cell in the grid
        */
        static void generate_3d_density_distribution(
            int dim_x, int dim_y, int dim_z, Vector3d offset, double cell_size, MatrixXd* V, MatrixXi* F,
            uint* densities
        ) {
            cout << "Generating 3d grid-based density distribution..." << endl;

            // Compute the barycenter of the mesh
            Vector3d mesh_barycent = V->colwise().mean();

            // Compute the size of each cell in the grid
            Vector3d grid_size = Vector3d(dim_x, dim_y, dim_z);

            // Compute vector to center of a grid cell from its corner
            Vector3d to_cell_center = Vector3d(0.5, 0.5, 0.5) * cell_size;

            // Initialize different ray directions, (this is a temporary fix for a bug whereby some cells
            // are not properly assigned a density of 1)
            vector<Vector3d> ray_directions = { Vector3d(0, 1.0, 0), Vector3d(1.0, 1.0, 1.0).normalized(), Vector3d(1.0, 0, 0) };

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
            int slices_done = 0;
#pragma omp parallel for
            for (int x = 0; x < dim_x; x++) {
                cout << "slices_done: " << slices_done << " / " << dim_x << endl;
                slices_done++;
                for (int y = 0; y < dim_y; y++) {
                    for (int z = 0; z < dim_z; z++) {
                        Cell cell;
                        Vector3d indices; indices << x, y, z;
                        cell.position = offset + indices * cell_size + to_cell_center;
                        cell.density = 0;
                        // Try three ray directions (temporary bug fix, see 'Initialize two different ray directions' above) 
                        for (int i = 0; i < 3; i++) {
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
                                //cout << "density 1" << endl;
                                cell.density = 1;
                                break;
                            }
                        }
                        densities[x * dim_x * dim_y + y * dim_y + z] = cell.density;
                    }
                }
            }
            cout << "Finished generating density distribution. Filtering out floating cells..." << endl;

            // Filter out floating cells that have no direct neighbors
#pragma omp parallel for
            for (int x = 0; x < dim_x; x++) {
                for (int y = 0; y < dim_y; y++) {
                    for (int z = 0; z < dim_z; z++) {
                        int filled = densities[x * dim_x * dim_y + y * dim_y + z];
                        if (!filled) continue;
                        int neighbor = 0;
                        for (int _x = -1; _x <= 1; _x++) {
                            for (int _y = -1; _y <= 1; _y++) {
                                for (int _z = -1; _z <= 1; _z++) {
                                    if (x+_x == dim_x || y + _y == dim_y || z+_z == dim_z) continue;
                                    if (x+_x == 0 || y+_y == 0 || z+_z == 0) continue;
                                    neighbor = densities[(x+_x) * dim_x * dim_y + (y+_y) * dim_y + (z+_z)];
                                    if (neighbor) break;
                                }
                                if (neighbor) break;
                            }
                            if (neighbor) break;
                        }
                        if (!neighbor) {
                            cout << "floating cell detected. Setting to 0" << endl;
                            densities[x * dim_x * dim_y + y * dim_y + z] = 0; // Set density to 0 if the cell has no neighbors
                        }
                    }
                }
            }
            cout << "Finished filtering floating cells." << endl;
        }

        static void filter_2d_density_distrib(uint* densities, int dim_x, int dim_y) {
            // Filter out floating cells that have no direct neighbors
            cout << "Filtering 2d floating cells..." << endl;
#pragma omp parallel for
            for (int x = 0; x < dim_x; x++) {
                for (int y = 0; y < dim_y; y++) {
                    int filled = densities[x * dim_y + y];
                    if (!filled) continue;
                    int neighbor = 0;
                    for (int _x = -1; _x <= 1; _x++) {
                        for (int _y = -1; _y <= 1; _y++) {
                            if (_y == 0 && _x == 0) continue;
                            if (x + _x == dim_x || y + _y == dim_y) continue;
                            if (x + _x == 0 || y + _y == 0) continue;
                            neighbor = densities[(x + _x) * dim_y + (y + _y)];
                            if (neighbor) break;
                        }
                        if (neighbor) break;
                    }
                    if (!neighbor) {
                        cout << "floating cell detected. Setting to 0" << endl;
                        densities[x * dim_y + y] = 0; // Set density to 0 if the cell has no neighbors
                    }
                }
            }
            cout << "Finished filtering floating cells." << endl;
        }

        static bool node_exists(std::map<int, int>* node_coords, int coords) {
            bool exists = !(node_coords->find(coords) == node_coords->end());
            //cout << "node with coords " << coords << " exists?  " << to_string(exists) << endl;
            //cout << "node coords: ";
            //help::print_map(node_coords);

            return exists;
        }

        static int add_node_if_not_exists(
            int x, int y, Vector3d offset, int dim_x, std::map<int, int>& node_coords, float cell_size,
            vector<string>& nodes, int& _node_idx)
        {
            dim_x += 1; // Number of nodes along an axis = number of cells + 1
            int node_idx = _node_idx;
            if (!node_exists(&node_coords, x * dim_x + y)) {
                uint node_coord = x * dim_x + y;
                node_coords[node_coord] = _node_idx;
                string node =
                    to_string(node_coord + 1) + " " + to_string(cell_size * x + offset(0)) + " "
                    + to_string(cell_size * y + offset(1)) + "\n";
                nodes.push_back(node);
                _node_idx++;
                //cout << "node for (" << x << ", " << y << " ) does not exist.Created new node " << _node_idx-1 << endl;
            }
            else {
                node_idx = node_coords[x * dim_x + y];
                //cout << "node for (" << x << ", " << y << " ) already existed" << endl;
            }

            return x * dim_x + y + 1;
        }

        /* Return the tag that corresponds to the given element index
        */
        static int get_tag(map<uint, uint>* bounds, uint element_idx, int type, int& tag) {
            int _tag = fessga::help::get_value(bounds, element_idx);
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
        static int find_unvisited_node(
            std::map<int, int>* boundary_nodes, std::vector<int>* ordered_boundary_nodes
        ) {
            for (auto [node_coord, node_idx] : *boundary_nodes) {
                bool visited = std::find(
                    ordered_boundary_nodes->begin(), ordered_boundary_nodes->end(), node_coord) != ordered_boundary_nodes->end();
                if (!visited) {
                    return node_coord;
                }
            }
            return -1; // Return -1 if no unvisited node was found
        }

        /* Generate a 2d Finite Element mesh from the given binary density distribution
        */
        static void generate_2d_FE_mesh(
            int dim_x, int dim_y, Vector3d offset, double cell_size, uint* densities, vector<string>& nodes,
            vector<vector<int>>& elements, map<uint, uint>* bounds
        );

        static void generate_msh_description(FEMesh2D fe_mesh, string& msh) {
            // Encode mesh data into .msh-description
            // -- Format section
            msh = {
                "$MeshFormat\n"
                "2.0 0 8\n"
                "$EndMeshFormat\n"
            };

            // -- Nodes section
            msh += "$Nodes\n";
            msh += to_string(fe_mesh.nodes.size()) + "\n";          // Number of nodes
            for (int i = 0; i < fe_mesh.nodes.size(); i++) {        // List of nodes, with each node encoded as <node_idx> <x> <y> <z>
                msh += fe_mesh.nodes.at(i);
            }
            msh += "$EndNodes\n";

            // -- Elements section
            msh += "$Elements\n";
            int no_elements = fe_mesh.elements.size();
            msh += to_string(no_elements) + "\n";

            // Iterate over all elements and add to description string
            // Each element is encoded as <elm-number> <elm-type> <number-of-tags> <physical entity> <tag> <node_1> ... <node_n>
            // Element types:
            //      1 : 2-node line
            //      2 : 3-node triangle
            //      3 : 4-node quad
            //      4 : 4-node tetrahedron
            //      5 : 8-node hexahedron (a cube is a regular hexahedron)
            //      15: 1-node point
            // For example: 3 1 2 0 1 1 4 encodes a line from node 1 to node 4 with boundary tag 1
            for (int i = 0; i < fe_mesh.elements.size(); i++) {
                vector<int> element = fe_mesh.elements[i];
                for (int j = 0; j < element.size(); j++) {
                    msh += to_string(element[j]) + " ";
                }
                msh += "\n";
            }

            // End elements section
            msh += "$EndElements\n";
        }

        /* Generate a grid-based description of a FE mesh that can be output as a .msh file
        Input:
            dim_x, dim_y, dim_z (int):  Number of cells along each dimension of the grid
            csize (float): Size of a single cell (size=width=height=depth)
            V (MatrixXd*): Pointer to the matrix of vertex positions for the given mesh
            F (MatrixXi*): Pointer to the matrix of faces for the given mesh
            msh (string): The generated description string in .msh-format
        */
        static void convert_to_FE_mesh(
            const int dim_x, const int dim_y, const int dim_z, const float cell_size, Vector3d offset, uint32_t* densities, FEMesh2D& fe_mesh)
        {
            map<uint32_t, uint32_t> line_bounds;
            vector<vector<int>> elements;
            vector<string> nodes;
            generate_2d_FE_mesh(dim_x, dim_y, offset, cell_size, densities, nodes, elements, &line_bounds);
            
            fe_mesh.elements = elements;
            fe_mesh.line_boundaries = line_bounds;
            fe_mesh.nodes = nodes;
        }
    };
}
