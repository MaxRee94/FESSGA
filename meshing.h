#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <map>
#include "densities.h"

#define TEST2D true;
#define TEST3D false;
#define TESTCROSSOVER true;

using namespace Eigen;
using namespace std;


namespace fessga {
    class msh
    {
    public:

        // Define the struct for a Surface mesh
        class SurfaceMesh {
        public:
            SurfaceMesh() = default;
            SurfaceMesh(MatrixXd* _V, MatrixXi* _F) {
                V = _V;
                F = _F;
                bounding_box.row(0) = V->colwise().minCoeff();
                bounding_box.row(1) = V->colwise().maxCoeff();
                diagonal = bounding_box.row(1) - bounding_box.row(0);
                offset = -0.5 * diagonal;
                cout << "diagonal: " << diagonal.transpose() << endl;
                cout << "offset: " << offset.transpose() << endl;
            }
            SurfaceMesh(Vector3d size) {
                offset = -0.5 * size;
                diagonal = size;
                bounding_box.row(0) = 0.5 * size;
                bounding_box.row(1) = -0.5 * size;
                cout << "diagonal: " << diagonal.transpose() << endl;
                cout << "offset: " << offset.transpose() << endl;
            }
            SurfaceMesh(Vector2d size) {
                diagonal = Vector3d(size(0), size(1), 0);
                offset = -0.5 * diagonal;
                bounding_box.row(0) = 0.5 * diagonal;
                bounding_box.row(1) = -0.5 * diagonal;
                cout << "diagonal: " << diagonal.transpose() << endl;
                cout << "offset: " << offset.transpose() << endl;
            }
            Vector3d offset;
            Vector3d diagonal;
            MatrixXd bounding_box = MatrixXd(2, 3);
            MatrixXd* V = 0;
            MatrixXi* F = 0;
        };

        // Define the struct for an element
        struct Element {
            int id = 0;
            int type = 0;
            int no_tags = 0;
            int body = 0;
            int tag = 0;
            int boundary_id = -1; // Default -1 indicates the element is not part of the boundary
            vector<int> nodes = {};
        };

        // Define the struct for a 2D Finite Element mesh
        struct FEMesh2D {
            vector<vector<double>> nodes;
            vector<Element> lines;
            vector<Element> surfaces;
        };

        static bool node_exists(std::map<int, int>* node_coords, int coords) {
            bool exists = !(node_coords->find(coords) == node_coords->end());
            return exists;
        }

        static int add_node_if_not_exists(
            int x, int y, Vector3d offset, int dim_y, std::map<int, int>& node_coords, Vector2d cell_size,
            vector<vector<double>>& nodes, int& _node_idx)
        {
            dim_y += 1; // Number of nodes along an axis = number of cells + 1
            int node_idx = _node_idx;
            if (!node_exists(&node_coords, x * dim_y + y)) {
                uint node_coord = x * dim_y + y;
                node_coords[node_coord] = _node_idx;
                vector<double> node = { (double)(node_coord + 1), cell_size(0) * x + offset(0), cell_size(1) * y + offset(1), 0.0};
                nodes.push_back(node);
                _node_idx++;
                //cout << "node for (" << x << ", " << y << " ) does not exist.Created new node " << _node_idx-1 << endl;
            }
            else {
                node_idx = node_coords[x * dim_y + y];
                //cout << "node for (" << x << ", " << y << " ) already existed" << endl;
            }

            return x * dim_y + y + 1;
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

        /* Generate a grid-based description of a FE mesh
        */
        static void create_FE_mesh(
            SurfaceMesh mesh, grd::Densities2d densities, FEMesh2D& fe_mesh, bool verbose = false
        );

        static string get_msh_element_description(Element element) {
            string description = "";
            description += {
                to_string(element.id) + " " + to_string(element.type) + " " + to_string(element.no_tags) + " "
                + to_string(element.body) + " " + to_string(element.tag)
            };
            for (int i = 0; i < element.nodes.size(); i++) {
                description += " " + to_string(element.nodes[i]);
            }
            description += "\n";

            return description;
        }

        // Encode FE mesh data into .msh-description
        static void create_msh_description(FEMesh2D* fe_mesh, string& msh) {
            // -- Format section
            msh = {
                "$MeshFormat\n"
                "2.0 0 8\n"
                "$EndMeshFormat\n"
            };

            // -- Nodes section
            msh += "$Nodes\n";
            msh += to_string(fe_mesh->nodes.size()) + "\n";          // Number of nodes
            for (int i = 0; i < fe_mesh->nodes.size(); i++) {        // List of nodes, with each node encoded as <node_idx> <x> <y> <z>
                vector<double> node = fe_mesh->nodes[i];
                for (int j = 0; j < node.size(); j++) {
                    if (j == 0) msh += to_string((int)node[j]) + " ";
                    else msh += to_string(node[j]) + " ";
                }
                msh += "\n";
            }
            msh += "$EndNodes\n";

            // -- Elements section
            msh += "$Elements\n";
            int no_elements = fe_mesh->lines.size() + fe_mesh->surfaces.size();
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
            for (int i = 0; i < fe_mesh->surfaces.size(); i++) {
                Element surface = fe_mesh->surfaces[i];
                msh += get_msh_element_description(surface);
            }
            for (int i = 0; i < fe_mesh->lines.size(); i++) {
                Element line = fe_mesh->lines[i];
                msh += get_msh_element_description(line);
            }

            // End elements section
            msh += "$EndElements\n";
        }

        // Encode the given FE mesh in a Gmesh .msh format and export it to the given output folder
        static void export_as_msh_file(FEMesh2D* fe_mesh, string base_folder) {
            // Encode the FE mesh in a .msh format
            string msh_description;
            msh::create_msh_description(fe_mesh, msh_description);

            // Export the .msh description to a .msh file
            string msh_output_path = base_folder + "/mesh.msh";
            IO::write_text_to_file(msh_description, msh_output_path);
        }

        // Generate a file description for an Elmer header file, based on the given FE mesh
        static void export_elmer_header(FEMesh2D* fe_mesh, string base_folder) {
            string header_description = {
                to_string(fe_mesh->nodes.size()) + " " + to_string(fe_mesh->surfaces.size()) + " " + to_string(fe_mesh->lines.size()) + "\n"
                + "2\n"
                + "404 " + to_string(fe_mesh->surfaces.size()) + "\n" // Number of surfaces
                + "202 " + to_string(fe_mesh->lines.size()) + "\n" // Number of lines
            };
            IO::write_text_to_file(header_description, base_folder + "/mesh.header");
        }

        // Generate a file description for an Elmer nodes file, based on the given FE mesh
        static void export_elmer_nodes(FEMesh2D* fe_mesh, string base_folder) {
            string nodes_description = "";
            for (int i = 0; i < fe_mesh->nodes.size(); i++) { // List of nodes, with each node encoded as <node_idx> <x> <y> <z>
                vector<double> node = fe_mesh->nodes[i];
                int idx = (int)node[0];
                string _node = to_string(idx) + " -1";
                for (int j = 0; j < 3; j++) {
                    _node += " " + to_string(node[j + 1]);
                }
                nodes_description += _node + "\n";
            }
            IO::write_text_to_file(nodes_description, base_folder + "/mesh.nodes");
        }

        // Generate a file description for an Elmer elements file, based on the given FE mesh
        static void export_elmer_elements(FEMesh2D* fe_mesh, string base_folder) {
            string elements_description = "";
            for (int i = 0; i < fe_mesh->surfaces.size(); i++) {
                Element surface = fe_mesh->surfaces[i];
                
                // Encode surface idx, body, and type
                string _surface = to_string(surface.id) + " " + to_string(surface.body + 1) + " 404";
                
                // Encode member node ids
                for (int j = 0; j < surface.nodes.size(); j++) {
                    _surface += " " + to_string(surface.nodes[j]);
                }

                elements_description += _surface + "\n";
            }
            IO::write_text_to_file(elements_description, base_folder + "/mesh.elements");
        }

        // Generate a file description for an Elmer boundary elements file, based on the given FE mesh
        static void export_elmer_boundary(FEMesh2D* fe_mesh, string base_folder) {
            string elements_description = "";
            for (int i = 0; i < fe_mesh->lines.size(); i++) {
                Element line = fe_mesh->lines[i];

                // Encode line idx, boundary id, p1, p2, and type
                int parent_id = line.id >> 2 + 1;
                string _line = to_string(line.id) + " " + to_string(line.boundary_id) + " " + to_string(parent_id) + " 0 202";

                // Encode member node ids
                for (int j = 0; j < line.nodes.size(); j++) {
                    _line += " " + to_string(line.nodes[j]);
                }

                elements_description += _line + "\n";
            }
            IO::write_text_to_file(elements_description, base_folder + "/mesh.boundary");
        }

        // Export FE mesh as elmer files (.header, .boundaries, .nodes, .elements)
        static void export_as_elmer_files(FEMesh2D* fe_mesh, string base_folder) {
            export_elmer_header(fe_mesh, base_folder);
            export_elmer_nodes(fe_mesh, base_folder);
            export_elmer_elements(fe_mesh, base_folder);
            export_elmer_boundary(fe_mesh, base_folder);
        }

        // Create batch file for running elmer and return its absolute path
        static string create_batch_file(string base_folder) {
            string batch_file = IO::get_fullpath(base_folder);
            batch_file += "/run_elmer.bat";
            IO::write_text_to_file("cd \"" + base_folder + "\"\nElmerSolver", batch_file);

            return batch_file;
        }

        // Parse the target boundaries, yielding a vector of boundary ids
        static void parse_target_boundaries(string line, string& prefix, vector<int>& bound_ids) {
            // Split the prefix and the numbers-part of the line
            vector<string> split_line;
            help::split(line, " = ", split_line);


            prefix = split_line[0] + " = ";

            // Split the numbers-part of the line
            //cout << "sec part of line: " << split_line[1] << endl;
            string numbers = split_line[1];
            vector<string> number_strings;
            help::split(numbers, " ", number_strings);

            // Convert each substring to an integer and add it to the boundary ids vector
            for (int i = 0; i < number_strings.size(); i++) {
                //cout << "num string: " << number_strings[i] << endl;
                bound_ids.push_back(stoi(number_strings[i]));
            }
        }

        // Parse the given case.sif file, extracting boundary ids and in-between sections of text
        static void read_boundary_conditions(
            map<string, vector<int>>& boundary_id_lookup, phys::FEACase* fea_case, phys::FEACaseManager* fea_casemanager
        ) {
            // Read content
            vector<string> case_content;
            IO::read_file_content(fea_case->path, case_content);

            // Get target boundaries and in-between sections
            bool target_boundaries = false;
            bool bound_name = false;
            string section = "";
            vector<int> bound_ids;
            fea_case->sections.clear();
            fea_case->names.clear();
            for (int i = 0; i < case_content.size(); i++) {
                string line = case_content[i];
                if (help::is_in(line, "Boundary Condition")) {
                    target_boundaries = true;
                }
                else if (target_boundaries) {
                    target_boundaries = false;
                    bound_name = true;
                    string prefix;
                    parse_target_boundaries(line, prefix, bound_ids);
                    section += prefix;
                    fea_case->sections.push_back(section);
                    section = "\n";
                    continue;
                }
                else if (bound_name) {
                    bound_name = false;
                    vector<string> split_line;
                    help::split(line, " = ", split_line);
                    string name = "";
                    int last_idx = -1;
                    for (string::iterator it = split_line[1].begin(); it != split_line[1].end(); it++) {
                        last_idx++;
                    }
                    for (int j = 1; j < last_idx; j++) {
                        name += split_line[1][j];
                    }
                    fea_case->names.push_back(name);
                    boundary_id_lookup[name] = bound_ids;
                    bound_ids.clear();
                }
                section += line + "\n";
            }
            fea_case->sections.push_back(section);
        }

        // Create a map from strings that encode the name of a boundary condition to vectors of boundary node coordinate pairs that
        // represent the edges to which those boundary conditions have been applied
        static void derive_boundary_conditions(
            phys::FEACase& fea_case, grd::Densities2d& densities, SurfaceMesh mesh, phys::FEACaseManager* fea_casemanager
        ) {
            // Read each boundary condition and the corresponding boundary node id's from the original case.sif file
            map<string, vector<int>> boundary_id_lookup;
            read_boundary_conditions(boundary_id_lookup, &fea_case, fea_casemanager);

            // Generate FE mesh
            FEMesh2D fe_mesh;
            create_FE_mesh(mesh, densities, fe_mesh);
            
            // For each boundary condition, add the vector of boundary node coordinates in the fe mesh to the bound_cond_lines map
            int q = 0;
            for (auto& [bound_name, bound_ids] : boundary_id_lookup) {
                vector<pair<int, int>> bound_lines;
                map<string, vector<int>> bound_cell_map;
                vector<int> bound_cells;
                vector<int> cutout_cells;
                vector<int> keep_cells;
                for (int i = 0; i < bound_ids.size(); i++) {
                    int bound_id = bound_ids[i];
                    Element line = fe_mesh.lines[bound_id - 1];
                    bound_lines.push_back(pair(line.nodes[0] - 1, line.nodes[1] - 1));
                    q++;

                    // Store the parent cell coordinates; this cell should not be removed during optimization
                    int cell_coord = (line.id - 1) >> 2;
                    if (!help::is_in(&fea_casemanager->cells_to_keep, cell_coord)) {
                        fea_casemanager->cells_to_keep.push_back(cell_coord);
                        bound_cells.push_back(cell_coord);
                    }

                    // Store the void cells adjacent to the boundary cell; these should remain empty so that the boundary line is maintained
                    int local_line_idx = line.id % 4;
                    vector<int> void_neighbors = densities.get_empty_neighbors(cell_coord, true);
                    for (auto& void_neighbor : void_neighbors) {
                        if (!help::is_in(&fea_casemanager->cutout_cells, void_neighbor)) {
                            fea_casemanager->cutout_cells.push_back(void_neighbor);
                            cutout_cells.push_back(void_neighbor);
                        }
                    }

                    // Store the filled cells neighboring the boundary cell as cells to keep
                    vector<int> neighbors = densities.get_neighbors(cell_coord);
                    for (auto& neighbor : neighbors) {
                        if (!help::is_in(&fea_casemanager->cells_to_keep, neighbor)) {
                            fea_casemanager->cells_to_keep.push_back(neighbor);
                            keep_cells.push_back(neighbor);
                        }
                    }
                }
                fea_case.bound_cond_lines[bound_name] = bound_lines;
                bound_cell_map["bound"] = bound_cells;
                bound_cell_map["cutout"] = cutout_cells;
                bound_cell_map["keep"] = keep_cells;
                fea_case.bound_cond_cells[bound_name] = bound_cell_map;
            }
            cout << "no boundary cells: " << fea_casemanager->cells_to_keep.size() << endl;
            cout << "no boundary lines: " << q << endl;
        }

        // Create map containing a vector of boundary ids corresponding to the given fe mesh for each boundary condition name
        static void create_bound_id_lookup(
            map<string, vector<pair<int, int>>>* bound_cond_lines, FEMesh2D* fe_mesh, map<string, vector<int>>& bound_id_lookup
        ) {
            for (auto& [bound_name, lines] : (*bound_cond_lines)) {
                // Get the boundary ids of all node coordinate pairs stored in the bound_cond_lines-map for the current boundary condition.
                vector<int> bound_ids;
                for (int i = 0; i < lines.size(); i++) {
                    pair<int, int> line = lines[i];
                    bool found = 0;
                    for (int j = 0; j < fe_mesh->lines.size(); j++) {

                        // Compare fe line node coordinates to the line coordinates stored in the boundary conditions map.
                        // Check both possible permutations of coordinates ({node1, node2} and {node2, node1}), because order does not matter.
                        int fe_line_node1 = fe_mesh->lines[j].nodes[0] - 1;
                        int fe_line_node2 = fe_mesh->lines[j].nodes[1] - 1;
                        if (
                            (line.first == fe_line_node1 && line.second == fe_line_node2) ||
                            (line.first == fe_line_node2 && line.second == fe_line_node1))
                        {
                            bound_ids.push_back(fe_mesh->lines[j].boundary_id);
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        cerr << "ERROR: Failed to find boundary id for line (" << line.first / 200 << ", " << line.first % 200 << "), ("
                            << line.second / 200 << ", " << line.second % 200 << ")" << endl;
                        throw std::runtime_error("Failed to find boundary id");
                    }
                }
                bound_id_lookup[bound_name] = bound_ids;
            }
        }

        // Re-assemble the case file's content by concatenating the sections and updated target boundaries
        static void assemble_fea_case(
            phys::FEACaseManager* fea_casemanager, phys::FEACase* fea_case, map<string, vector<int>>* bound_id_lookup
        ) {
            fea_case->content = "";
            for (int i = 0; i < fea_case->names.size(); i++) {
                string name = fea_case->names[i];
                vector<int> bound_ids = bound_id_lookup->at(name);
                string bound_ids_string = help::join_as_string(bound_ids, " ");
                fea_case->content += fea_case->sections[i] + bound_ids_string;
            }
            fea_case->content += fea_case->sections[fea_case->names.size()];

            // Replace 'Output File Name' default setting of 'case' with the fea_case's specific case name.
            fea_case->content = help::replace_occurrences(
                fea_case->content, "Output File Name = case\n", "Output File Name = " + fea_case->name + "_\n");
        }

        static void get_vtk_paths(phys::FEACaseManager* fea_casemanager, string output_folder, vector<string>& vtk_paths) {
            for (auto& fea_case : fea_casemanager->active_cases) {
                vtk_paths.push_back(output_folder + "/" + fea_case.name + "_0001.vtk");
            }
        }

        static void init_fea_cases(
            phys::FEACaseManager* fea_casemanager, string case_folder, vector<string> case_names,
            grd::Densities2d* densities, SurfaceMesh* mesh, bool source_cases = true
        ) {
            vector<phys::FEACase>* fea_cases;
            if (source_cases) fea_cases = &fea_casemanager->sources;
            else fea_cases = &fea_casemanager->targets;
            for (auto& case_name : case_names) {
                string case_path = case_folder + "/" + case_name + ".sif";
                phys::FEACase fea_case(
                    case_path, densities->dim_x, densities->dim_y, fea_casemanager->max_stress_threshold
                );
                fea_case.name = case_name;
                msh::derive_boundary_conditions(fea_case, *densities, *mesh, fea_casemanager);
                fea_cases->push_back(fea_case);
            }
        }

    };
}
