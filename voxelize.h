#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>


using namespace Eigen;
using namespace std;

// Define the size of the grid
const int GRID_SIZE = 5;

// Define the struct for a cell in the grid
struct Cell {
    Vector3d position;
    int density;
};


void print_density_distrib(vector<vector<vector<int>>>* DensityDistrib) {
    int z_idx = 2;
    for (int x = 0; x < DensityDistrib->at(z_idx).size(); x++) {
        for (int y = 0; y < DensityDistrib->at(z_idx)[0].size(); y++) {
            cout << DensityDistrib->at(0)[x][y];
        }
        cout << endl;
    }
}


/* Generate a binary density distribution on the grid based on the given mesh file
Input:
Returns (by reference):
    DensityDistrib (MatrixXi*): Matrix which contains a binary density value for each cell in the grid
*/
void generate_density_distribution(
    int dim_x, int dim_y, int dim_z, float cell_size, MatrixXd* V, MatrixXi* F,
    vector<vector<vector<int>>>& DensityDistrib, vector<vector<float>>& nodes, vector<vector<int>>& elements,
    vector<uint64_t>* bounds
) {
    // Compute the bounding box of the mesh
    Vector3d margin; margin << 0.5, 0.5, 0.5;
    Vector3d minPoint = V->colwise().minCoeff();
    Vector3d maxPoint = V->colwise().maxCoeff();
    minPoint -= margin;
    maxPoint += margin;

    // Compute the size of each cell in the grid
    Vector3d gridSize = (maxPoint - minPoint) / static_cast<float>(GRID_SIZE);

    cout << "gridSize: " << gridSize << endl;
    cout << "minPoint: " << minPoint << endl;
    cout << "maxPoint: " << maxPoint << endl;

    // Assign density values to cells in the grid
    for (int x = 0; x < GRID_SIZE; x++) {
        vector<vector<int>> values_2d;
        for (int y = 0; y < GRID_SIZE; y++) {
            vector<int> values_1d;
            for (int z = 0; z < GRID_SIZE; z++) {
                Cell cell;
                Vector3d indices; indices << x, y, z;
                cell.position = minPoint + indices.cwiseProduct(gridSize);
                cell.density = 0;
                float min_dist = 10e8;
                int closest_face = 0;

                // Find closest face
                for (int face_idx = 0; face_idx < F->rows(); face_idx ++) {
                    Vector3d v0 = V->row(F->coeff(face_idx, 0));
                    Vector3d v1 = V->row(F->coeff(face_idx, 1));
                    Vector3d v2 = V->row(F->coeff(face_idx, 2));

                    // Check against smallest found distance.
                    // If smaller, update smallest distance and face_idx of closest face.
                    float dist_cell2face = (((v0 + v1 + v2) / 3) - cell.position).norm();
                    if (dist_cell2face < min_dist) {
                        min_dist = dist_cell2face;
                        closest_face = face_idx;
                    }
                }

                // Check if normal of closest face points in the same direction as the vectors from the cell to each of the vertices
                Vector3d v0 = V->row(F->coeff(closest_face, 0));
                Vector3d v1 = V->row(F->coeff(closest_face, 1));
                Vector3d v2 = V->row(F->coeff(closest_face, 2));
                Vector3d barycent = (v0 + v1 + v2) / 3;
                Vector3d normal = (v1 - v0).cross(v2 - v0);
                if ((barycent - cell.position).dot(normal) > 0.0) {
                    cell.density = 1;
                    cout << endl;
                    cout << "v0: " << v0.transpose() << endl;
                    cout << "v1: " << v1.transpose() << endl;
                    cout << "v2: " << v2.transpose() << endl;
                    cout << "vec to barycent: " << (barycent - cell.position).transpose() << endl;
                    cout << "normal: " << normal.transpose() << endl;
                    cout << "dens=1 for cell position " << cell.position.transpose() << " and face " << ((v0 + v1 + v2) / 3).transpose() << endl;
                }
                else {
                    //cout << "dens=0 for cell position " << cell.position.transpose() << " and face " << ((v0+v1+v2)/3).transpose() << endl;
                }
                values_1d.push_back(cell.density);
            }
            values_2d.push_back(values_1d);
        }
        DensityDistrib.push_back(values_2d);
    }
}


//int main() {
//    // Read in the mesh file
//    std::vector<glm::vec3> vertices;
//    std::vector<int> indices;
//
//    // Assign density values to cells in the grid
//    std::vector<Cell> cells;
//    assignDensityValues(cells, vertices, indices);
//
//    // Print the density values of the cells
//    for (const auto& cell : cells) {
//        std::cout << cell.density << " ";
//    }
//    std::cout << std::endl;
//
//    return 0;
//}
