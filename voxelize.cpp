#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>


using namespace Eigen;
using namespace std;

// Define the size of the grid
const int GRID_SIZE = 64;

// Define the struct for a cell in the grid
struct Cell {
    Vector3d position;
    int density;
};

/* Generate a binary density distribution on the grid based on the given mesh file
Input:
Returns (by reference):
    DensityDistrib (MatrixXi*): Matrix which contains a binary density value for each cell in the grid
*/
void generate_density_distribution(
    int dim_x, int dim_y, int dim_z, float cell_size, MatrixXd* V, MatrixXi* F,
    vector<vector<int>>& DensityDistrib, vector<vector<float>>& nodes, vector<vector<int>>& elements,
    vector<uint64_t>* bounds
) {
    // Compute the bounding box of the mesh
    Vector3d minPoint = V->colwise().minCoeff();
    Vector3d maxPoint = V->colwise().maxCoeff();

    // Compute the size of each cell in the grid
    Vector3d gridSize = (maxPoint - minPoint) / static_cast<float>(GRID_SIZE);

    // Assign density values to cells in the grid
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            for (int k = 0; k < GRID_SIZE; k++) {
                Cell cell;
                cell.position = minPoint + Vector3d(i, j, k) * gridSize;
                cell.density = 0;
                for (int index = 0; index < F->rows(); index ++) {
                    Vector3d v0 = V->row(F->coeff(index, 0));
                    Vector3d v1 = V->row(F->coeff(index, 1));
                    Vector3d v2 = V->row(F->coeff(index, 2));
                    Vector3d normal = (v1 - v0).cross(v2 - v0);
                    if (
                        (v0 - cell.position).dot(normal.transpose()) > 0.0 &&
                        (v1 - cell.position).dot(normal.transpose()) > 0.0 &&
                        (v2 - cell.position).dot(normal.transpose()) > 0.0)
                    {
                        cell.density = 1;
                        break;
                    }
                }
                cells.push_back(cell);
            }
        }
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
