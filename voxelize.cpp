#include <iostream>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>

// Define the size of the grid
const int GRID_SIZE = 64;

// Define the struct for a cell in the grid
struct Cell {
    glm::vec3 position;
    int density;
};

// Define the function to assign density values to cells in the grid
void assignDensityValues(std::vector<Cell>& cells, const std::vector<glm::vec3>& vertices, const std::vector<int>& indices) {
    // Compute the bounding box of the mesh
    glm::vec3 minPoint(std::numeric_limits<float>::max());
    glm::vec3 maxPoint(std::numeric_limits<float>::min());
    for (const auto& vertex : vertices) {
        minPoint = glm::min(minPoint, vertex);
        maxPoint = glm::max(maxPoint, vertex);
    }

    // Compute the size of each cell in the grid
    glm::vec3 gridSize = (maxPoint - minPoint) / static_cast<float>(GRID_SIZE);

    // Assign density values to cells in the grid
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            for (int k = 0; k < GRID_SIZE; k++) {
                Cell cell;
                cell.position = minPoint + glm::vec3(i, j, k) * gridSize;
                cell.density = 0;
                for (int index = 0; index < indices.size(); index += 3) {
                    glm::vec3 v0 = vertices[indices[index]];
                    glm::vec3 v1 = vertices[indices[index + 1]];
                    glm::vec3 v2 = vertices[indices[index + 2]];
                    glm::vec3 normal = glm::cross(v1 - v0, v2 - v0);
                    if (glm::dot(cell.position - v0, normal) > 0.0f && glm::dot(cell.position - v1, normal) > 0.0f && glm::dot(cell.position - v2, normal) > 0.0f) {
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
