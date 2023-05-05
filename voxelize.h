#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>


using namespace Eigen;
using namespace std;

// Define the size of the grid
const int GRID_SIZE = 6;

// Define the struct for a cell in the grid
struct Cell {
    Vector3d position;
    int density;
};


void print_density_distrib(vector<vector<vector<int>>>* DensityDistrib) {
    int x_idx = 2;
    for (int y = 0; y < GRID_SIZE; y++) {
        for (int z = 0; z < GRID_SIZE; z++) {
            cout << DensityDistrib->at(x_idx)[y][z];
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

bool cast_ray(const Ray& ray, const std::vector<Triangle>& triangles, Vector3d& hitPoint, Vector3d& hit_normal) {
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

        if (dist > 1e-8 && dist < closestDist) {
            closestDist = dist;
            hit = true;
            hitPoint = ray.origin + ray.direction * dist;
            hit_normal = e1.cross(e2).normalized();
            /*cout << "hit normal: " << hit_normal.transpose() << endl;
            cout << "ray origin: " << ray.origin.transpose() << endl;
            cout << "ray direction: " << ray.direction.transpose() << endl;
            cout << "hit point: " << hitPoint.transpose() << endl;*/
            if (ray.origin[1] > 1.0) cout << "hit from point above cube (" << ray.origin.transpose() << ")" << endl;
        }
    }

    return hit;
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
    Vector3d margin; margin << 0.7, 0.7, 0.7;
    Vector3d minPoint = V->colwise().minCoeff();
    Vector3d maxPoint = V->colwise().maxCoeff();
    minPoint -= margin;
    maxPoint += margin;

    // Compute the barycenter of the mesh
    Vector3d mesh_barycent = V->colwise().mean();

    // Compute the size of each cell in the grid
    Vector3d gridSize = (maxPoint - minPoint) / static_cast<float>(GRID_SIZE);

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
    for (int x = 0; x < GRID_SIZE; x++) {
        vector<vector<int>> values_2d;
        for (int y = 0; y < GRID_SIZE; y++) {
            vector<int> values_1d;
            for (int z = 0; z < GRID_SIZE; z++) {
                Cell cell;
                Vector3d indices; indices << x, y, z;
                cell.position = minPoint + indices.cwiseProduct(gridSize);
                cell.density = 0;

                // Cast ray and check for hits
                Ray ray;
                ray.origin = cell.position;
                ray.direction << 0.0, 1.0, 0.0;
                Vector3d hitPoint;
                Vector3d hit_normal;
                bool hit = cast_ray(ray, triangles, hitPoint, hit_normal);

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
