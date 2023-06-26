#pragma once
#include <Eigen/Core>
#include <vector>
#include <map>
#include <cstdlib>
#include <set>

using namespace Eigen;


namespace fessga {
	class tracer {
	public:
        struct Ray {
	        Vector3d origin;
	        Vector3d direction;
        };

        struct Triangle {
	        Vector3d v0 = Vector3d(0.0, 0.0, 0.0);
	        Vector3d v1 = Vector3d(0.0, 0.0, 0.0);
	        Vector3d v2 = Vector3d(0.0, 0.0, 0.0);
        };

        struct Cell3D {
	        Vector3d position;
	        int density;
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
	};
}
