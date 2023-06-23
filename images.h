#pragma once
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "meshing.h"

#define STB_IMAGE_IMPLEMENTATION
#include "lib/stb_image.h"


using namespace std;
using namespace fessga;


namespace fessga{

	class img {
	public:

		static void load_distribution_from_image(const char* filename, uint* densities, mesher::Grid3D grid) {
			// Load image
			int width, height, numChannels;
			unsigned char* image = stbi_load(filename, &width, &height, &numChannels, 0);
			if (!image) {
				std::cerr << "Failed to load image." << std::endl;
			}

			// Extract red values
			int numPixels = width * height;
			int pixels_per_cell = height / grid.y;
			int* redValues = new int[numPixels];
			for (int i = 0; i < numPixels; ++i) {
				int index = i * numChannels;
				redValues[i] = image[index];
			}
			stbi_image_free(image);

			// Sample red values with intervals
			for (int y = 0; y < grid.y; y++) {
				for (int x = 0; x < grid.x; x++) {
					int img_idx = y * width * pixels_per_cell + x * pixels_per_cell + (int)(0.5 * width * pixels_per_cell);
					densities[(x + 1) * grid.y - y - 1] = redValues[img_idx] / 255;
				}
			}
		}
	};
}
