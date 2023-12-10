from argparse import ArgumentParser
import numpy as np
import os


class Distribution:
    def __init__(self, filepath):
        self.path = filepath
        self.name = filepath
        self.values = None
        self.stats = None
        self.dim_x = None
        self.dim_y = None
        self.domain_size = None

    def create_array(self, dim_x, dim_y):
        self.dim_x = dim_x
        self.dim_y = dim_y
        self.domain_size = dim_x * dim_y
        self.values = np.zeros((dim_x, dim_y))

    def __getitem__(self, coord):
        return self.values[coord]

    def __setitem__(self, coord, value):
        self.values[coord] = value

    def count(self):
        return int(self.values.sum())

    def compute_stats(self):
        self.stats = {
            "count": self.count(),
            "relative area": self.count() / self.domain_size,
        }

    def print_stats(self):
        print("Statistics for distribution", self.name, ":")
        for stat, value in self.stats.items():
            print("  " + stat, ":", value)

    def compare(self, distr2):
        overlap = 0
        for x in range(self.dim_x):
            for y in range(self.dim_y):
                overlap += int(self.values[x, y] == distr2.values[x, y])

        print("-" * 30 + "\nDistance metrics\n" + "-" * 30)
        print("Absolute number of overlapping cells:", overlap, "/", self.domain_size)
        print("Relative area similarity:", overlap / self.domain_size)


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("-distribution1", "--d1", type=str)
    parser.add_argument("-distribution2", "--d2", type=str)
    parser.add_argument("-verbose", action="store_true")
    args = parser.parse_args()
    return args


def load_distr(distr_file, verbose=True):
    print("\nLoading distribution", distr_file) if verbose else None
    distribution = Distribution(distr_file)
    with open(distr_file, "r") as f:
        j = 0
        for i, line in enumerate(f.readlines()):
            line = line.replace("\n", "").replace(" ", "")
            if i == 0:
                continue
            elif i == 1:
                print("dim x: ", line) if verbose else None
                dim_x = int(line)
            elif i == 2:
                dim_y = int(line)
                print("dim y: ", line) if verbose else None
                distribution.create_array(dim_x, dim_y)
            else:
                for char in line:
                    distribution[int(j / dim_y), j % dim_y] = int(char)
                    j += 1

    distribution.compute_stats()
    print(f"Finished loading distribution {distr_file}") if verbose else None
    print(distribution) if verbose else None

    return distribution


def main():
    BASELINE_MORPHOLOGY = (
        "E:/Development/FESSGA/data/EVOMA/trex_100_elements/distribution2d.dens"
    )

    # Setup
    args = parse_args()
    verbose = args.verbose
    print("Distribution1 file: ", args.d1, "\n") if verbose else None
    print("Distribution2 file: ", args.d2, "\n") if verbose and args.d2 else None

    # Load distributions
    distr1 = load_distr(args.d1, verbose=verbose)
    if args.d2:
        distr2 = load_distr(args.d2, verbose=verbose)
    else:
        distr2 = load_distr(BASELINE_MORPHOLOGY, verbose=verbose)

    # Print stats
    distr1.print_stats()
    distr2.print_stats() if args.d2 else None

    # Get similarity measures
    distr1.compare(distr2)


if __name__ == "__main__":
    main()
