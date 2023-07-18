import os


def main():
    boundconds = [
        "306 402",
    ]
    bound_conds_amended = []
    no_edges_per_boundcond = []
    for bounds in boundconds:
        _bounds = bounds.split(" ")
        _bounds_amended = []
        for i in range(int(len(_bounds) / 2)):
            start = int(_bounds[i * 2])
            end = int(_bounds[i * 2 + 1])
            diff = end - start
            for j in range(diff + 1):
                _bounds_amended.append(str(start + j))
        bound_conds_amended.append(" ".join(_bounds_amended))
        no_edges_per_boundcond.append(len(_bounds_amended))

    for no_edges, bounds in zip(no_edges_per_boundcond, bound_conds_amended):
        print("\nno edges =", no_edges, "\n", bounds)


if __name__ == "__main__":
    main()
