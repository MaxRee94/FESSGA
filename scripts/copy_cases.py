import os
from glob import glob
import shutil

base = "./"

if "EVOMA" in os.getcwd():
    print("Detected EVOMA.")
    cases = list(glob("./best_solutions/*/constrained_pull_and_bite_0001.vtk"))
else:
    cases = list(glob("./iteration_*/constrained_pull_and_bite_0001.vtk"))

if not os.path.isdir("./cases"):
    os.makedirs("cases")

print("Copying cases...")
for f in cases:
    f_target = f.replace("\\", "/").split("/")
    iteration_folder = f_target.pop(-2)
    if "EVOMA" in os.getcwd():
        f_target.pop(-2)
    f_target.insert(-1, "cases")
    f_target = "/".join(f_target)
    f_target = f_target.replace(".vtk", " " + iteration_folder + ".vtk").replace(
        " ", "_"
    )
    print("file target: ", f_target)
    shutil.copy(f, f_target)

print("Super position copying finished.")
