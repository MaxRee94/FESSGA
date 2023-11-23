import os
from glob import glob
import shutil

base = "./"

if "EVOMA" in os.getcwd():
    print("Detected EVOMA.")
    superpositions = list(glob("./best_solutions/*/SuperPosition.vtk"))
else:
    superpositions = list(glob("./iteration_*/SuperPosition.vtk"))

if not os.path.isdir("./superpositions"):
    os.makedirs("superpositions")

print("Copying superpositions...")
for f in superpositions:
    f_target = f.replace("\\", "/").split("/")
    iteration_folder = f_target.pop(-2)
    if "EVOMA" in os.getcwd():
        f_target.pop(-2)
    f_target.insert(-1, "superpositions")
    f_target = "/".join(f_target)
    f_target = f_target.replace(".vtk", " " + iteration_folder + ".vtk").replace(
        " ", "_"
    )
    print("file target: ", f_target)
    shutil.copy(f, f_target)

print("Super position copying finished.")
