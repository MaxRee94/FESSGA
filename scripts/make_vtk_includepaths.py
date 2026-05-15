basedir = r"F:\Development\FESSGA\build\_deps\vtk\VTK-8.2.0\build\lib\Release"
import os
from glob import glob

libs = []
for path in glob(basedir + "/*"):
    if path.endswith(".lib"):
        libs.append(os.path.basename(path))

total = ";".join(libs)
print(total)
exit(0)

includepaths = []
for _dir in glob(basedir + "/*"):
    if _dir.endswith(".vs") or _dir.endswith("x64") or _dir.endswith("lib") or _dir.endswith("CMakeFiles") or _dir.endswith("Release"):
        continue
    for _subdir in glob(_dir + "/*"):
        _subdir = _subdir.replace("\\", "/").replace("F:/Development/FESSGA/build/", "")
        if "." in _subdir.split("/")[-1]:
            continue
        includepaths.append(_subdir)

total = ";".join(includepaths)
print(total)



