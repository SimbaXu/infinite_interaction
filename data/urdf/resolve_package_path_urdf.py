"""This script processes a ros urdf file so that all paths specified
using ros package:// syntax are resolved properly, allowing one to use
the urdf file with pybullet directly.

Usage
  python resolve_package_path_urdf.py [filename]

"""
import sys
import rospkg
import re

rospack = rospkg.RosPack()
filename = sys.argv[1]
match = re.match(r'([a-zA-Z_]*).urdf', filename, re.M)
if match is None:
    print("given file {:} is not an urdf file".format(filename))
filename1 = match.group(1)

with open(filename, 'r') as f:
    doc = f.read()


while True:
    # match package://[package_name]/
    match = re.search(r'package://([A-Za-z_]*)', doc, re.M)
    if match is None:
        break
    package_name = match.group(1)
    print("Found package: {:}".format(package_name))
    package_path = rospack.get_path(package_name)
    doc = doc.replace(match.group(), package_path)

filename_new = "{:}_res.urdf".format(filename1)
with open(filename_new, 'w') as f:
    f.write(doc)
    print("Wrote a new urdf file at: {:}".format(filename_new))
