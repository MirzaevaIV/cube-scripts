# cube-scripts

The project contains scripts for working with Gaussian cube files and related.

cube-find-voxel.py
--------------------
The script finds the values of the property in the set of given points. Points are read from file.

Requirements: numpy, argparse

usage: cube-find-voxel.py [-h] [-c CUBE] [-p POINTS] [-o OUTPUT]

optional arguments: \
  -h, --help   --         show this help message and exit \
  -c CUBE, --cube CUBE  --  Cube file name \
  -p POINTS, --points POINTS  --   Name of the file with the list of points coordinates \
  -o OUTPUT, --output OUTPUT  --   Output filename


vtk2cube.py
--------------------
The script reads binary vtk file from adf and xyz file with atomic coordinates and writes a cube file.

Requirements: numpy, vtk, argparse

usage: vtk2cube.py [-h] [-v VTK] [-c COORDS] [-o OUTPUT]

optional arguments: \
  -h, --help     --       show this help message and exit \
  -v VTK, --vtk VTK  --   VTK file name \
  -c COORDS, --coords COORDS  --    XYZ atomic coordinates file name \
  -o OUTPUT, --output OUTPUT  --    Output cube file name \
  \
  -f, --fractional      points coordinates are given in fractional coordinates (default) \
  -a, --angstrom        points coordinates are given in angstroms \
  -b, --bohr            points coordinates are given in bohrs



Known problems
--------------------
vtk2cube.py doesn't work with anaconda python on windows

