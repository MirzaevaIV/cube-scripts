# cube-scripts

The project contains scripts for working with Gaussian cube files and related.

cube-find-voxel.py
--------------------
The script finds the values of the property in the set of given points. Points are read from file.

Requirements: numpy, argparse

usage: cube-find-voxel.py [-h] [-c CUBE] [-p POINTS] [-o OUTPUT] [-f | -a | -b]

optional arguments: \
  -h, --help   --         show this help message and exit \
  -c CUBE, --cube CUBE  --  Cube file name \
  -p POINTS, --points POINTS  --   Name of the file with the list of points coordinates \
  -o OUTPUT, --output OUTPUT  --   Output filename \
  -f, --fractional      points coordinates are given in fractional coordinates (default) \
  -a, --angstrom        points coordinates are given in angstroms \
  -b, --bohr            points coordinates are given in bohrs


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


multiplicate-cube.py
------------------------
usage: multiplicate-cube.py [-h] [-c CUBE] [-o OUTPUT] [-r Nx Ny Nz]

optional arguments: \
  -h, --help     show this help message and exit \
  -c CUBE, --cube CUBE  Cube file name \
  -o OUTPUT, --output OUTPUT    Output filename \
  -r Nx Ny Nz, --repite Nx Ny Nz      Number of repetitions. Nx - along the first translation vector, \
  Ny - along the second translation vector, Nz - along the third translation vector. 


Known problems
--------------------
vtk2cube.py doesn't work with anaconda python on windows \
multiplicate-cube.py works only with rectangular mesh

