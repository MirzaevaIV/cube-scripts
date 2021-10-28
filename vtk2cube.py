#!/usr/bin/python3

import numpy as np
from vtk import vtkStructuredPointsReader
from vtk.util import numpy_support as VN

def readVTK(fname):
    data = {}
    reader = vtkStructuredPointsReader()
    reader.SetFileName(fname)
    reader.Update()
    data0 = reader.GetOutput()
    data['N_steps'] = data0.GetDimensions()
    data['Steps'] = data0.GetSpacing()
    data['Origin'] = data0.GetOrigin()
    vox = VN.vtk_to_numpy(data0.GetPointData().GetScalars())
    vox = vox.reshape(data['N_steps'], order='F')
    vox = vox.flatten('C')    
    data['voxels'] = vox
    return data

def readXYZ(fname):
    coords = []
    with open(fname, 'r') as xyz:
        coords.append(xyz.readline().strip().split())
        N_at = int(coords[0][0])
        xyz.readline()
        for i in range(N_at):
           coords.append(xyz.readline().strip().split())
    return coords

def writeCube(fname, data, xyz, l):
    with open(fname, 'w') as f:
        f.write("# converted from vtk \n")
        f.write("# converted from vtk \n")
        entry = "%8s %12s %12s %12s \n" % (xyz[0][0], data['Origin'][0], data['Origin'][1], data['Origin'][2])
        f.write(entry)
        entry = "%8s %12.8f %12.8f %12.8f \n" % (data['N_steps'][0], data['Steps'][0], 0.0, 0.0)
        f.write(entry)
        entry = "%8s %12.8f %12.8f %12.8f \n" % (data['N_steps'][1], 0.0, data['Steps'][1], 0.0)
        f.write(entry)
        entry = "%8s %12.8f %12.8f %12.8f \n" % (data['N_steps'][2], 0.0, 0.0, data['Steps'][2])
        f.write(entry)
        for i in range(int(xyz[0][0])):
            el = l.index(xyz[i+1][0])+1
            x = float(xyz[i+1][1])/0.52917720859
            y = float(xyz[i+1][2])/0.52917720859
            z = float(xyz[i+1][3])/0.52917720859
            entry =  "%8s %10.6f %12.8f %12.8f %12.8f \n" % (el, 0.0, x, y, z)
            f.write(entry)
        entry = ""
        k = 0
        for i in range(data['N_steps'][0]*data['N_steps'][1]*data['N_steps'][2]):
            ent = "%16.12f" % float(data['voxels'][i])
            entry += ent
            k += 1
            if k == 9:
                entry += "\n"
                f.write(entry)
                entry = ""
                k = 0
        entry +="\n\n"
        f.write(entry)

if __name__ == '__main__':
   
    el_list = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si',
 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',
 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo',
 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba',
 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
 'Po', 'At', 'Rn', 'Fr',  'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm',
 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']
    vtkfile = input("Enter vtk file name: ")
    coordsfile = input("Enter atomic coordinates file name: ")
    cubefile = input("Enter resulting cube file name: ")
    cubedata = readVTK(vtkfile)
    coords = readXYZ(coordsfile)
    writeCube(cubefile, cubedata, coords, el_list)


