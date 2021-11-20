#!/usr/bin/python3

import numpy as np
import argparse
import sys

el_list = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si',
 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',
 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo',
 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba',
 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
 'Po', 'At', 'Rn', 'Fr',  'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm',
 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']
a_B = 1/0.52917720859


class cubeParser():
    
    def readCube(fname):
        """ 
        Read file of Gaussian cube format into dictionary    
        params:  (fname: filename of cube file)       
        returns: (data: dictionary)
        """        
        data = {}
        with open(fname, 'r') as cube:
            cube.readline(); cube.readline()  # ignore two first lines
            data['Origin'] = list(np.float_(cube.readline().strip().split()))
            data['N_atoms'] = int(data['Origin'].pop(0))
            data['da'] = list(np.float_(cube.readline().strip().split())) # shift along lattice vector a
            data['N_asteps'] = int(data['da'].pop(0))
            data['db'] = list(np.float_(cube.readline().strip().split())) # shift along lattice vector b
            data['N_bsteps'] =int(data['db'].pop(0))
            data['dc'] = list(np.float_(cube.readline().strip().split())) # shift along lattice vector c
            data['N_csteps'] = int(data['dc'].pop(0))
            data['atoms'] = [list(np.float_(cube.readline().strip().split())) for i in range(data['N_atoms'])]
            Nvox = data['N_asteps']*data['N_bsteps']*data['N_csteps']
            data['voxels'] = np.zeros((Nvox))
            i = 0 # count voxels
            for line in cube:
                for val in line.strip().split():
                    data['voxels'][i] = float(val)
                    i +=1                               
        print('Done reading cube file.')
        data['voxels'] = np.reshape(data['voxels'], (data['N_asteps'],data['N_bsteps'],data['N_csteps']))        
        return data
    
def writeCube(fname, data, xyz):
    '''Takes dictionary with vtk data and list of atomic coordinates and writes a cube file'''
    with open(fname, 'w') as f:
        f.write("# converted from vtk \n")
        f.write("# We wanted the best, but it turned out like always \n")
        entry = "%8s %12s %12s %12s \n" % (xyz[0][0], data['Origin'][0], data['Origin'][1], data['Origin'][2])
        f.write(entry)
        entry = "%8s %12.8f %12.8f %12.8f \n" % (data['N_steps'][0], data['Steps'][0], 0.0, 0.0)
        f.write(entry)
        entry = "%8s %12.8f %12.8f %12.8f \n" % (data['N_steps'][1], 0.0, data['Steps'][1], 0.0)
        f.write(entry)
        entry = "%8s %12.8f %12.8f %12.8f \n" % (data['N_steps'][2], 0.0, 0.0, data['Steps'][2])
        f.write(entry)
        for i in range(int(xyz[0][0])):
            el = el_list.index(xyz[i+1][0])+1
            x = float(xyz[i+1][1])/0.52917720859
            y = float(xyz[i+1][2])/0.52917720859
            z = float(xyz[i+1][3])/0.52917720859
            entry =  "%8s %10.6f %12.8f %12.8f %12.8f \n" % (el, 0.0, x, y, z)
            f.write(entry)
        entry = ""
        k = 0
        for i in range(data['N_steps'][0]*data['N_steps'][1]*data['N_steps'][2]):
            ent = "%.12f " % float(data['voxels'][i])
            entry += ent
            k += 1
            if k == 9:
                entry += "\n"
                f.write(entry)
                entry = ""
                k = 0
        entry +="\n\n"
        f.write(entry)


def multiplicate(data, Nx, Ny, Nz):
    m_data = {}
    m_data['voxels'] = np.zeros((data['N_asteps']*Nx, data['N_bsteps']*Ny, data['N_csteps']*Nz))
    m_data['Steps'] = [data['da'][0], data['db'][1], data['dc'][2]]
    m_data['N_steps'] = [data['N_asteps']*Nx, data['N_bsteps']*Ny, data['N_csteps']*Nz]
    m_data['Origin'] = data['Origin']
    coords = []
    coords.append([data['N_atoms']])
    for i in range(coords[0][0]):
        el = el_list[int(data['atoms'][i][0])-1]
        coords.append([el, data['atoms'][i][2]/a_B, data['atoms'][i][3]/a_B, data['atoms'][i][4]/a_B])

    for i in range(Nx):
        for j in range(data['N_asteps']):
            for k in range(Ny):
                for l in range(data['N_bsteps']):
                    for m in range(Nz):
                        for n in range(data['N_csteps']):
                            m_data['voxels'][(data['N_asteps']*i)+j][(data['N_bsteps']*k)+l][(data['N_csteps']*m)+n] = data['voxels'][j][l][n]
    
    m_data['voxels'] = m_data['voxels'].flatten('C')
    return m_data, coords

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cube', help = "Cube file name")
    parser.add_argument('-o', '--output', help = "Output filename")
    parser.add_argument('-r', '--repite', nargs = 3, type = int, metavar = ('Nx','Ny', 'Nz'), \
         default = ('1', '1', '1'), help = "Number of repetitions. Nx - along the first translation vector, \
                        Ny - along the second translation vector, Nz - along the third translation vector. ")
    args = parser.parse_args()

    if args.cube == None:
        cubefile = input("Enter cube file name: ")
    else:
        cubefile = args.cube
   
    if args.output == None:
        resfile = input("Enter results file name: ")
    else:
        resfile = args.output

    a = cubeParser.readCube(cubefile)

    cubedata, coords = multiplicate(a, int(args.repite[0]), int(args.repite[1]), int(args.repite[2]))

    writeCube(resfile, cubedata, coords)
    
