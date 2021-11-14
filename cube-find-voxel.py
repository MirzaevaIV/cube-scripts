#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 12:36:38 2021

@author: dairdre
"""

import numpy as np

class cubeParser():
    
    def readCube(fname):
        """ 
        Read file of Gaussian cube format into dictionary    
        params:  (fname: filename of cube file)       
        returns: (data: dict)
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
    
    def getVoxel(vec, cube):
        """ 
        Get the value of the property in the point  
        params:  (vec: np.array of cartesian coordinates of the point in bohrs; 
                  cube: dictionary with property cube data)       
        returns: (val: the value of the property)
        """
        a = np.array(cube['da'])
        b = np.array(cube['db'])
        c = np.array(cube['dc'])
        Tmatrix = cubeParser.getTransformMatrix(a, b, c)
        internalCoords = np.dot(Tmatrix, vec-np.array(cube['Origin']))
        val = cube['voxels'][int(internalCoords[0]), int(internalCoords[1]), int(internalCoords[2])]
        return val

    def getVoxelInternal(vec, cube):
        """ 
        Get the value of the property in the point  
        params:  (vec: np.array of the internal coordinates of the point; 
                  cube: dictionary with property cube data)       
        returns: (val: the value of the property)
        """
        a = int(vec[0]*cube['N_asteps'])
        b = int(vec[1]*cube['N_bsteps'])
        c = int(vec[2]*cube['N_csteps'])
        val = cube['voxels'][a, b, c]
        return val
        
        
    def getTransformMatrix(vec1, vec2, vec3):
        """ 
        Defines lattice angles for lattice vectors given in cartesian coordinates  
        params:  (vec1: np.array of a coordinates;
                  vec2: np.array of b coordinates;
                  vec3: np.array of c coordinates)       
        returns: (Tmatrix: transformation matrix in np.array form)
        """
        a = np.linalg.norm(vec1)
        b = np.linalg.norm(vec2)
        c = np.linalg.norm(vec3)
        angles = cubeParser.getAngles(vec1, vec2, vec3)
        cosa = np.cos(angles[0])
        cosb = np.cos(angles[1])
        cosg = np.cos(angles[2])
        sing = np.sin(angles[2])
        volume = 1.0 - cosa**2.0 - cosb**2.0 - cosg**2.0 + 2.0 * cosa * cosb * cosg
        volume = np.sqrt(volume)
        Tmatrix = np.zeros((3, 3))
        Tmatrix[0, 0] = 1.0 / a
        Tmatrix[0, 1] = -cosg / (a * sing)
        Tmatrix[0, 2] = (cosa * cosg -  cosb) / (a * volume * sing)
        Tmatrix[1, 1] = 1.0 / (b * sing)
        Tmatrix[1, 2] = (cosb * cosg - cosa) / (b * volume * sing)
        Tmatrix[2, 2] = sing / (c * volume)
        return Tmatrix
        
        
    def getAngles(vec1, vec2, vec3):
        """ 
        Defines lattice angles for three lattice vectors given in cartesian coordinates  
        params:  (vec1: np.array of a lattice vector cartesian coordinates;
                  vec2: np.array of b lattice vector cartesian coordinates;
                  vec3: np.array of c lattice vector cartesian coordinates)       
        returns: (vecA: np.array of angles in  radians)
        """
        norm1 = np.linalg.norm(vec1)
        norm2 = np.linalg.norm(vec2)
        norm3 = np.linalg.norm(vec3)
        alpha = np.arccos(np.dot(vec2, vec3)/(norm2*norm3))
        beta = np.arccos(np.dot(vec3, vec1)/(norm3*norm1))
        gamma = np.arccos(np.dot(vec1, vec2)/(norm1*norm2))
        vecA = np.array([alpha, beta, gamma])        
        return vecA
    
def readPoints(fname):
    """ 
    Read text file the list of points (cartesian coordinates in bohrs)   
    params:  (fname: filename of the file with the list of points)       
    returns: (points: the list of lists, cartesian coordinates of the points in bohrs)
    """
    points = []
    with open(fname, 'r') as fpoints:
        for line in fpoints:
            points.append(list(np.float_(line.strip().split())))
    return points
    


if __name__ == '__main__':
    
    cubefile = input("Enter cube file name: ") #"RHO-man.cube"
    pointsfile = input("Enter points list file name: ")
    resfile = input("Enter results file name: ")
    plist = readPoints(pointsfile)
    a = cubeParser.readCube(cubefile)
    results = []
    for i in plist:
        results.append(cubeParser.getVoxelInternal(np.array(i), a))
    with open(resfile, 'w') as f:
        for item in results:
            f.write("%s \n" % item)
    
    
    
