
# coding: utf-8

from __future__ import unicode_literals
"""
Class for analyzing workfunction by reading LOCPOT, vasprun.xml
"""
__author__ = "Zhe Luo"
__version__ = "0.0"
__maintainer__ = "Zhe Luo"
__email__ = "luozhe0629@sjtu.edu.cn"
__status__ = "Development"
__date__ = "Aug 16, 2016"


import warnings
import os
import numpy as np
from pymatgen.io.vasp import Vasprun, Locpot
from pymatgen import structure
from monty.json import MSONable
from pymatgen.core import PeriodicSite
from copy import deepcopy


# Check the Periodic Boundary Conditions, Make Sure that All the Frac Coords are in the range of (0,1)
def _wrap_frac_coords(frac_coords):
    t_frac = deepcopy(frac_coords)
    for i,coord in enumerate(t_frac):
        for j,x in enumerate(coord):
            while t_frac[i][j] < 0:
                t_frac[i][j] += 1
            while t_frac[i][j] > 1:
                t_frac[i][j] -= 1
    return t_frac

# Parse the structure to get the direction 
# ***********
# For convinience the structure of the cell must satisfy vacuum-layer-vacuum,layer-vacuum or vacuum-layer. Layer-vacuum-layer cannot be corrected analyzed and the progrom will NOT CHECK IT FOR YOU.
# 
# We strongly recommend that the direction along the vacuum region be perpendicular to the surface
# 
# ***********

def _parse_structure(struc, v_tolerance = 5, layer_tol = 0.5):
    
    """
    Args: 
        struc: Structure Object from pymatgen
        v_tolerance: If the layer distance is greater than 5 Ang (by Default), it will be regards as vacuum region
        layer_tol: If the distance (in the direction normal to the surface) between an atom and its next nearest 
            neighbor is larger than 0.5 than we assume that they are not in the same layer. 
    Return:
        index: the lattice parameter index where vacuum region is build
        proj: the cos value of the angle between the lattice parmameter and that normal to the surface
        layers: the dict storing the each layer.
    """
    
    #First Find the Vacuum Region
    lat = struc.lattice  
    mat = lat.matrix
    abc = tuple((np.abs((lat.a, lat.b, lat.c))))
    v_mat = [np.cross(mat[1],mat[2]),np.cross(mat[0],mat[2]),np.cross(mat[0],mat[1])]
    proj_cos = list()
    for i,vector in enumerate(v_mat):
        proj_cos.append(np.dot(vector,mat[i])/np.linalg.norm(vector)/np.linalg.norm(mat[i]))           
    frac_coords = struc.frac_coords
    frac_coords = _wrap_frac_coords(frac_coords)   
    max_frac = frac_coords.max(axis = 0)
    min_frac = frac_coords.min(axis = 0)  
    max_gap = (1- (max_frac - min_frac))*(abc)
    max_gap = abs(max_gap*proj_cos)

    
    flags = (max_gap > v_tolerance).tolist()
    flags_sum = sum(flags)
    
    if flags_sum == 1:
        for index,flag in enumerate(flags):
            if flag:break            
    #If only one vacuum region is found, which is reasonable go on parse the structure            
        proj = proj_cos[index]
        v_frac = (frac_coords).transpose()[index]
        order_v_frac = list(zip(range(len(v_frac)),v_frac))
        orderred_v_frac = sorted(order_v_frac, key = lambda couple:couple[1])
        order,v_frac = zip(*orderred_v_frac)
        
        
        gaps = np.abs(np.ediff1d(np.abs(np.array(v_frac)*abc[index])))
        layer_index = 0
        layers = {0:[order[0]]}
        for i,gap in enumerate(gaps):
            if (gap < layer_tol): 
                layers[layer_index].append(order[i+1])
            else:
                layer_index += 1
                layers[layer_index] = [(order[i+1])]
        return index, proj, layers     
        
    elif flags_sum>1:
        print "Multiple Vacuum Regions are Found. "
        print "Please Make Sure the Structure is Slab rather than Quantum Dots or Nanowire and the Tolerance is proper"
        pass
    
    else:
        print "No Vacuum Region, Please Check the Tolerance or Structure"
        pass
    


# Parse the sorted_structure and return the information of each layers(atoms type, coordinate)
class SingleLayer(MSONable):
    
    """
    Class storing the information of a single layer
    Args: Collectioon of Sites, List type or tuple
    
    """
    
    def __init__(self, sites, vdir):
    
        #Collect the coords information in the direction we care
        lat = sites[0].lattice.abc
        mat = sites[0].lattice.matrix
        v_mat = [np.cross(mat[1],mat[2]),np.cross(mat[0],mat[2]),np.cross(mat[0],mat[1])]
       
        v_frac_coords = list()  
        for site in sites:
            v_frac_coords.append(site.frac_coords[vdir])
        for i,z in enumerate(v_frac_coords):
            while v_frac_coords[i] < 0:
                v_frac_coords[i] += 1
            while v_frac_coords[i] > 1:
                v_frac_coords[i] -= 1
        v_frac_coords = np.array(v_frac_coords)
        v_cart_coords = v_frac_coords*lat[vdir]
    
        #Collect the number of atoms per species
        atoms = dict()
        for site in sites:
            specie = site.specie
            if not (specie in atoms.keys()):
                atoms[specie] = 1
            else:
                atoms[specie] += 1
        
        proj_cos = list()
        for i,vector in enumerate(v_mat):
            proj_cos.append(np.dot(vector,mat[i])/np.linalg.norm(vector)/np.linalg.norm(mat[i]))
        proj = abs(proj_cos[vdir])
        
        self._v_frac_coords = v_frac_coords
        self._v_cart_coords = v_cart_coords
        self._atoms = atoms
        self._proj = proj
        
    @property
    def v_frac_coords(self):
        return self._v_frac_coords
        
    @property
    def v_cart_coords(self):
        return self._v_cart_coords
        
    @property
    def layer_frac_loc(self):
        return np.mean(self.v_frac_coords)
        
    @property
    def layer_cart_loc(self):
        return np.mean(self.v_cart_coords)*self._proj
        
    @property
    def cart_std(self):
        return np.std(self.v_cart_coords)*self._proj
        
    @property
    def atoms(self):
        return self._atoms
        
    @property
    def total_number(self):
        return len(sites)
        
    @property
    def max_width(self):
        return (max(self.v_cart_coords) - min(self.v_cart_coords))*self._proj
        



class Layers(MSONable):
    
    def __init__(self, struc, v_tol = 5, layer_tol = 0.5):
        self._struc = deepcopy(struc)
        self._vdir, self._proj, self._layers = _parse_structure(struc, v_tol, layer_tol)
        
        
    @property
    def structure(self):
        return self._struc
    
    @property
    def layer_dict(self):
        return self._layers
    
    @property
    def layer_number(self):
        return len(self._layers)
    
    def layer_sites(self, layer_index):
        site_collection = list()
        for site_index in self.layer_dict[layer_index]:
            site_collection.append(struc[site_index])
        return site_collection

    @property
    def layer_frac_loc(self):
        loc = list()
        l = self.layer_number
        for i in range(l-1):
            loc.append(np.mean(SingleLayer(self.layer_sites(i), self._vdir).v_frac_coords))
        return loc
    
    def distance(self, i, j):
        ilayer = SingleLayer(self.layer_sites(i), self._vdir)
        jlayer = SingleLayer(self.layer_sites(j), self._vdir)
        lat = self.structure.lattice
        p = abs(((lat.a,lat.b,lat.c)[self._vdir]))
        dis = abs(ilayer.layer_cart_loc - jlayer.layer_cart_loc)
        return min(dis, (p-dis))
        
    @property
    def distances(self):
        dis = list()
        l = self.layer_number
        for i in range(l-1):
            dis.append(self.distance(i,i+1))
        return dis
    
    def layer_info(self, layer_index):
        slayer = SingleLayer(self.layer_sites(layer_index), self._vdir)
        return {"atoms":slayer.atoms,"loc":slayer.layer_cart_loc,"std":slayer.cart_std,"max_width":slayer.max_width}
        
    def move_layers(self, layer_indexes, vector):
        """
        Move all the sites in the same layer
        
        Args:
            layer_indexes: Sequence of the layers to move
            vector: Fractional Vector to move
        
        """
        for layer_index in layer_indexes:
            self._struc.translate_sites(self.layer_dict[layer_index],vector)
            
    def remove_layers(self, layer_indexes):
        """
        Remove all the sites in the same layer
        
        Args:
            layer_indexes: Sequence of the layers to delete
        """
        for layer_index in layer_indexes:
            self._struc.remove_sites(self.layer_dict[layer_index])
  
                
                

