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
from layers import Layers


# Extract the fermi level, final structure from vasprun.xml
def _parse_vasprun(filename="vasprun.xml"):
    vasprun = Vasprun(filename)
    return vasprun.efermi, vasprun.structures[-1]


# Parse the LOCPOT file to get the potential in the direction along vacuum region, perpendicular to the surfae

def _parse_locpot(vdir,filename="LOCPOT"):
    locpot = Locpot(filename)
    return locpot.get_average_along_axis(vdir)


class WFquickAnalyzer(MSONable):
    
    def __init__(self, struc, pot, efermi):
        self._layers = Layers(struc, v_tol, layer_tol)
        self._pot = pot
        self._efermi = efermi
        
    @staticmethod
    def from_file(vasprun = "vasprun.xml", locpot = "LOCPOT"):
        efermi, struc = _parse_vasprun(vasprun_file)
        layers = Layers(struc)
        pot = _parse_locpot(locpot_file)
        return WFquickAnalyzer(struc, pot, efermi)
    
    @staticmethod
    def from_dir(path):
        vasprun = os.path.join(path, "vasprun.xml")
        locpot = os.path.join(path, "LOCPOT")
        return WFquickAnalyzer.from_file(vasprun, locpot)
    

class WorkfunctionAnalyzer(MSONable):
    
    """
    Class for analyzing workfunction by reading LOCPOT, vasprun.xml
    
    ****IMPORTANT******
    For convinience the structure of the cell must satisfy vacuum-layer-vacuum
    ,layer-vacuum or vacuum-layer. Layer-vacuum-layer cannot be corrected
    analyzed and the progrom will NOT CHECK IT FOR YOU. 
    *******************
    
    Args:
        structure: final relaxed structure(last step of the relaxation, or from the static run)
        vdir = 
    
    Return:
        *************
        For accuarcy we would recommend you use the reference way to calculate work
        function rather than just by the fermi level and vacuum level if
        possible. For more information see the reference.( Ramprasad R, von
        Allmen P, Fonseca L R C. Contributions to the work function: A
        density-functional study of adsorbates at graphene ribbon edges[J].
        Physical Review B, 1999, 60(8): 6023.)
        *************
        vacuum_level: The vacuum level. By default We use the 1/4 of the vacuum
            length near the vacuum center to calculate vacuum potential. You can
            also specify the data you want to use by inputting the list of
            index of all data.
    """
    
    def __init__(self, struc, efermi, pot, v_tol = 5, layer_tol = 0.5):
        
        self._layers = Layers(struc, v_tol, layer_tol)
        self._pot = pot
        self._ref_bulk = ref_bulk
        self._efermi = efermi
        
    @staticmethod
    def from_file(, locpot = "LOCPOT", vasprun = "vasprun.xml", v_tol = 5, layer_tol = 0.5):
        efermi, struc = _parse_vasprun(vasprun)
        layers = Layers(struc, v_tol, layer_tol)
        pot = _parse_locpot(locpot)
        return WorkfunctionAnalyzer.(struc, efermi, pot, v_tol = 5, layer_tol = 0.5)
        
        
        

