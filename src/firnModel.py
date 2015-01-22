'''
Created on Dec 10, 2014

@author: huongnvo
'''
import Tkinter
import tkFileDialog
from Tkinter import Menu
from Tkinter import *
import os
import csv
import json
import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.sparse.linalg as splin
from scipy.sparse import spdiags
from scipy.sparse.linalg import lsqr
from scipy.integrate import cumtrapz
import math
import ModelParameters.Gasses as MPG
import ModelParameters.Sites as MPS
import ModelParameters.Plotting as plots
import ModelParameters.Diffusivity as MPD
import ModelParameters.density as MPRHO
import time

class firnModel(Tkinter.Tk):
    def __init__(self, parent):
        Tkinter.Tk.__init__(self, parent)
        self.grid()
        self.parent = parent
        self.initialize()     
    
    def initialize(self):
        pass