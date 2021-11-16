#RRT implementation for reorientation with collision-free
from openravepy import *
from pylab import *

import time
import string
import numpy as np
import copy
import random
import os

import lie

import Utils
import Heap

import TOPP
from TOPP import TOPPpy
from TOPP import Trajectory



# global variables for RRTPlanners
FW = 0
BW = 1
REACHED = 0
ADVANCED = 1
TRAPPED = 2

# gobal variables for collision checking 
INCOLLISION = -1
OK = 1

# def vectfromquat
# def eulerfromquat
# def quatfromuler
# def quatfromvect

class Config():
    """Attributes:
         q   -- quaternion vector
         qs  -- angular velocity
    """
    def __init__(self, q, qs = None, qss = None):
        self.q = q
        if (qs == None):
            self.qs = zeros(3)
        else:
            self.qs = qs


class Vertex():
    """Attributes:
         config     -- s