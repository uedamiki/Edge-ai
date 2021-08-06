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

class Config():
    """Attributes:
         q   -- quaternion vector
         qs  -- angular velocity

         qt  -- translation vector
         qts -- translational velocity
    """
    def __init__(self, q, qt, qs = None, qts = None, qss = None, qtss = None):
        self.q = q
        if (qs == None):
            self.qs = zeros(3)
        else:
            self.qs = qs

        self.qt = qt
        if (qts == None):
            self.qts = zeros(3)
        else:
            self.qts = qts


class Vertex():
    """Attributes:
         config     -- stores a Config obj
         parent     -- the parent for FW vertex, the child for BW vertex
         trajstring -- a trajectory fro