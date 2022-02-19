from openravepy import *
from pylab import *

from numpy import *

import lie
import time 

import string
import StringIO

import TOPP
from TOPP import TOPPpy
from TOPP import TOPPbindings
from TOPP import Trajectory
from TOPP import Utilities

import pdb ########################
import matplotlib.pyplot as plt


def QuatDistance(quat0, quat1): 
    rotationweight = 1
    innerProduct = dot(quat0,quat1)
    quatDistance = rotationweight*(1-abs(innerProduct))
    return quatDistance

def SO3Distance(R0, R1): # bi-invariance
    return linalg.norm(lie.logvect(dot(R0.T,R1)))

def R3Distance(b0, b1):
    return linalg.norm(b1-b0)

def SE3Distance(X0, X1, c = None, d = None): # left invariance
    R0 = X0[:3,:3]
    R1 = X1[:3,:3]
    b0 = X0[:3,3]
    b1 = X1[:3,3]
    if (c == None):
        c = 1
    else: c = c
    if (d == None):
  