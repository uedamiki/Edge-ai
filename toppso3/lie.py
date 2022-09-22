# Interpolation in SO(3) following Park and Ravani
import time

import TOPP
from TOPP import Trajectory
import bisect
from pylab import *
from numpy import *
import pdb
import matplotlib.pyplot as plt




Eps = zeros((3,3,3))
Eps[0,2,1] = -1
Eps[0,1,2] = 1
Eps[1,0,2] = -1
Eps[1,2,0] = 1
Eps[2,1,0] = -1
Eps[2,0,1] = 1


class LieTraj():
    def __init__(self, Rlist,trajlist):
        self.Rlist = Rlist
        self.trajlist = trajlist
        self.duration = 0
        self.trajcumulateddurationslist = []
        for t in trajlist:
            self.trajcumulatedduratio