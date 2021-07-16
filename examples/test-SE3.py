
from openravepy import *
from numpy import *

from toppso3.SE3RRT import *
from toppso3 import Utils
from toppso3 import lie

import time

import TOPP
from TOPP import TOPPpy
from TOPP import Trajectory
from pylab import *
import scipy.optimize
from mpl_toolkits.mplot3d import Axes3D

import pdb

ion()



env = Environment()
# This model was downloaded from http://nasa3d.arc.nasa.gov/models/printable
env.Load("../MESSENGER/messengerWithEnvSE3.xml")
env.SetViewer('qtcoin')

robot = env.GetBodies()[0]

