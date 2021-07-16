
from openravepy import *
from numpy import *

from toppso3.SO3RRT import *
from toppso3 import Utils
from toppso3 import lie

import time
import TOPP
from TOPP import TOPPpy
from TOPP import TOPPbindings
from TOPP import Trajectory
from pylab import *
import scipy.optimize
from mpl_toolkits.mplot3d import Axes3D


ion()



env = Environment()
# This model was downloaded from http://nasa3d.arc.nasa.gov/models/printable
env.Load("../MESSENGER/messengerWithEnv.xml")
env.SetViewer('qtcoin')

robot = env.GetBodies()[0]


phi = pi

R0 = eye(3)
q0 = quatFromRotationMatrix(R0)
omega0 = zeros(3)

q1 = array([cos(phi/2.),0,0,sin(phi/2.)])
omega1 = zeros(3)

taumax = ones(3)
vmax = ones(3)
inertia = eye(3)
################################## BiRRT planner #################################

vertex_beg = Vertex(Config(q0,omega0), FW)
vertex_end = Vertex(Config(q1,omega1), BW)