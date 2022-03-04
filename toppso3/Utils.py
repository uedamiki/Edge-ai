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
        d = 1
    else: d = d
    return sqrt(c*(SO3Distance(R0,R1)**2) + d*(R3Distance(b0,b1)**2))
    
################## interpolate translation ####################################
def TrajString3rdDegree(q_beg, q_end, qs_beg, qs_end, duration):
    trajectorystring = ''
    ndof = len(q_beg)
    
    trajectorystring += "%f\n%d"%(duration, ndof)

    for k in range(ndof):
        a, b, c, d = Utilities.Interpolate3rdDegree(q_beg[k], q_end[k], qs_beg[k], qs_end[k], duration)
        trajectorystring += "\n%f %f %f %f"%(d, c, b, a)
    return trajectorystring

#################### SE3 traj ##################################################

def SE3TrajFromTransandSO3(transtraj, rtraj): # same chunk.duration
    #return duration-dimension-trans polynomial- rot polynomial
    if len(transtraj.chunkslist) != len(rtraj.chunkslist):
        print 'error'
        return 0
    clist = []
    for c in transtraj.chunkslist:
        plist = []
        for i in range(3):
            plist.append(c.polynomialsvector[i])
        for i in range(3):
            rc = rtraj.chunkslist[len(clist)]
            plist.append(rc.polynomialsvector[i])
        chunk = Trajectory.Chunk(c.duration, plist)
        clist.append(chunk)
    return Trajectory.PiecewisePolynomialTrajectory(clist)

###################### Decompose SE3 traj to ROT and Trans traj ################
def TransRotTrajFromSE3Traj(SE3traj):
    transclist = []
    rclist = []
    for c in SE3traj.chunkslist:
        transchunk = Trajectory.Chunk(c.duration, c.polynomialsvector[:3])
        transclist.append(transchunk)
        rchunk = Trajectory.Chunk(c.duration, c.polynomialsvector[3:])
        rclist.append(rchunk)
    transtraj = Trajectory.PiecewisePolynomialTrajectory(transclist)
    rtraj = Trajectory.PiecewisePolynomialTrajectory(rclist)
    return transtraj, rtraj

##########################SE3 constraint ########################################
def ComputeSE3Constraints(SE3traj, taumax, fmax, discrtimestep, I = None, m = None):
    ndiscrsteps = int((SE3traj.duration + 1e-10) / discrtimestep) + 1
    a = zeros((ndiscrsteps,12))
    b = zeros((ndiscrsteps,12))
    c = zeros((ndiscrsteps,12))
    transtraj, rtraj = TransRotTrajFromSE3Traj(SE3traj)
    for i in range(ndiscrsteps):
        #rotconstraints
        t = i * discrtimestep
        r = rtraj.Eval(t)
        rd = 