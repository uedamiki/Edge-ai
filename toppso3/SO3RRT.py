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
         config     -- stores a Config obj
         parent     -- the parent for FW vertex, the child for BW vertex
         trajstring -- a trajectory from its parent (or child)
         # sdmin      -- minimum reachable sd
         # sdmax      -- maximum reachable sd
         level      -- its level from the root of the tree (0 for the root)
         # drawn      -- True if this vertex has been plotted via Vertex::Plot
    """
    def __init__(self, config, vertextype = FW):
        self.config = config
        self.vertextype = vertextype
        self.parent = None # to be assigned when added to a tree
        self.traj = '' # to be assigned when added to a tree
        self.level = 0
        # self.drawn = False 


class Tree():
    """Attributes:
         verticeslist -- stores all vertices added to the tree
         treetype     -- FW or BW    
    """
    def __init__(self, treetype = FW, vroot = None):
        if (vroot == None):
            self.verticeslist = []
        else:
            self.verticeslist = [vroot]
        self.treetype = treetype

    def __len__(self):
        return len(self.verticeslist)

    def __getitem__(self, index):
        return self.verticeslist[index]        
                    
    def AddVertex(self, parent, traj, vnew):
        vnew.parent = parent
        vnew.traj = traj
        self.verticeslist.append(vnew)

    def GenTrajList(self):
        trajlist = []
        if (self.treetype == FW):
            vertex = self.verticeslist[-1]
            parent = vertex.parent
            while (vertex.parent != None):
                trajlist.append(vertex.traj)
                vertex = parent
                if (vertex.parent != None):
                    parent = vertex.parent
            trajlist = trajlist[::-1]
        else:
            vertex = self.verticeslist[-1]
            while (vertex.parent != None):
                trajlist.append(vertex.traj)
                if (vertex.parent != None):
                    vertex = vertex.parent
        return trajlist
    
    def GenRotationMatList(self):
        RotationMatList = []
        if (self.treetype == FW):
            vertex = self.verticeslist[-1]
            RotationMatList.append(rotationMatrixFromQuat(vertex.config.q))
            parent = vertex.parent
            while (vertex.parent != None):
                RotationMatList.append(rotationMatrixFromQuat(parent.config.q))
                vertex = parent
                if (vertex.parent != None):
                    parent = vertex.parent
            RotationMatList =  RotationMatList[::-1]
        else:
            vertex = self.verticeslist[-1]
            RotationMatList.append(rotationMatrixFromQuat(vertex.config.q))                       
            while (vertex.parent != None):
                RotationMatList.append(rotationMatrixFromQuat(vertex.parent.config.q))
                if (vertex.parent != None):
                    vertex = vertex.parent
        return RotationMatList


class RRTPlanner():
    """Base class for RRT planners"""
    REACHED = 0
    ADVANCED = 1
    TRAPPED = 2

    def __init__(self, vertex_start, vertex_goal, robot):
        """Initialize a planner. RRTPlanner always has two trees. For a unidirectional planner, 
        the treeend will not be extended and always has only one vertex, vertex_goal.        
        """        
        # np.random.seed(np.random.randint(0, 10))
        ## need more unpredictable sequence than that generated from np.random
        self.RANDOM_NUMBER_GENERATOR = random.SystemRandom()
        
        self.treestart = Tree(FW, vertex_start)
        self.treeend = Tree(BW, vertex_goal)
        self.connectingtraj = []
        self.runningtime = 0.0
        self.nn = -1
        self.iterations = 0
        self.result = False
        
        # DEFAULT PARAMETERS
        self.STEPSIZE = 0.8
        self.INTERPOLATIONDURATION = 0.5
        # self.PLOT = False
        
        #Openrave paras
        self.robot = robot
        
        self.discrtimestep = 1e-2 ## for collision checking, etc.

    def __str__(self):
        ret = "Total running time :" + str(self.runningtime) + "sec.\n"
        ret += "Total number of iterations :" + str(self.iterations)
        return ret

    def RandomConfig(self):
        """RandomConfig samples a random configuration uniformly from the quaternion unit sphere in four dimensions."""
        q_rand = lie.RandomQuat()
        vellowerlimit = -5 ##
        velupperlimit = 5  ##
        qs_rand = np.array([1e-1,1e-1,1e-1 ])
        # for i in range(3):
        #     qs_rand[i] = self.RANDOM_NUMBER_GENERATOR.uniform(vellowerlimit,velupperlimit) 
        return Config(q_rand,qs_rand)

    def Extend(self, c_rand):
        if (np.mod(self.iterations - 1, 2) == FW):
            ## treestart is to be extended
            res = self.ExtendFW(c_rand)
        else:
            ## treeend is to be extended
            res = self.ExtendBW(c_rand)
        return res

    def ExtendFW(self, c_rand):
        nnindices = self.NearestNeighborIndices(c_rand, FW)
        for index in nnindices:
            v_near = self.treestart.verticeslist[index]
            q_beg = v_near.config.q
            qs_beg = v_near.config.qs
            
            ## check if c_rand is too far from vnear
            ## if the new ramdonly-chose node is close, it's safer . Or in another words, the interpolated path will have more chances that it won't collide with the obstacles
            delta = self.Distance(v_near.config, c_rand)
            if (delta <= self.STEPSIZE):
                q_end = c_rand.q
                STATUS = REACHED
            else:
                q_end = q_beg + self.STEPSIZE*(c_rand.q - q_beg)/np.sqrt(delta)
                q_end /= np.linalg.norm(q_end)
                STATUS = ADVANCED
            qs_end = c_rand.qs
            c_new = Config(q_end, qs_end)

            ## check feasibility of c_new
            if (not self.IsFeasibleConfig(c_new)):
                # print "status : TRAPPED (infeasible configuration)"
                STATUS = TRAPPED
                continue                        
            ## interpolate a trajectory
            #trajectory = lie.InterpolateSO3ZeroOmega(rotationMatrixFromQuat(q_beg),rotationMatrixFromQuat(q_end),self.INTERPOLATIONDURATION)
            trajectory = lie.InterpolateSO3(rotationMatrixFromQuat(q_beg),rotationMatrixFromQuat(q_end),qs_beg,qs_end,self.INTERPOLATIONDURATION)
            ## check feasibility ( collision checking for the trajectory)
            result = self.IsFeasibleTrajectory(trajectory, q_beg, FW)
            if (result[0] == OK):
                  ## extension is now successful
                v_new = Vertex(c_new, FW)
                # v_new.sdmin = result[1]
                # v_new.sdmax = result[2]
                v_new.level = v_near.level + 1
                self.treestart.AddVertex(v_near, trajectory, v_new)
                return STATUS
            else:
                STATUS = TRAPPED  #trajecory doesnt satify the collision-free constraint
        return STATUS

    def ExtendBW(self, c_rand):
        # Implement NearestneiborIndices return the list of nodes in order of increasing distance
        nnindices = self.NearestNeighborIndices(c_rand, BW)
        for index in nnindices:
            v_near = self.treeend.verticeslist[index]
            q_end = v_near.config.q
            qs_end = v_near.config.qs

            ## check if c_rand is too far from vnear
            ## if the new ramdonly-chose node is close, it's safer . Or in another words, the interpolated path will have more chances that it won't collide with the obstacles
            delta = self.Distance(v_near.config, c_rand)
            if (delta <= self.STEPSIZE):
                q_beg = c_rand.q
                STATUS = REACHED
            else:
                q_beg = q_end + self.STEPSIZE*(c_rand.q - q_end)/np.sqrt(delta)
                q_beg /= np.linalg.norm(q_beg)
                STATUS = ADVANCED
            qs_beg = c_rand.qs
            c_new = Config(q_beg, qs_beg)
            
            ## check feasibility of c_new
            if (not self.IsFeasibleConfig(c_new)):
                # print "status : TRAPPED (infeasible configuration)"
                STATUS = TRAPPED
                continue            

            ## interpolate a trajectory
            #trajectory = lie.InterpolateSO3ZeroOmega(rotationMatrixFromQuat(q_beg),rotationMatrixFromQuat(q_end),self.INTERPOLATIONDURATION)
            trajectory = lie.InterpolateSO3(rotationMatrixFromQuat(q_beg),rotationMatrixFromQuat(q_end),qs_beg,qs_end,self.INTERPOLATIONDURATION)
            ## check feasibility ( collision checking for the trajectory)
            result = self.IsFeasibleTrajectory(trajectory, q_beg, BW)
            if (result[0] == OK):
                ## extension is now successful
                v_new = Vertex(c_new, BW)
                # v_new.sdmin = result[1]
                # v_new.sdmax = result[2]
                v_new.level = v_near.level + 1
                self.treeend.AddVertex(v_near, trajectory, v_new)
                return STATUS
            else:
                STATUS = TRAPPED  #trajecory doesnt satify the collision-free constraint
        return STATUS


    def Connect(self):
        if (np.mod(self.iterations - 1, 2) == FW):
            ## treestart has just been extended
            res = self.ConnectBW()
        else:
            ## treeend has just been extended
            res = self.ConnectFW()
        return res
        
    def ConnectFW(self):
        v_test = self.treeend.verticeslist[-1]
        nnindices = self.NearestNeighborIndices(v_test.config, FW)
        for index in nnindices:
            v_near = self.treestart.verticeslist[index]
            
            q_beg = v_near.config.q
            qs_beg = v_near.config.qs
            
            q_end = v_test.config.q
            qs_end = v_test.config.qs
            
             ## interpolate a trajectory
            #trajectory = lie.InterpolateSO3ZeroOmega(rotationMatrixFromQuat(q_beg),rotationMatrixFromQuat(q_end),self.INTERPOLATIONDURATION)
            trajectory = lie.InterpolateSO3(rotationMatrixFromQuat(q_beg),rotationMatrixFromQuat(q_end),qs_beg,qs_end,self.INTERPOLATIONDURATION)
             ## check feasibility ( collision checking for the trajectory)
            result = self.IsFeasibleTrajectory(trajectory, q_beg, FW)
            if (result[0] == 1):
                 ## conection is now successful
                self.treestart.verticeslist.append(v_near)
                self.connectingtraj = trajectory
                return REACHED
        return TRAPPED

    def ConnectBW(self):
        v_test = self.treestart.verticeslist[-1]
      