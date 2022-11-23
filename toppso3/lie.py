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
            self.trajcumulateddurationslist.append(self.duration)
            self.duration += t.duration

    def FindTrajIndex(self, s):
        if s == 0:
            s = 1e-10
        i = bisect.bisect_left(self.trajcumulateddurationslist, s) - 1
        remainder = s - self.trajcumulateddurationslist[i]
        return i, remainder

    # Rotation
    def EvalRotation(self,s):
        i, remainder = self.FindTrajIndex(s)
        return(dot(self.Rlist[i],expmat(self.trajlist[i].Eval(remainder))))

    # Velocity in body frame
    def EvalOmega(self,s):
        i, remainder = self.FindTrajIndex(s)
        r = self.trajlist[i].Eval(remainder)
        rd = self.trajlist[i].Evald(remainder)
        return dot(Amat(r),rd)

    # Acceleration in body frame
    def EvalAlpha(self,s):
        i, remainder = self.FindTrajIndex(s)
        r = self.trajlist[i].Eval(remainder)
        rd = self.trajlist[i].Evald(remainder)
        rdd = self.trajlist[i].Evaldd(remainder)
        return dot(Bmat(r),rdd) + dot(rd,tensordot(Ctensor(r),rd,([2],[0])))

    # Torques
    def EvalTorques(self,s,I):
        i, remainder = self.FindTrajIndex(s)
        r = self.trajlist[i].Eval(remainder)
        rd = self.trajlist[i].Evald(remainder)
        rdd = self.trajlist[i].Evaldd(remainder)
        omega = dot(Amat(r),rd)
        alpha =  dot(Bmat(r),rdd) + dot(rd,tensordot(Ctensor(r),rd,([2],[0])))
        return dot(I,alpha) + cross(omega,dot(I,omega))

    
    def Plot(self,dt=0.01,figstart=0,vmax=[],accelmax=[],taumax=[],I=None):

        tvect = arange(0, self.duration + dt, dt)
        omegavect = array([self.EvalOmega(t) for t in tvect])
        figure(figstart)
        clf()
        
        plt.plot(tvect,omegavect[:,0],'--',label = '$\omega^1$',linewidth = 2)
        plt.plot(tvect,omegavect[:,1],'-.',label = '$\omega^2$',linewidth = 2)
        plt.plot(tvect,omegavect[:,2],'-',label ='$\omega^3$',linewidth = 2)
        plt.legend()
        for v in vmax:
            plt.plot([0, self.duration],[v, v], '-.',color = 'k')
        for v in vmax:
            plt.plot([0, self.duration],[-v, -v], '-.',color = 'k')
        ylabel('Angular velocities (rad/s)')
        xlabel('Time (s)')

        alphavect = array([self.EvalAlpha(t) for t in tvect])
        figure(figstart+1)
        clf()
        plt.plot(tvect,alphavect[:,0],'--',label = '$\dot \omega^1$',linewidth = 2)
        plt.plot(tvect,alphavect[:,1],'-.',label = '$\dot \omega^2$',linewidth = 2)
        plt.plot(tvect,alphavect[:,2],'-',label = '$\dot \omega^3$',linewidth = 2)      
        plt.legend()
        for a in accelmax:
            plt.plot([0, self.duration],[a, a], '-.',color = 'k')
        for a in accelmax:
            plt.plot([0, self.duration],[-a, -a], '-.',color = 'k')
        ylabel('Angular accelerations (rad/s^2)')
        xlabel('Time (s)')

        if not(I is None):
            torquesvect = array([self.EvalTorques(t,I) for t in tvect])
            figure(figstart+2)
            clf()
            plt.plot(tvect,torquesvect[:,0],'--',label = r'$\tau^1$',linewidth = 2)
            plt.plot(tvect,torquesvect[:,1],'-.',label = r'$\tau^2$',linewidth = 2)
            plt.plot(tvect,torquesvect[:,2],'-',label = r'$\tau^3$',linewidth = 2)
            plt.legend()
            
            for tau in taumax:
                plt.plot([0, self.duration],[tau, tau], '-.',color = 'k')
            for tau in taumax:
                plt.plot([0, self.duration],[-tau, -tau], '-.',color = 'k')
            ylabel('Torques (N.m)')
            xlabel('Time (s)')
        
def SplitTraj(Rlist,traj):
    trajlist = []
    chunkindex = 0
    clist = []
    for i in range(len(Rlist)-1):
        while chunkindex <  len(traj.chunkslist):
            chunkcur = traj.chunkslist[chunkindex]
            chunknext = traj.chunkslist[chunkindex+1]
            clist.append(chunkcur)
            chunkindex += 1
            if(norm(dot(Rlist[i],expmat(chunkcur.Eval(chunkcur.duration)))-dot(Rlist[i+1],expmat(chunknext.Eval(0)))))< 1e-8:
                trajlist.append(Trajectory.PiecewisePolynomialTrajectory(clist))
                clist = []
                break
    # Last traj
    clist = []
    while chunkindex < len(