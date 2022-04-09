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
        rd = rtraj.Evald(t)
        rdd = rtraj.Evaldd(t)
        nr = linalg.norm(r)
        nr2 = nr*nr
        nr3 = nr2*nr
        nr4 = nr3*nr
        nr5 = nr4*nr
        R = lie.skewfromvect(r)

        snr = sin(nr)
        cnr = cos(nr)
        rcrd = cross(r,rd)
        rdrd = dot(r,rd)

        Amat =  eye(3) - (1-cnr)/nr2*R + (nr-snr)/nr3*dot(R,R)
        C1 = (nr-snr)/nr3 * cross(rd,rcrd)
        C2 = -(2*cnr+nr*snr-2)/nr4 * rdrd*rcrd
        C3 = (3*snr-nr*cnr - 2*nr)/nr5 * rdrd*cross(r,rcrd)
        C = C1+C2+C3

        Ard = dot(Amat,rd)
        if I is None:            
            at = Ard
            bt = dot(Amat,rdd) + C
        else:
            at = dot(I,Ard)
            bt = dot(I,dot(Amat,rdd)) + dot(I,C) + cross(Ard,dot(I,Ard))
        
        a[i,3:6] = at
        a[i,9:12] = -at
        b[i,3:6] = bt
        b[i,9:12] = -bt
        c[i,3:6] = -taumax
        c[i,9:12] = -taumax

        #transconstraints
        td = transtraj.Evald(t)
        tdd = transtraj.Evaldd(t)
        if m is None:
            at = td
            bt = tdd

        a[i,:3] = at
        a[i,6:9] = -at
        b[i,:3] = bt
        b[i,6:9] = -bt
        c[i,:3] = -fmax
        c[i,6:9] = -fmax
    return a, b, c

######################## se3 traj collision checking ########################

def CheckCollisionSE3Traj( robot, transtraj, rtraj, R_beg,  checkcollisiontimestep = 1e-3):
    """CheckCollisionSE3Traj accepts a robot and trans, rot trajectory object as its inputs.
       (checkcollisiontimestep is set to 1e-3 as a default value)
       It returns True if any config along the traj is IN-COLLISION.
    """
    env = robot.GetEnv()
    for s in np.arange(0, transtraj.duration, checkcollisiontimestep):
        with robot:
            transformation = eye(4)
            transformation[0:3,0:3] = lie.EvalRotation(R_beg, rtraj, s)
            transformation[0:3,3] = transtraj.Eval(s)
            robot.SetTransform(transformation)           
            isincollision = (env.CheckCollision(robot, CollisionReport()))
            #print  "s =", s, " ", isincollision
            if (isincollision):
                return True
    with robot:
        robot.SetTransform(transformation)
        isincollision = (env.CheckCollision(robot, CollisionReport()))
        if (isincollision):
            return True
        else:
            return False


######################### SE3 shortcutting ##################################
def SE3Shortcut(robot, taumax, fmax, vmax, se3traj, Rlist, maxiter, expectedduration = -1,  meanduration = 0, upperlimit = -1, plotdura = None):
    if plotdura == 1:
        plt.axis([0, maxiter, 0, se3traj.duration])
        plt.ion()
        plt.show()
        ylabel('Trajectory duration (s)')
        xlabel('Iteration')

    t_sc_start = time.time()
    originalduration =  se3traj.duration
    #return shortcuted traj
    if upperlimit < 0:
        dur = se3traj.duration
        upperlimit = se3traj.duration
    else:
        dur = upperlimit
    attempt = 0

    ## for shortcutting
    integrationtimestep = 1e-2             
    reparamtimestep = 1e-2                  
    passswitchpointnsteps = 5                
    discrtimestep = 1e-2                    
    assert(dur > 10.0*discrtimestep)
    
    ncollision = 0
    nnotretimable = 0 
    nnotshorter = 0
    
    transtraj, rtraj = TransRotTrajFromSE3Traj(se3traj)
    lietraj = lie.SplitTraj(Rlist, rtraj)
   

    for it in range(maxiter):
        if plotdura == 1:
            plt.scatter(it, se3traj.duration)
            plt.draw()
        #transtraj, rtraj = TransRotTrajFromSE3Traj(se3traj)
        #lietraj = lie.SplitTraj2(Rlist, rtraj)
        if (expectedduration > 0): # check, if newlietraj.duration is short enough, stop SHORTCUTING
            if (se3traj.duration < expectedduration):
                print "\033[1;32mTrajectory's duration is already shorter than expected time --> stop shortcuting\033[0m"
                break
        if (dur < discrtimestep):
            print "[Utils::Shortcut] trajectory duration is less than discrtimestep.\n"
            break ## otherwise, this will cause an error in TOPP        
        
        ## select an interval for shortcutting
        t0 = random.rand()* dur
        
        if meanduration == 0:
            meanduration = dur - t0
            
        T = random.rand()*min(meanduration,dur - t0)
        t1 = t0 + T

        while (T < 2.0*discrtimestep):
            t0 = random.rand()*dur
            if meanduration == 0:
                meanduration = dur - t0
                
            T = random.rand()*min(meanduration, dur - t0)
            t1 = t0 + T

            if t1 > upperlimit:
                t1 = upperlimit
                if (t1 < t0):
                    temp = t0
                    t0 = t1
                    t1 = temp
                    T = t1 - t0

        # print "\n\nShortcutting iteration", it + 1
        # print t0, t1, t1- t0       
        # interpolate from t0 to t1
        R_beg = lietraj.EvalRotation(t0)
        R_end = lietraj.EvalRotation(t1)
        omega0 = lietraj.EvalOmega(t0)
        omega1 = lietraj.EvalOmega(t1)
        shortcutrtraj = lie.InterpolateSO3(R_beg,R_end,omega0,omega1, T)

        t_beg = transtraj.Eval(t0)
        t_end = transtraj.Eval(t1)
        v_beg = transtraj.Evald(t0)
        v_end = transtraj.Evald(t1)
        
        shortcuttranstraj = Trajectory.PiecewisePolynomialTrajectory.FromString(TrajString3rdDegree(t_beg,t_end,v_beg,v_end, T))
        
        shortcutse3traj = SE3TrajFromTransandSO3(shortcuttranstraj, shortcutrtraj)
        #check feasibility only for the new portion
        
        isincollision = CheckCollisionSE3Traj(robot, shortcuttranstraj, shortcutrtraj, R_beg, discrtimestep)
        if (not isincollision):
            a,b,c = ComputeSE3Constraints(shortcutse3traj, taumax, fmax, discrtimestep)
            topp_inst = TOPP.QuadraticConstraints(shortcutse3traj, discrtimestep, vmax, list(a), list(b), list(c))
            x = topp_inst.solver
            ret = x.RunComputeProfiles(1,1) 
            if (ret == 1):
                x.resduration
                ## check whether the new one has shorter duration
                if (x.resduration + 0.1 < T): #skip if not shorter than 0.1 s
                    
                    x.ReparameterizeTrajectory()
                    x.WriteResultTrajectory()
                    TOPPed_shortcutse3traj = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
                    TOPPed_shortcuttranstraj, TOPPed_shortcutrtraj = TransRotTrajFromSE3Traj(TOPPed_shortcutse3traj)
                    newlietraj = ReplaceTrajectorySegment(lietraj, TOPPed_shortcutrtraj , t0, t1)
                   
                    newtranstraj = ReplaceTransTrajectorySegment(transtraj, TOPPed_shortcuttranstraj, t0,t1)

                    #####################################################
                    newrtraj = Trajectory.PiecewisePolynomialTrajectory.FromString(TrajStringFromTrajList(newlietraj.trajlist))
                    newse3traj = SE3TrajFromTransandSO3(newtranstraj,newrtraj)

                    Rlist = newlietraj.Rlist
                    rtraj = newrtraj
                    transtraj = newtranstraj
                    lietraj = newlietraj
                    se3traj = newse3traj
                    
                    dur = se3traj.duration

                    #print "*******************************************"
                    print "Success at iteration",it + 1,":", t0, t1,"Deta_t:", t1 - t0 - x.resduration
                    attempt += 1
                    #print "T:", nnotretimable, "; S:", nnotshorter , "; C:", ncollision , "; OK:", attempt
                    #print "*******************************************"
                else:
                    # print "Not shorter"
                    nnotshorter += 1
            else: 
                # print "Not retimable"
                nnotretimable += 1
        else:
            # print "Collision"
            ncollision += 1

            # print "T:", nnotretimable, "; S:", nnotshorter , "; C:", ncoll