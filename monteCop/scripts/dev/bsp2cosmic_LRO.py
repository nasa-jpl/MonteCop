#!/usr/bin/env mpython_q

# ===========================================================================
# Section 392 Navigation and Mission Design
#
# Copyright (C) 2021, California Institute of Technology.
# U.S. Government Sponsorship under NASA Contract NAS7-03001 is acknowledged.
#
# ===========================================================================

# ============================================================================
# Imports here:
# ============================================================================
'''
    Generate a Cosmic tl from bsp for LRO mission.

    Options:  -to 1 -> to skip Earth departure and quick velocity change
              -tl  10 -> to reduce timeline to 10 searchdays
              -o  3 -> to plot dv disc. history (use -o 2  for vervose)

    RUN:  >> bsp2cosmic_LRO.py lro_tli2loi.bsp -to 1 -tl 10 -o 3 -ov

'''


import Monte as M
from mpy.units import *
import mpy.io.data as defaultData
from mpy.opt.cosmic import Manager
from mpy.units import s, hour, m, minute, km, day, deg, rad, kg

import matplotlib.pyplot as plt

import ntpath
import argparse
import json
import os

from time import process_time

import monteCop.utils.cosmicUtils as  mcpUtil
import mpylab

# ============================================================================

# Default Data: (OVERWIRTE BY JSON CONFIG)
boaPlanets = '/nav/common/import/ephem/de430.boa'
cosmicTempPath  = '/home/ricgomez/lib/monteCop/templates/'
#cosmicTemp  = cosmicTempPath + 'lateEncnf5.py'
cosmicTemp  = cosmicTempPath + 'cosmicTemp_EM.py'


# ============================================================================
# Parse User Inputs:
# ============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("inputBSPfile", metavar="inputFile.bsp", help="BSP Input file")
parser.add_argument("-n", "--outputName",
                     default="", help = "Optional output file name. Default:... ")
parser.add_argument("-id","--spiceID", default= "-85", help="Spice ID. defautl = '' (use non-Body id on bsp) ")
parser.add_argument("-to","--tiniOffset", default="0", help="t0 offset- days. Default='' (use t0 from bsp)")
parser.add_argument("-tl","--tlInterval", default="0", help="Trajectory time line interval from t0_offset. Default='' (use bsp time span)")
parser.add_argument('-o', "--outputLevel", default = 2,
                    help='outputLevel: outputLevel = 1 -> msgs; outputLevel = 2 -> save JSON files) outputLevel = 3 -> Plots;')
parser.add_argument('-ov', action='store_true', help='Files Overwrite -> Overwrite data.json files')



# TODO: Add S/c Name: e.g. 'encnf5'
scName = 'lro'

args = parser.parse_args()

# ============================================================================
# PARSE INPUTS ::

trajBSP = args.inputBSPfile
baseName = trajBSP.replace('.bsp','')

if args.outputName in '':
    outputName = ntpath.basename(trajBSP).replace('.bsp','_b2m.py')
else:
    outputName = args.outputName
outputFolder = baseName + '_TMP'

periEventFile = outputFolder + '/periEvents_out.json'
periNoFBsEventFile = outputFolder + '/periNoFBsEvents_out.json'
apoEventFile = outputFolder + '/apoEvents_out.json'
dvDiscFile = outputFolder + '/dvDiscEvents_out.json'

# Insert sc_name and spiceID
scID = int(args.spiceID)

#SpiceName.bodyInsert(scID,'mySC')
SpiceName.bodyInsert(scID,scName)

outputLevel = int(args.outputLevel)

if outputLevel >= 2:
    if not  os.path.exists(outputFolder):
        os.makedirs( outputFolder )

# ============================================================================
# Params Setup
# ============================================================================

#------------------------------------------------------
# MOVE to Json:
#------------------------------------------------------

#Apsis Search Params:
#apsisSearchStep = 0.1*day

#DV Disc search Params:
SminDVSrch = 18.0*m/sec  # Tunned for LRO
dvSrchtimeStep =  10.0*sec
dvSrchDt1 = 5*minute  #  width of DV pulse - expected < dvSearDT. Use + and - dvSrchDt
dvSrchDt2 = 40*minute  # Use for exceptions: i.w. Finite Burn LOI DV.
dvSrchCenter = 'Moon'
#dvSrchFrame = 'EMO2000'
dvSrchFrame = 'IAU Moon Fixed'


dvSearch = [
    {'centerBody' : dvSrchCenter,
     'minDV' : SminDVSrch,
     'timeStep' :  dvSrchtimeStep ,
     'searchType' : 'dvDisc',
     'searchFrame' : 'EMO2000',
    }
]

addDVs = [
    {
    'source' : dvDiscFile,
    'name' : 'dvDisc',
    'defualts' :
        {'centerBody' : dvSrchCenter,
         'Frame' : 'EMO2000',
         'params' : 'params',
        },
    }
]

#------------------------------------------------------

# ============================================================================
# Initializations:
# ============================================================================

#-----------------------------------------------------------------------------
# -> BOA, create and load data/ephems :
#-----------------------------------------------------------------------------

# Load boa
boa = M.BoaLoad()
boa.load( boaPlanets )
#boa.load( boaSats )
boa.load( trajBSP )

# Load defautls
defaultData.loadInto(boa,["frame","body","frame/IAU 2000","frame/inertial"])

#------------------------------------------------------
# Trajectory Interval set up:
#------------------------------------------------------

if args.tiniOffset: t0_offset = args.tiniOffset
if args.tlInterval: tl_interval = args.tlInterval

t0_offset =  float(args.tiniOffset)*day  # Default = 0
tl_interval = float(args.tlInterval)*day  # Default = 0: If == 0 -> use tl_end()

trajInterval = M.TrajSetBoa.read(boa).intervals(scName)
traj_t0 = trajInterval[0].begin() + t0_offset
if not tl_interval:
    traj_tf = trajInterval[0].end()
else:
    traj_tf = traj_t0 +tl_interval

if tl_interval == 0:
    traj_tf = trajInterval[0].end()
searchInterval = TimeInterval(traj_t0,traj_tf)

if outputLevel >= 2:
    print('... Time Interval for Scaning:')
    print('    t0 = ' + str(traj_t0))
    print('    tf = ' + str(traj_tf))
    print('    dt = ' + str(traj_tf - traj_t0))
#------------------------------------------------------

# ============================================================================
# Short Functions:
# ============================================================================
def geDVfromSMAs(sma1, sma2, radius, centralBody):
        mu = GmBoa.read(boa,centralBody).gm()
        dv_out = sqrt(mu*(2/radius - 1/sma1))-sqrt(mu*(2/radius - 1/sma2))
        return dv_out

# ============================================================================
# SCAN trajectory:
# ============================================================================

#-----------------------------------------------------------------------------
#Extract central Bodies:
naturalBodies =  TrajSetBoa.read(boa).getAll()
naturalBodies.remove(scName)

#-----------------------------------------------------------------------------
# Find DVs:
#-----------------------------------------------------------------------------

# # ------------------------------------------------------
# #Find DV Method 1:  Find Dv. disct, compute DV = Vi+1-Vi
# #if findDvDiscon and not os.path.exists(dvDiscFile):
# findDvDiscon = False
# if findDvDiscon and os.path.exists(dvDiscFile) and not args.ov:
#     print('... Using DV Disc. from: ' + dvDiscFile )
#     with open(dvDiscFile, 'r') as jsonInput:
#        dvSearchDic = json.load(jsonInput)
#
# if findDvDiscon and (args.ov or not os.path.exists(dvDiscFile)):
#     t1_cpu = process_time()
#     print( '... Searching for DV discontinuities: ')
#     dvSearchDic=[]
#     querySat = TrajQuery( boa,scName,dvSrchCenter,dvSrchFrame)
#     timeStep = dvSrchtimeStep
#     tList =Epoch.range(traj_t0,traj_tf, timeStep)
#     numDvDisc = 0
#     dvTot = 0
#     dvMag_list = []
#     for tt in tList[:-1]:
#         #dvMag = (querySat.state(tt+timeStep).vel() - querySat.state(tt).vel()).mag()/querySat.state(tt).vel().mag()
#         dvMag = (querySat.state(tt+timeStep).vel() - querySat.state(tt).vel()).mag()
#         dvMag_list.append(dvMag)
#         if dvMag > SminDVSrch.value():
#             dv_vec = querySat.state(tt+timeStep).vel() - querySat.state(tt).vel()
#             dvSearchDic.append({
#                 'time' : str(tt),
#                 'center' : dvSrchCenter ,
#                 'frame' : dvSrchFrame,
#                 'eventType' : 'dvDisc',
#                 'value' : [dv_vec[0],
#                            dv_vec[1],
#                            dv_vec[2]],
#                 'dv_mag' : dvMag,
#                 'units' : 'km/sec',
#                  })
#             dvTot += dvMag
#             numDvDisc += 1
#
#     print('     Num of DV Disc. : ' + str(numDvDisc))
#     print('     DV Total : ' + str(dvTot) + 'km/s')
#     t2_cpu = process_time()
#     print(f"     time: {(t2_cpu - t1_cpu)} sec")

# --------------------------------------------------------------------
#Find DV Method 2:  Find Dv. disct, compute ∆sma = sma(t+dt) -sma(t-dt)
#                   Compute |DV|= DV(∆sma). -> Add along-track velocity.
#                   From the vis-viva equation: |∆v| = () -()
#                   TODO: Start drwing diagrams, and equations for paper
#                         This is a fun problem.
# --------------------------------------------------------------------
findDvDiscon = True
#findDvDiscon = False
if findDvDiscon and os.path.exists(dvDiscFile) and not args.ov:
    print('... Using DV Disc. from: ' + dvDiscFile )
    with open(dvDiscFile, 'r') as jsonInput:
       dvSearchDic = json.load(jsonInput)

if findDvDiscon and (args.ov or not os.path.exists(dvDiscFile)):
    t1_cpu = process_time()
    print( '... Searching for DV discontinuities: ')
    dvSearchDic=[]
    velFrame = 'VUW_'+dvSrchCenter
    M.BodyVelDirFrame( boa, velFrame ,'EMO2000',TimeInterval(),scName,dvSrchCenter)
    querySat = TrajQuery( boa,scName,dvSrchCenter,dvSrchFrame)
    timeStep = dvSrchtimeStep
    tList =Epoch.range(traj_t0,traj_tf, timeStep)
    numDvDisc = 0
    dvTot = 0
    dvMag_list = []
    dvSrchDtPost = dvSrchDt1
    skipBellowThisTime =  Epoch.distantPast()
    for tt in tList[:-1]:
        #NOTE: dvSrchDt : width of DV pulse - expected < dvSearDT. Use + and - dvSrchDt
        #      When pulse/dv detected,SKIP skip search till tt > skipTimeMax (till tt > DV_time + dvSrchDt)
        dvMag = (querySat.state(tt+timeStep).vel() - querySat.state(tt).vel()).mag()
        dvMag_list.append(dvMag)
        if tt < skipBellowThisTime:
            continue
        if dvMag > SminDVSrch.value():
            skipBellowThisTime = tt + dvSrchDt1
            # to fix finite burn at LOI
            dvSrchDtPre = dvSrchDt2 if numDvDisc == 0 else dvSrchDtPre
            #  TODO: Discrivbe this logic here!!  -45min  AND +5min for LOI1, -5 to 5 min for LOI3...5
            sma1 = M.Conic.semiMajorAxis(querySat.state(tt-dvSrchDtPre))
            sma2 = M.Conic.semiMajorAxis(querySat.state(tt+dvSrchDtPost))
            range = M.Conic.radius(querySat.state(tt))
            dvMag = geDVfromSMAs(sma2, sma1, range, dvSrchCenter)
            print(f"dv={dvMag} epoch:{tt}")
            dv_vec =  [0, 1*dvMag.value(),  0 ]
            dvSearchDic.append({
                'time' : str(tt),
                'center' : dvSrchCenter ,
                'frame' : velFrame,
                'eventType' : 'dvDisc',
                'value' : [dv_vec[0],
                           dv_vec[1],
                           dv_vec[2]],
                'dv_mag' : abs(dvMag.value()),
                'units' : 'km/sec',
                 })
            dvTot += dvMag
            numDvDisc += 1

    print('     Num of DV Disc. : ' + str(numDvDisc))
    print('     DV Total : ' + str(dvTot) + 'km/s')
    t2_cpu = process_time()
    print(f"     time: {(t2_cpu - t1_cpu)} sec")

# -----------------
#Print or Save Data:
findDvDiscon = True
if findDvDiscon and (not os.path.exists(dvDiscFile) or args.ov):
    # Data Save:
    saveDvData = True
    if saveDvData:
        with open(dvDiscFile, 'w+' ) as outfile:
           json.dump( dvSearchDic, outfile, indent = 4, separators=(',', ': ') )
        print('     ' + ntpath.basename(dvDiscFile) + ' Saved!' )

    if outputLevel >= 3:
        mpylab.plot(tList[:-1],dvMag_list)
        mpylab.show()

# ============================================================================
# Create Cosmic Timeline:
# ============================================================================
# Create Manger:
print('... Creating Manager.')
mgr = Manager(boa)
# quiet manager to avoid: WARNING COSMIC/Output:0 control points, when loadInput
mgr.quiet = True
mgr.loadInput(cosmicTemp)
mgr.quiet = False

# ----------------------------------------------------------------------------
# ADD CP and DV: At tl_0: mgr.tl.begin()
# ----------------------------------------------------------------------------
# #ADD CP at  bsp_traj_t0
print('... Adding fix CP at tl_t0: mgr.tl.begin().')
cpDepCenter= 'Earth'
cpDepFrame =  'EMO2000'
stDepQuery= M.TrajQuery( boa, scName,cpDepCenter,cpDepFrame)
cpName = "CP00"
cpTime = M.TrajSetBoa.read(boa).intervals(scName)[0].begin()
cpState = [*stDepQuery.state(cpTime).pos()*km ,
           *stDepQuery.state(cpTime).vel()*km/sec]
mcpUtil.appendCpCart(
    mgr,
    cpName ,
    cpTime,
    cpState,
    center=cpDepCenter,
    frame=cpDepFrame,
    body=scName,
    propagator="DIVA",
    fixCP = True
)

# ADD DV_null At t0 and tf of Timeline:
print('    Adding DV at tl_t0: mgr.tl.begin()')
dvName = 'DV00'
mcpUtil.addCpBurn(
    mgr, dvName, cpName,
    tDelta = 1.0*sec,
    frame=cpDepFrame,
    dvel=M.Dbl3Vec([1e-5] * 3),
    dvBound=0.1*km/s,
    )

# ----------------------------------------------------------------------------
# ADD DV's CPs at DV_Disc (at LOIs)
# ----------------------------------------------------------------------------
print('... Adding CPs and DVs at dvDisc.:')
cpPeriFile = addDVs[0]['source']
# Check if file exist:
print('    Loading: ' + dvDiscFile )
with open(dvDiscFile, 'r') as jsonInput:
   dvDiscDic = json.load(jsonInput)

stQuery= M.TrajQuery( boa, scName,dvSrchCenter,dvSrchFrame)

for ii, dv in enumerate(dvDiscDic):
    dvName ='DV'+str(ii+1).zfill(2)
    cpName = 'CP-'+ dvName
    cpTime = dv['time']
    dvSrchDtPre = dvSrchDt2 if ii == 0 else dvSrchDt1
    mcpUtil.appendCpCoe_V3(
        mgr,
        cpName,
        str(M.Epoch(cpTime)+dvSrchDtPost),
        stQuery,
        mass=1000*kg,
        center=dvSrchCenter,
        frame=dvSrchFrame,
        body=scName,
        propagator="DIVA",
        #controls =[]
    )
    # if ii == 0
    # ... Add Cart...

    mcpUtil.addCpBurn(
        mgr, dvName,cpName,
        tDelta = -1*dvSrchDtPre,   # Use post instead  os Post to perfomr DV at DV_Disc (See CP is added at + dvSrchDt1)
        frame=dv['frame'],
        dvel=M.Dbl3Vec(dv['value']),
        dvBound=2*dv['dv_mag']*km/s,
    )
print('    -> '+str(len(dvDiscDic))+' DVs Added!')

# ----------------------------------------------------------------------------
# ADD CP at some tf: (fixed)
# ----------------------------------------------------------------------------
# #ADD CP at  bsp_traj_tf
print('... Adding fix CP at tl_t0: mgr.tl.begin().')
cpName = "CP-END"
cpTime = M.Epoch(traj_tf)
cpState = [*stQuery.state(cpTime).pos()*km ,
           *stQuery.state(cpTime).vel()*km/sec]
mcpUtil.appendCpCoe_V3(
    mgr,
    cpName,
    cpTime,
    stQuery,
    mass=1000*kg,
    center=dvSrchCenter,
    frame=dvSrchFrame,
    body=scName,
    propagator="DIVA",
    #controls =[]
)


#TODO:  For LOI-DV, use dvSearDT = 40(-) and 5(+) minutes, for others 5(+/-) min.
#       AND: use tt-DT to apply maneuver, otherwise pre-mvr state is going to be already an updated state
#       -commit
#       -Opt only  capture

# ============================================================================
# Save Files and Data:
# ============================================================================

# Propagate  and Save Trajectory:
mgr.tl.createTraj(boa, mgr.problem)
mgr.saveChkPt(outputName, allowOverwrite = True)


# ============================================================================
# ============================================================================
#  Just scan BSP  file:

# -> find Natural bodies
# -> Find Peri/Apos
# -> Find DV discontinuities
#
# ============================================================================
# ============================================================================
