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
    Generate PLots for Veloctity and Speed discontinuities

    Examples:
        >> bsp2cosmic.py gen_LLO_to_NRHO_imp_ext7d_BSP.bsp -tl 1 -dt 10 -dv 20

    Note:

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
cosmicTemp  = cosmicTempPath + 'cosmicTemp_EM.py'


# ============================================================================
# Parse User Inputs:
# ============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("inputBSPfile", metavar="inputFile.bsp", help="BSP Input file")
parser.add_argument("-n", "--outputName",
                     default="", help = "Optional output file name. Default:... ")
parser.add_argument("-id","--spiceID", default= "-30100", help="Spice ID. defautl = '' (use non-Body id on bsp) ")
parser.add_argument("-to","--tiniOffset", default="0", help="t0 offset- days. Default='' (use t0 from bsp)")
parser.add_argument("-tl","--tlInterval", default="0", help="Trajectory time line interval from t0_offset. Default='' (use bsp time span)")
parser.add_argument("-dv","--minDvSrch", default="100",
                    help="Min velocity Discontinuity (in m/s) to account as a ∆V. Default: 100  \
                     (Warning, usually a lower values is better, e.g. 10 m/s, \
                     but some time it can produce  large number of ficticios maneuvers \
                     due to fast angular rotation near peripasis/flybys )")
parser.add_argument("-dt","--dtDvSrch", default="10", help="Velocity Discontinuity Search Step Size (in sec). Default: 10 (recommended 1)")


# TODO: Add S/c Name: e.g. 'encnf5'
scName = 'mySC'

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

# periEventFile = outputFolder + '/periEvents_out.json'
# apoEventFile = outputFolder + '/apoEvents_out.json'
dvDiscFile = outputFolder + '/dvDiscEvents_out.json'

# Insert sc_name and spiceID
scID = int(args.spiceID)

#SpiceName.bodyInsert(scID,'mySC')
SpiceName.bodyInsert(scID,scName)

#outputLevel = int(args.outputLevel)

tOffsetDays  = 2.0/86400   # Cut ini and end of BSP to avoid conflics reading data

# ============================================================================
# Params Setup
# ============================================================================

#------------------------------------------------------
# TODO: MOVE to Json:
#------------------------------------------------------
#DV Disc search Params:
SminDVSrch = float(args.minDvSrch)*m/sec
dvSrchtimeStep = float(args.dtDvSrch)*sec

dvSrchDt = 5*dvSrchtimeStep  #  width of DV pulse - expected < dvSearDT. Use + and - dvSrchDt
dvSrchCenter = 'Moon'
dvSrchFrame = 'EMO2000'
#dvSrchFrame = 'IAU Moon Fixed'
#------------------------------------------------------

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

#---------------- --------------------------------------
# Trajectory Interval set up:
#------------------------------------------------------

#if args.tiniOffset: t0_offset = args.tiniOffset
#if args.tlInterval: tl_interval = args.tlInterval

t0_offset =  float(args.tiniOffset)*day +  tOffsetDays*day  # Default = 0
tl_interval = float(args.tlInterval)*day  -2*tOffsetDays*day  # Default = 0: If == 0 -> use tl_end()

if M.TrajSetBoa.read(boa).intervals(scName):
    trajInterval = M.TrajSetBoa.read(boa).intervals(scName)
else:
    print(" *** ERROR: Incorrect scID ")
    exit()

traj_t0 = trajInterval[0].begin() + t0_offset
if not tl_interval:
    traj_tf = trajInterval[0].end()
else:
    traj_tf = traj_t0 +tl_interval

if tl_interval == 0:
    traj_tf = trajInterval[0].end()
searchInterval = TimeInterval(traj_t0,traj_tf)

print('... Time Interval for Scaning:')
print('    t0 = ' + str(traj_t0))
print('    tf = ' + str(traj_tf))
print('    dt = ' + str(traj_tf - traj_t0))

#------------------------------------------------------

# ============================================================================
# Short Functions:
# ============================================================================
# def geDVfromSMAs(sma1, sma2, radius, centralBody):
#         mu = GmBoa.read(boa,centralBody).gm()
#         dv_out = sqrt(mu*(2/radius - 1/sma1))-sqrt(mu*(2/radius - 1/sma2))
#         return dv_out

# ============================================================================
# SCAN trajectory:
# ============================================================================

#-----------------------------------------------------------------------------
#Extract central Bodies:
# TODO: Read naturalBodies on bsp, i.e. -> before loading defaultdata and ephem
naturalBodies =  TrajSetBoa.read(boa).getAll()
naturalBodies.remove(scName)

#-----------------------------------------------------------------------------
# Find DVs:
#-----------------------------------------------------------------------------

# ------------------------------------------------------
#Find DV Method 1:  Find Dv. disct, compute DV = Vi+1-Vi
speedDiscPLot = True
if speedDiscPLot:
    # print('... Using DV Disc. from: ' + dvDiscFile )
    # with open(dvDiscFile, 'r') as jsonInput:
    #    dvSearchDic = json.load(jsonInput)

    t1_cpu = process_time()
    print( '... Searching for DV discontinuities: ')
    dvSearchDic=[]
    querySat = TrajQuery( boa,scName,dvSrchCenter,dvSrchFrame)
    timeStep = dvSrchtimeStep
    tList =Epoch.range(traj_t0,traj_tf, timeStep)
    numDvDisc = 0
    dvTot = 0
    dvMag_list = []
    for tt in tList[:-1]:
        #dvMag = (querySat.state(tt+timeStep).vel() - querySat.state(tt).vel()).mag()/querySat.state(tt).vel().mag()
        dvMag = (querySat.state(tt+timeStep).vel() - querySat.state(tt).vel()).mag()
        dvMag_list.append(dvMag)
        if dvMag > SminDVSrch.value():
            dv_vec = querySat.state(tt+timeStep).vel() - querySat.state(tt).vel()
            dvSearchDic.append({
                'time' : str(tt),
                'center' : dvSrchCenter ,
                'frame' : dvSrchFrame,
                'eventType' : 'dvDisc',
                'value' : [dv_vec[0],
                           dv_vec[1],
                           dv_vec[2]],
                'dv_mag' : dvMag,
                'units' : 'km/sec',
                 })
            dvTot += dvMag
            numDvDisc += 1

    print('     Num of DV Disc. : ' + str(numDvDisc))
    print('     DV Total : ' + str(dvTot) + 'km/s')
    t2_cpu = process_time()
    print(f"     time: {(t2_cpu - t1_cpu)} sec")

    mpylab.ylabel( " Velocity Discontinuities(km/sec)" )
    mpylab.plot(tList[:-1],dvMag_list)
    mpylab.show()


# # ============================================================================
# # Create Cosmic Timeline:
# # ============================================================================
# # Create Manger:
# print('... Creating Manager.')
# mgr = Manager(boa)
# # quiet manager to avoid: WARNING COSMIC/Output:0 control points, when loadInput
# mgr.quiet = True
# mgr.loadInput(cosmicTemp)
# mgr.quiet = False
#
# # ============================================================================
# # Create traj Query:
# # ============================================================================
#
# stQuery= M.TrajQuery( boa, scName,dvSrchCenter,dvSrchFrame)
#
# # ----------------------------------------------------------------------------
# # ADD CP00 (at tl_0: mgr.tl.begin())
# # ----------------------------------------------------------------------------
# # #ADD CP at  bsp_traj_t0
# print('... Adding fix CP at traj_t0')
# cpDepCenter= 'Moon'
# cpDepFrame =  'EMO2000'
# cpName = "CP00"
# cpTime = M.Epoch(traj_t0)
# cpState = [*stQuery.state(cpTime).pos()*km ,
#            *stQuery.state(cpTime).vel()*km/sec]
# mcpUtil.appendCpCart(
#     mgr,
#     cpName ,
#     cpTime,
#     cpState,
#     center=cpDepCenter,
#     frame=cpDepFrame,
#     body=scName,
#     propagator="DIVA",
#     fixCP = True
# )
#
# # # ----------------------------------------------------------------------------
# # # ADD DV's CPs at DV_Disc (at LOIs)
# # # ----------------------------------------------------------------------------
# print('... Adding CPs and DVs at dvDisc.:')
# cpPeriFile = addDVs[0]['source']
# # Check if file exist:
# print('    Loading: ' + dvDiscFile )
# with open(dvDiscFile, 'r') as jsonInput:
#    dvDiscDic = json.load(jsonInput)
#
# for ii, dv in enumerate(dvDiscDic):
#     dvName ='DV'+str(ii+1).zfill(2)
#     cpName = 'CP-'+ dvName
#     cpTime = dv['time']
#     mcpUtil.appendCpCoe_V3(
#         mgr,
#         cpName,
#         str(M.Epoch(cpTime)+dvSrchDt),
#         stQuery,
#         mass=1000*kg,
#         center=dvSrchCenter,
#         frame=dvSrchFrame,
#         body=scName,
#         propagator="DIVA",
#         #controls =[]
#     )
#
#     mcpUtil.addCpBurn(
#         mgr, dvName,cpName,
#         #tDelta = -1*dvSrchDtPre,   # Use post instead  os Post to perfomr DV at DV_Disc (See CP is added at + dvSrchDt1)
#         tDelta = -1*dvSrchtimeStep,
#         frame=dv['frame'],
#         dvel=M.Dbl3Vec(dv['value']),
#         dvBound=2*dv['dv_mag']*km/s,
#     )
# print('    -> '+str(len(dvDiscDic))+' DVs Added!')
# #
# # ----------------------------------------------------------------------------
# # ADD CP at some tf: (fixed)
# # ----------------------------------------------------------------------------
# # #ADD CP at  bsp_traj_tf
# print('... Adding fix CP at traj_tf')
# cpName = "CP-END"
# cpTime = M.Epoch(traj_tf)
# cpState = [*stQuery.state(cpTime).pos()*km ,
#            *stQuery.state(cpTime).vel()*km/sec]
# mcpUtil.appendCpCart(
#     mgr,
#     cpName ,
#     cpTime,
#     cpState,
#     center=dvSrchCenter,
#     frame=dvSrchFrame,
#     body=scName,
#     propagator="DIVA",
#     fixCP = True
# )
#
#
# # ============================================================================
# # Save Files and Data:
# # ============================================================================
#
# # Propagate  and Save Trajectory:
# mgr.tl.createTraj(boa, mgr.problem)
# mgr.saveChkPt(outputName, allowOverwrite = True)


# ============================================================================
# ============================================================================
#  Just scan BSP  file:

# -> find Natural bodies
# -> Find Peri/Apos
# -> Find DV discontinuities
#
# ============================================================================
# ============================================================================

# =======================================================================
# END
# =======================================================================
