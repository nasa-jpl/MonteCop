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

import Monte as M
from mpy.units import *
import mpy.io.data as defaultData
from mpy.opt.cosmic import Manager
from mpy.units import s, hour, m, km, day, deg, rad, kg

import matplotlib.pyplot as plt

import ntpath
import argparse
import json
import os

from time import process_time

# ============================================================================


#-----------------------------------------------------------------------------
# Parse User Inputs:
#-----------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("inputBSPfile", metavar="inputFile.bsp", help="BSP Input file")
parser.add_argument("-n", "--outputName",
                     default="", help = "Optional output file name. Default:... ")
parser.add_argument("-id","--spiceID", default= "-303", help="Spice ID. defautl = '' (use non-Body id on bsp) ")
parser.add_argument("-to","--tiniOffset", default="0", help="t0 offset- days. Default='' (use t0 from bsp)")
parser.add_argument("-tl","--tlInterval", default="0", help="Trajectory time line interval from t0_offset. Default='' (use bsp time span)")
parser.add_argument('-o', "--outputLevel", default = 2,
                    help='outputLevel: outputLevel = 1 -> msgs; outputLevel = 2 -> save JSON files) outputLevel = 3 -> Plots;')

# parser.add_argument('-b','--bodyList', nargs='+', help='<Required> Set flag')
# parser.add_argument('-c', '--bodyCenter',help = 'Body Center. \
#                     The defualt option "Natural body" will use the SPK natural body for the state at a given time')
# parser.add_argument('-f','--frame', help='Visualization frame (default: "J2000")')
# parser.add_argument("jsonConfig", default = None, nargs='?',
#                     help='json config file with predefined setup information')
# parser.add_argument('-s', action='store_true', help='Save json file solution')

args = parser.parse_args()

#-----------------------------------------------------------------------------
# Params Setup
#-----------------------------------------------------------------------------

if args.tiniOffset: t0_offset = args.tiniOffset
if args.tlInterval: tl_interval = args.tlInterval

# t0_offset =  80.0*day  # Default = 0
# tl_interval = 5.0*day  # Default = 0: If == 0 -> use tl_end()
t0_offset =  float(args.tiniOffset)*day  # Default = 0
tl_interval = float(args.tlInterval)*day  # Default = 0: If == 0 -> use tl_end()


#Apsis Search Params:
apsisSearchStep = 0.1*day

#DV Disc search Params:
minDVSearch = 1.0*m/sec
dvSearchtimeStep = 100.0*sec
dvSearchCenter = 'Saturn'
dvSearchFrame = 'EMO2000'

# Data Save:
saveDvData = False

# ----------------------------------------------------------------------------
# -> USER UNPUTS ::
trajBSP = args.inputBSPfile
baseName = trajBSP.replace('.bsp','')

if args.outputName in '':
    outputName = ntpath.basename(trajBSP).replace('.bsp','_out.py')
else:
    outputName = args.outputName
outputFolder = baseName + '_TMP'

periEventFile = outputFolder + '/periEvents_out.json'
periNoFBsEventFile = outputFolder + '/periNoFBsEvents_out.json'
apoEventFile = outputFolder + '/apoEvents_out.json'
dvDiscFile = outputFolder + '/dvDiscEvents_out.json'

# Insert sc_name and spiceID
scID = int(args.spiceID)
SpiceName.bodyInsert(scID,'mySC')

outputLevel = int(args.outputLevel)

# -> Default data:
cosmicTemp  = 'inputs/cosmicTemp.py'    # IF NOT in local folder, use default (lib/Templates)
boaPlanets = '/nav/common/import/ephem/de430.boa'
boaSats    = '/nav/common/import/ephem/sat375l.boa'

# Load boa
boa = M.BoaLoad()
#boa.load( boaPlanets )
#boa.load( boaSats )
boa.load( trajBSP )


#raise Exception('exit')

#-----------------------------------------------------------------------------
# Initializations:
#-----------------------------------------------------------------------------

#------------------------------------------------------
# Trajectory Scan time interval:
trajInterval = M.TrajSetBoa.read(boa).intervals('mySC')
traj_t0 = trajInterval[0].begin() + t0_offset
if not tl_interval:
    traj_tf = trajInterval[0].end()
else:
    traj_tf = traj_t0 +tl_interval

if tl_interval == 0:
    traj_tf = trajInterval[0].end()
searchInterval = TimeInterval(traj_t0,traj_tf)

if outputLevel >= 2:
    print(' ... Time Interval for Scaning:')
    print('     t0 = ' + str(traj_t0))
    print('     tf = ' + str(traj_tf))
    print('     dt = ' + str(traj_tf - traj_t0))

# ------------------------------------------------------

#raise Exception('exit')

#Extract central Bodies:
naturalBodies =  TrajSetBoa.read(boa).getAll()
naturalBodies.remove('mySC')

#Load additinoal ephemeris kernels
boa.load( boaPlanets )
boa.load( boaSats )

# Load defautls
defaultData.loadInto(boa,["frame","body"])

if outputLevel >= 2:
    if not  os.path.exists(outputFolder):
        os.makedirs( outputFolder )

#Bodies to read Close Approach: Need to come from user inputFile/user inputs
searchBodies = [
    {'bodyName' : 'Enceladus',
     'minAlt' : 10000.0*km,
     'searchType' : 'Peri',
     'searchFrame' : 'EMO2000',
     },
     {'bodyName' : 'Saturn',
      'minAlt' : 30000.0*km, #(??)
      'searchType' : 'periNoFB',
      'searchFrame' : 'EMO2000',
      },
    {'bodyName' : 'Saturn',
     'minAlt' : 10000.0*km, #(??)
     'searchType' : 'Apo',
     'searchFrame' : 'EMO2000',
     }
]

#-----------------------------------------------------------------------------
# Seach for Peris at Enceladus:
findPerisEnc = True
if findPerisEnc:
    # TODO: Replace by search intervals below minAlt. then search for Peris.
    # Strategy, search per searchBodies, and do the follow:
    ev = ApsisEvent( TrajQuery( boa, "mySC", "Enceladus" ), ApsisEvent.PERIAPSIS )
    r = ev.search(searchInterval, apsisSearchStep )
    #Filter Periapsis by Min Alt:
    periEvents = [xx  for xx in r if xx.value() < searchBodies[0]['minAlt']]
    #Save to Json
    # Epoch, rpRange, Center, Frame
    periEventDic=[]
    numPeris = 0
    for pp in periEvents:
        periEventDic.append({
            'time' : str(pp.time()),
            'value': str(pp.value()),
            'center' : searchBodies[0]['bodyName'],
            'frame' : searchBodies[0]['searchFrame'],
            'eventType' : searchBodies[0]['searchType'],
             })
        numPeris += 1

    if outputLevel >= 2:
        #TODO: add periEventDic Info!
        # ONLY IF Not existing file
        with open(periEventFile, 'w+' ) as outfile:
           json.dump( periEventDic, outfile, indent = 4, separators=(',', ': ') )
        print(' ... ' + ntpath.basename(periEventFile) + ' Saved!' )
        print('     Peri Events Found: ' + str(numPeris))
    #encBFframe='IAU Enceladus Fixed'

#-----------------------------------------------------------------------------
# Seach for Peris of Sturn with not Enceladus Flyby:
findPeriNoFBs = False
if findPeriNoFBs:
    # TODO: Replace by search intervals below minAlt. then search for Peris.
    # Strategy, search per searchBodies, and do the follow:
    ev = ApsisEvent( TrajQuery( boa, "mySC", "Saturn" ), ApsisEvent.PERIAPSIS )
    r = ev.search( searchInterval, apsisSearchStep )
    encQuery = TrajQuery( boa, "mySC", "Enceladus" )
    #check that peri is far from Enceladus (above minAlt)
    #encAlt = encQuery.state(x.time(),'EMO2000',3).posMag()
    #Filter Periapsis by Min Alt:
    periNoFBsEvents = [xx  for xx in r if encQuery.state(xx.time(),searchBodies[1]['searchFrame'],3).posMag() > searchBodies[1]['minAlt']]
    periNoFBsEventDic=[]
    numPerisNoFBs = 0
    for pp in periNoFBsEvents:
        periNoFBsEventDic.append({
            'time' : str(pp.time()),
            'value': str(pp.value()),
            'center' : searchBodies[1]['bodyName'],
            'frame' : searchBodies[1]['searchFrame'],
            'eventType' : searchBodies[1]['searchType'],
             })
        numPerisNoFBs += 1

    if outputLevel >= 2:
        #TODO: add periNoFBsEventDic Info!
        # ONLY IF Not existing file
        with open(periNoFBsEventFile, 'w+' ) as outfile:
           json.dump( periNoFBsEventDic, outfile, indent = 4, separators=(',', ': ') )
        print(' ... ' + ntpath.basename(periNoFBsEventFile) + ' Saved!' )
        print('     Peri w No FBs Events Found: ' + str(numPerisNoFBs))
    #encBFframe='IAU Enceladus Fixed'


#-----------------------------------------------------------------------------
# Seach for Apoapsis of 'searchBodies':
findApos = True
if findApos:
    ev = ApsisEvent( TrajQuery( boa, "mySC", "Enceladus" ), ApsisEvent.APOAPSIS)
    r = ev.search( searchInterval, apsisSearchStep )
    #Filter Periapsis by Min Alt:
    #periEvents = [xx  for xx in r if xx.value() < searchBodies[1]['minAlt']]
    apoEvents = r

    #Save to Json
    # Epoch, rpRange, Center, Frame
    apoEventDic=[]
    numApos = 0
    for pp in apoEvents:
        apoEventDic.append({
            'time' : str(pp.time()),
            'value': str(pp.value()),
            'center' : searchBodies[2]['bodyName'],
            'frame' : searchBodies[2]['searchFrame'],
            'eventType' : searchBodies[2]['searchType'],
             })
        numApos += 1

    if outputLevel >= 2:
        #TODO: add periEventDic Info!
        # ONLY IF Not existing file
        with open(apoEventFile, 'w+' ) as outfile:
           json.dump( apoEventDic, outfile, indent = 4, separators=(',', ': ') )
        print(' ... ' + ntpath.basename(apoEventFile) + ' Saved!' )
        print('     Apo Events Found: ' + str(numApos))



# ----------------------------------------------------------------------------
#Find DV discontinuities:
#findDvDiscon = True
findDvDiscon = False
if findDvDiscon and os.path.exists(dvDiscFile):
    print(' ... Using DV Disc. from: ' + dvDiscFile )
    with open(dvDiscFile, 'r') as jsonInput:
       dvSearchDic = json.load(jsonInput)

if findDvDiscon and not os.path.exists(dvDiscFile):
    t1_cpu = process_time()
    print( ' ... Searching for DV discontinuities: ')
    dvSearchDic=[]
    querySat = TrajQuery( boa,'mySC',dvSearchCenter,dvSearchFrame)
    timeStep = dvSearchtimeStep
    tList =Epoch.range(traj_t0,traj_tf, timeStep)
    numDvDisc = 0
    dvTot = 0
    dvMag_list = []
    for tt in tList[:-1]:
        #dvMag = (querySat.state(tt+timeStep).vel() - querySat.state(tt).vel()).mag()/querySat.state(tt).vel().mag()
        dvMag = (querySat.state(tt+timeStep).vel() - querySat.state(tt).vel()).mag()
        dvMag_list.append(dvMag)
        if dvMag > minDVSearch.value():
            dv_vec = querySat.state(tt+timeStep).vel() - querySat.state(tt).vel()
            dvSearchDic.append({
                'time' : str(tt),
                'center' : dvSearchCenter ,
                'frame' : dvSearchFrame,
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

    if outputLevel >= 3:
        xx = range(0,len(dvMag_list))
        xx = [rr*(timeStep/86400/sec)  for rr in xx]
        plt.plot(xx,dvMag_list)
        plt.show()

    if saveDvData:
        with open(dvDiscFile, 'w+' ) as outfile:
           json.dump( dvSearchDic, outfile, indent = 4, separators=(',', ': ') )
        print('     ' + ntpath.basename(dvDiscFile) + ' Saved!' )

# ----------------------------------------------------------------------------

# NEXT:

# create a bsp2cosmic_EncTour.py:
#
# - Read json, create control points at Peris
# - Add velocities
# - [OPT] add CU mvrs ~ 5 hours after C/A


#-----------------------------------------------------------------------------
# Load bsp
#-----------------------------------------------------------------------------


# ============================================================================
#  Just scan BSP  file:

# -> find Natural bodies
# -> Find Peri/Apos
# -> Find DV discontinuities
#
# Then:
#
# -> add CPs at Peris (at <  minAlt)   (Use natural center or user defined center)
# -> Add CUs and APO Mvrs  (at, or t_event, or at t_delta from event )
# -> Only add constraints on minAlt, to and tf times
#
# Save data (Apo/Peri, DVs) to JASON, outputName_TMP/.

# ============================================================================
