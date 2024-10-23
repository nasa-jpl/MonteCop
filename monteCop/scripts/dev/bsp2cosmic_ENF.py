#!/usr/bin/env mpython_q

# ===========================================================================
# Section 392 Navigation and Mission Design
#
# Copyright (C) 2021, California Institute of Technology.
# U.S. Government Sponsorship under NASA Contract NAS7-03001 is acknowledged.
#
# ===========================================================================

# Examples:
#    >> bsp2cosmic_ENF.py Enceladus_Tour0.bsp -id -303 -to 1 -tl 30 -dt 10
#

# ============================================================================
# Imports here:
# ============================================================================

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

import monteCop.utils.cosmicUtils as  mcpUtil

# ============================================================================


# Default Data: (OVERWIRTE BY JSON CONFIG)
cosmicTemp  = 'inputs/cosmicTemp.py'    # IF NOT in local folder, use default (lib/Templates)
boaPlanets = '/nav/common/import/ephem/de430.boa'
boaSats    = '/nav/common/import/ephem/sat375l.boa'

cosmicTemp  = '/home/ricgomez/lib/monteCop/templates/cosmicTemplateEncnf5.py'

# ============================================================================
# Parse User Inputs:
# ============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("inputBSPfile", metavar="inputFile.bsp", help="BSP Input file")
parser.add_argument("-n", "--outputName",
                     default="", help = "Optional output file name. Default:... ")
parser.add_argument("-id","--spiceID", default= "-303", help="Spice ID. defautl = '' (use non-Body id on bsp) ")
parser.add_argument("-to","--tiniOffset", default="0", help="t0 offset- days. Default='' (use t0 from bsp)")
parser.add_argument("-tl","--tlInterval", default="0", help="Trajectory time line interval from t0_offset. Default='' (use bsp time span)")
parser.add_argument('-o', "--outputLevel", default = 3,
                    help='outputLevel: outputLevel = 1 -> msgs; outputLevel = 2 -> save JSON files) outputLevel = 3 -> Plots; Defualt:3')
parser.add_argument("-dv","--minDvSrch", default="100",
                    help="Min velocity Discontinuity (in m/s) to account as a ∆V. Default: 100  \
                     (Warning, usually a lower values is better, e.g. 10 m/s, \
                     but some time it can produce  large number of ficticios maneuvers \
                     due to fast angular rotation near peripasis/flybys )")
parser.add_argument("-dt","--dtDvSrch", default="10", help="Velocity Discontinuity Search Step Size (in sec). Default: 10 (recommended 1)")



# TODO: Add S/c Name: e.g. 'encnf5'
scName = 'encnf5'

# parser.add_argument('-b','--bodyList', nargs='+', help='<Required> Set flag')
# parser.add_argument('-c', '--bodyCenter',help = 'Body Center. \
#                     The defualt option "Natural body" will use the SPK natural body for the state at a given time')
# parser.add_argument('-f','--frame', help='Visualization frame (default: "J2000")')
# parser.add_argument("jsonConfig", default = None, nargs='?',
#                     help='json config file with predefined setup information')
# parser.add_argument('-s', action='store_true', help='Save json file solution')

args = parser.parse_args()

# ============================================================================
# PARSE INPUTS ::

trajBSP = args.inputBSPfile
baseName = trajBSP.replace('.bsp','')

if args.outputName in '':
    outputName = ntpath.basename(trajBSP).replace('.bsp','_B2M.py')
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
apsisSearchStep = 0.1*day

#DV Disc search Params:
minDVSearch = float(args.minDvSrch)*m/sec
dvSearchtimeStep = float(args.dtDvSrch)*sec
dvSearchCenter = 'Saturn'
dvSearchFrame = 'EMO2000'

dvSearch = [
    {'centerBody' : 'Saturn',
     'minDV' : 0.1*km/sec,
     'timeStep' :  100.0*sec,
     'searchType' : 'dvDisc',
     'searchFrame' : 'EMO2000',
    }
]

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

# CPs config
addCPs = [
    {
    'source' : periEventFile,
    'name' : 'cpEncPeri',
    'defualts' :
        {'centerBody' : 'Enceladus',
         'Frame' : 'EMO2000',
         'stateParams' : 'flybyState',
        },
    }
]

addDVs = [
    {
    'source' : dvDiscFile,
    'name' : 'dvDisc',
    'defualts' :
        {'centerBody' : 'Saturn',
         'Frame' : 'EMO2000',
         'params' : 'params',
        },
    }
]

addApoDVs = [
    {
    'source' : apoEventFile,
    'name' : 'apoDVs',
    'defualts' :
        {'centerBody' : 'Saturn',
         'Frame' : 'EMO2000',
         'params' : 'params',
        },
    }
]

addPeriDVs = [
    {
    'source' : periNoFBsEventFile,
    'name' : 'periDVs',
    'defualts' :
        {'centerBody' : 'Saturn',
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
boa.load( boaSats )
boa.load( trajBSP )

# Load defautls
defaultData.loadInto(boa,["frame","body"])

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
# SCAN trajectory:
# ============================================================================

#-----------------------------------------------------------------------------
#Extract central Bodies:
naturalBodies =  TrajSetBoa.read(boa).getAll()
naturalBodies.remove(scName)

#-----------------------------------------------------------------------------
# Seach for Peris at Enceladus:
findPerisEnc = True
if findPerisEnc:
    # TODO: Replace by search intervals below minAlt. then search for Peris.
    # Strategy, search per searchBodies, and do the follow:
    ev = ApsisEvent( TrajQuery( boa, scName, "Enceladus" ), ApsisEvent.PERIAPSIS )
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
        print('... ' + ntpath.basename(periEventFile) + ' Saved!' )
        print('    Peri Events Found: ' + str(numPeris))
    #encBFframe='IAU Enceladus Fixed'

#-----------------------------------------------------------------------------
# Seach for Peris of Sturn with not Enceladus Flyby:
findPeriNoFBs = True
if findPeriNoFBs:
    # TODO: Replace by search intervals below minAlt. then search for Peris.
    # Strategy, search per searchBodies, and do the follow:
    ev = ApsisEvent( TrajQuery( boa, scName, "Saturn" ), ApsisEvent.PERIAPSIS )
    r = ev.search( searchInterval, apsisSearchStep )
    encQuery = TrajQuery( boa, scName, "Enceladus" )
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
        print('... ' + ntpath.basename(periNoFBsEventFile) + ' Saved!' )
        print('    Peri w No FBs Events Found: ' + str(numPerisNoFBs))
    #encBFframe='IAU Enceladus Fixed'


#-----------------------------------------------------------------------------
# Seach for Apoapsis of 'searchBodies':
findApos = True
if findApos:
    ev = ApsisEvent( TrajQuery( boa, scName, "Enceladus" ), ApsisEvent.APOAPSIS)
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
        print('... ' + ntpath.basename(apoEventFile) + ' Saved!' )
        print('    Apo Events Found: ' + str(numApos))

# ----------------------------------------------------------------------------
#Find DV discontinuities:
findDvDiscon = True
#findDvDiscon = False
if findDvDiscon and os.path.exists(dvDiscFile):
    print('... Using DV Disc. from: ' + dvDiscFile )
    with open(dvDiscFile, 'r') as jsonInput:
       dvSearchDic = json.load(jsonInput)

if findDvDiscon and not os.path.exists(dvDiscFile):
    t1_cpu = process_time()
    print( '... Searching for DV discontinuities: ')
    dvSearchDic=[]
    querySat = TrajQuery( boa,scName,dvSearchCenter,dvSearchFrame)
    timeStep = dvSearchtimeStep
    tList =Epoch.range(traj_t0,traj_tf, timeStep)
    numDvDisc = 0
    dvTot = 0
    dvMag_list = []
    for tt in tList[:-1]:
        #dvMag = (querySat.state(tt+timeStep).vel() - querySat.state(tt).vel()).mag()/querySat.state(tt).vel().mag()
        #dvMag = abs(querySat.state(tt+timeStep).vel().mag() - querySat.state(tt).vel().mag())
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

    # Data Save:
    saveDvData = True
    if saveDvData:
        with open(dvDiscFile, 'w+' ) as outfile:
           json.dump( dvSearchDic, outfile, indent = 4, separators=(',', ': ') )
        print('     ' + ntpath.basename(dvDiscFile) + ' Saved!' )


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
# ADD CP's: At t0 and tf of Timeline:
# ----------------------------------------------------------------------------
#TMP: traj_t0 (READ from CONFIG if)
AddCPatT0 = False
if AddCPatT0:
    print('... Adding CP  at traj_t0 (FB type).')
    stEncQuery= M.TrajQuery( boa, scName,'Enceladus','EMO2000',)
    mcpUtil.appendCpFBs(
        mgr,'ETF-01', traj_t0, stEncQuery,
        center='Enceladus',
        frame='EMO2000',
        body='encnf5'
        )
else:
    print('... Ignoring CP at traj_t0.')

#TMP: traj_tf (READ from CONFIG if)

# ----------------------------------------------------------------------------
# ADD CP's: At Peri
# ----------------------------------------------------------------------------
print('... Adding CPs: -> Peri')
#addCPs[0]['source']
cpPeriFile = addCPs[0]['source']
# Check if file exist:
print('... Loading: ' + cpPeriFile )
with open(cpPeriFile, 'r') as jsonInput:
   cpPeriDic = json.load(jsonInput)

#TODO  In JSON files: keep Center and Frame common to all events.
#      create new searches if new valueas are needed
#      Improve the call namne -> searchBodies[0]

#Create query:
stQuery= M.TrajQuery( boa, scName,
                      searchBodies[0]['bodyName'],
                      searchBodies[0]['searchFrame'])

if AddCPatT0:
    EFTini = 2
else:
    EFTini = 1

for ii, cp in enumerate(cpPeriDic):
    cpName = searchBodies[0]['bodyName'][0] + 'TF-' + str(ii+EFTini).zfill(2) # Start at ETF02 (Europa Flyby 02 - ETF: Eenceladus Tour Flyby)
    if M.Epoch(cp['time']) > M.Epoch(traj_tf):
        print('... WARNING: Skipping ' + cpName + '. CP_time > traj_tf')
        continue
    mcpUtil.appendCpFBs(
        mgr,cpName,cp['time'],stQuery,
        center=cp['center'],
        frame=cp['frame'],
        body=scName
        )

# # TODO: query and time. stateParams names, and bounds (numerically)
#
# ----------------------------------------------------------------------------
# ADD DV's at DV_Disc
# ----------------------------------------------------------------------------
print('... Adding DVs:')
cpPeriFile = addDVs[0]['source']
# Check if file exist:
print('... Loading: ' + dvDiscFile )
with open(dvDiscFile, 'r') as jsonInput:
   dvDiscDic = json.load(jsonInput)

for ii, dv in enumerate(dvDiscDic):
    dvName ='DV'+str(ii).zfill(2)
    mcpUtil.addTimeBurn(
        mgr, dvName, dv['time'],
        frame=dv['frame'],
        dvel=M.Dbl3Vec(dv['value']),
        dtBound=5*3600*s,
        dvBound=0.04*km/s,
        )
print('    -> '+str(len(dvDiscDic))+' DVs Added!')

# ----------------------------------------------------------------------------
# ADD DV's at Apo
# ----------------------------------------------------------------------------
print('... Adding DVs at APO (every 10 Apos):')
# Check if file exist:
print('... Loading: ' + addApoDVs[0]['source'] )
with open(addApoDVs[0]['source'], 'r') as jsonInput:
   dvApoDic = json.load(jsonInput)

dvApoEvery= 5
dvsAdded = 0
for ii, apo in enumerate(dvApoDic):
    if ii%dvApoEvery == 0:
        dvsAdded +=1
        dvName ='DVApo'+str(dvsAdded).zfill(2)
        mcpUtil.addTimeBurn(
            mgr, dvName, apo['time'],
            frame=apo['frame'],
            dvel=M.Dbl3Vec([1e-5] * 3),
            dtBound=5*3600*s,
            dvBound=0.01*km/s,
            )
print('    -> '+str(dvsAdded)+' DVs at Apo Added!')


# ----------------------------------------------------------------------------
# ADD DV's at Peri no Enc
# ----------------------------------------------------------------------------
print('... Adding DVs at Peri (every 10 Apos):')
# Check if file exist:
print('... Loading: ' + addPeriDVs[0]['source'] )
with open(addPeriDVs[0]['source'], 'r') as jsonInput:
   dvPeriDic = json.load(jsonInput)

dvPeriEvery= 5
dvsAdded = 0
for ii, peri in enumerate(dvPeriDic):
    if ii%dvPeriEvery == 0:
        dvsAdded +=1
        dvName ='DVPeri'+str(dvsAdded).zfill(2)
        mcpUtil.addTimeBurn(
            mgr, dvName, peri['time'],
            frame=peri['frame'],
            dvel=M.Dbl3Vec([1e-5] * 3),
            dtBound=5*3600*s,
            dvBound=0.01*km/s,
            )
print('    -> '+str(dvsAdded)+' DVs at Peris Added!')

# ----------------------------------------------------------------------------
# ADD CP's: At tf of Timeline:
# ----------------------------------------------------------------------------
#TMP: traj_tf (READ from CONFIG if)
AddCPatTf = True
if AddCPatTf:
    print('... Adding CP  at traj_tf (EOI).')
    stEncQuery= M.TrajQuery( boa, scName,'Enceladus','EMO2000',)
    mcpUtil.appendCpFBs(
        mgr,'EOI', traj_tf, stEncQuery,
        center='Enceladus',
        frame='EMO2000',
        body='encnf5'
        )
else:
    print('... Ignoring CP at traj_tf.')


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
