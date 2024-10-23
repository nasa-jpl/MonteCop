#!/usr/bin/env mpython_q

# ===========================================================================
# Section 392 Navigation and Mission Design
#
# Copyright (C) 2021, California Institute of Technology.
# U.S. Government Sponsorship under NASA Contract NAS7-03001 is acknowledged.
#
# ===========================================================================

'''
MonteCop script:

  Decription:  Generate a Copernicus ideck that load and plot the trajectory
               ...
  Usage:  ...


  last upfate: July 20 2020.
  by:          Ricardo L. Restrepo (ricardo.l.restrepo@jpl.nasa.gov)
'''

#==============================================================================
#imports here:

import Monte as M
import os
import json
import argparse
from copy import deepcopy
from mpy.opt.cosmic import Manager

import monteCop.src.CopPy510.robocoppy as rcpy
from monteCop.src.CopPy510.robocoppy import ColorEnum as copColor

import monteCop.utils.copUtils as mcp
from   monteCop.utils.cosmic2json_mcp import cosmic2json

#For debugging:
from importlib import reload

#==============================================================================


#==============================================================================
# Create and Populate parser:
#==============================================================================
parser    = argparse.ArgumentParser(description =
                    'Generate a Copernicus ideck from an cosimc Timeline. ')
parser.add_argument("inputFile", metavar="cosmicFile.py",
                    help='Cosmic Input file (Timeline)')
parser.add_argument("jsonConfig", default = None, nargs='?',
                    help='json config file with predefined setup information')
parser.add_argument('-o', '--outputDirPath',  metavar='OutputPath',
                    default='./', help = 'Optional output directory. Default: .')
parser.add_argument('-n', '--ideckName',  metavar='ideckName',
                     help = 'Output ideck file name. If none, cosmicFile name is used')
parser.add_argument('-s', action='store_true', help='Save json file solution')

args = parser.parse_args()

#==============================================================================
#Load and read files :
#==============================================================================
cosmicFile=args.inputFile
baseName=os.path.basename(cosmicFile)
caseName=baseName[:-3]
outputDir = args.outputDirPath
saveJson = args.s


#==============================================================================
# Paths and Files Name:
#==============================================================================
mcpDir         = '/home/ricgomez/lib/monteCop/'
mcpUtils       = mcpDir+'utils/'
mcpTemp        = mcpDir+'templates/'
jsonFileOut    = caseName+'_mcp.json'


jsonConfigFile = args.jsonConfig if args.jsonConfig else mcpTemp+'jsonConfig_mcp.json'
ideckFileName  = args.ideckName if args.ideckName else caseName+'.ideck'

#==============================================================================
#Global Params:
#==============================================================================

#SPICE ID's :
euclip_ID  =-159
Europa_ID  = 502
Jupiter_ID = 599

#Jupiter_ID = mcp.passBodyID('Jupiter')
#Europa_ID = mcp.passBodyID('Europa')
#euclip_ID = mcp.passBodyID('Europa Clipper')

# bp_timeFactor -> where to add the bp between cp_i and cp_i+1
bp_tFactor = 0.5
bp_tFactorBack = 1.0 - bp_tFactor

#==============================================================================
# Problem initialization:
#==============================================================================

#-----------------------------------------------------------------------
# Generate JSON format from cosmic file:
#Load cosmic timeline and mcpConfig  into a dictionary: solJson
solJson, mgr= cosmic2json(cosmicFile, jsonConfigFile, saveToFile = saveJson)
copSetup = solJson['Copernicus']
cpList    = solJson['ControlPoints']
impulseMvrs = solJson['ImpulseBurns']

#Set First and final CP's:
cpBeginID = copSetup['cpBegin'] if copSetup['cpBegin'] else 0
cpEndID   = copSetup['cpEnd'] if copSetup['cpEnd'] else (len( solJson['ControlPoints'])-1)

#Set Timeline Epochs:
tlBegin = M.Epoch(solJson['ControlPoints'][cpBeginID]['Time'])
tlEnd   = M.Epoch(solJson['ControlPoints'][cpEndID]['Time'])

#--- Print for debugging---->
print(tlBegin)
print(tlEnd)
print('-------------------------------------------------------')
#Set ideckEpoch (t0) of first control point
#TODO: consider user option.
ideck_t0_JD=mgr.tl.controlPoint(0).time().julianDate('ET')
ideck_tf_JD=mgr.tl.interval().end().julianDate('ET')

# More set up from config File:
maxColorsToUse = min(len(mcp.colorList),copSetup['maxColorsToUse'])

#==============================================================================
# Utils functions (That should go on mcpUtils, but no working)
#==============================================================================
# create Copernicus frame from monte name:
def m2cFrame(monteFrame, mainBodyId="Earth", auxBody="Moon"):
    newFrame = rcpy.Frame()
    newFrame.frametype_id = mcp.passFrameID(monteFrame)
    newFrame.mainbody.id  = mcp.passBodyID(mainBodyId)
    newFrame.auxbody.id   = mcp.passBodyID(auxBody)
    newFrame.framecenter_id = 1
    return newFrame

def setSegment(segIn,centralBodyId,paramsSet,stateFrame,stateIn=[]):
    segIn.def_central_body.id=centralBodyId
    segIn.state=deepcopy(paramsSet)
    segIn.state.frame=stateFrame
    if stateIn:
        for ii in range(6):
            segIn.state[ii].value=stateIn[ii]

def inheritState(segInhering,segInhereted,nodeIn):
    for xx in range(6):
        segInhering.state[xx].assume      = 2  #inherit ( value(1),inherit(2),DValue(3))
        segInhering.state[xx].node        = mcp.stateNodes[nodeIn]  #(t0m,t0p,tfm,tfp)
        segInhering.state[xx].inherit_seg = segInhereted
    return segInhering

# def inheritTime(segInhering,segInhereted,nodeIn):
#     for xx in range(6):
#         segInhering.state[xx].assume      = 2  #inherit ( value(1),inherit(2),DValue(3))
#         segInhering.state[xx].node        = mcp.stateNodes[nodeIn]  #(t0m,t0p,tfm,tfp)
#         segInhering.state[xx].inherit_seg = segInhereted
#     return segInhering

def getStateParams(cpObj):
    if cpObj['StateParamsType'] == 'apoStateParams':
        paramsSet = apoStateParams
    elif cpObj['StateParamsType'] == 'flybyStateParams':
        paramsSet = flybyState
    else:
        sys.exit('ERROR: Incorrect stateParamsType for Control Point: ' +
                 str(cpObj['Name']))
    stateOut=[]
    if cpObj['StateParamsType'] == 'apoStateParams':
            stateOut= [
                M.UnitDbl(cpObj['State'][0]).convert('sec'),  # period(sec)
                M.UnitDbl(cpObj['State'][1]).convert('km'),   # rp(km) periapsisRange
                M.UnitDbl(cpObj['State'][2]).convert('deg'),  # inc(deg)
                M.UnitDbl(cpObj['State'][3]).convert('deg'),  # raan(deg)
                M.UnitDbl(cpObj['State'][4]).convert('deg'),  # aop(deg)
                M.UnitDbl(cpObj['State'][5]).convert('deg'),  # ta (deg) trueAnomaly
            ]
    elif cpObj['StateParamsType'] == 'flybyStateParams':
            stateOut= [
                M.UnitDbl(cpObj['State'][0]).convert('km/sec'),  # vin (km/sec)
                M.UnitDbl(cpObj['State'][1]).convert('deg'),  # ra_vinf(deg)
                M.UnitDbl(cpObj['State'][2]).convert('deg'),  # dec_vinf (deg)
                M.UnitDbl(cpObj['State'][3]).convert('deg'),  # btheta (deg)
                M.UnitDbl(cpObj['State'][4]).convert('km'),   # rp (km) periapsisRange
                M.UnitDbl(cpObj['State'][5]).convert('deg'),  # ta (deg) trueAnomaly
            ]

    return paramsSet, stateOut

#==============================================================================
# Define Param Sets   (to be move to mcpUtils)
#==============================================================================

#------> State OE (Classic Orbtila Elements):
apoStateParams = rcpy.State()
apoStateParams.param.params_id = [
            rcpy.Param1Enum.period,
            rcpy.Param2Enum.rp,
            rcpy.Param3Enum.inc,
            rcpy.Param4Enum.raan,
            rcpy.Param5Enum.aop,
            rcpy.Param6Enum.ta,
            ]
apoStateParams.param.angle_unit = rcpy.AngleUnits.deg  # (need ot be defined)

#------> State OEH (Hyperbolic Orbtila Elements):
flybyState = rcpy.State()
flybyState.param.params_id = [
            rcpy.Param1Enum.vin,        #vInfinity
            rcpy.Param2Enum.ra_vinf,     #inboundRA
            rcpy.Param3Enum.dec_vinf,   #inboundDec
            rcpy.Param4Enum.btheta,     #bPlaneTheta
            rcpy.Param5Enum.rp,         #periapsisRange
            rcpy.Param6Enum.ta,         #TRUE-ANOMALY (careful, NOT timeFromPeriapsis)
            ]
flybyState.param.angle_unit = rcpy.AngleUnits.deg  # (need ot be defined)


#==============================================================================
# Create Ideck and set it up: (Maybe move down)
#==============================================================================

myIdeck = rcpy.Ideck()
myIdeck.filename=ideckFileName


# add *.bsp files here
# myIdeck.spice.append(trajbsp)

# Ideck force Frame
ideckFrame = m2cFrame(copSetup['forceFrame'],copSetup['forceCenter'])

# Force_vars
myIdeck.force_vars = {
    'force_date_choice': 2,  #Set to Julian Date
    'force_frame':  ideckFrame,
    'force_julian': ideck_t0_JD,
    'force_time_system_cal': 'TDB',
    'force_time_system_jd': 'JDTDB',
    'force_time_system_et': 'TDB',
    'force_bodies': copSetup['GravityBodies'],
    }

#Assign frame and bodies to Graphic:
myIdeck.graphics.vis_frame      = ideckFrame
myIdeck.graphics.bodies_to_plot = copSetup['PlotBodies']

#Graphics (defaults):
myIdeck.graphics={
    'ogl_enableprintingmode': False,  #use Printer friendly mode color
    'ogl_maxbodytimestep': 0.01,      #(days),times step to plot body traj.
    'grap_axis_length': 1.0E+05,
    }

#==============================================================================
# Convert Control Points into Segments:
#==============================================================================

#-------------------------------------------------------------------------
# Set iSeg (Use p++ for forward propagation, p-- for backward propagation)

seg_list=[]
for cp_ID in range(cpBeginID,cpEndID + 1):
    # --> Instanciate Segment:
    segNameBase = str(cpList[cp_ID]['Name']) + '_p++'
    iSeg = rcpy.Segment(segNameBase)

    #--> Set time:
    iSeg_t0_JD = Epoch(cpList[cp_ID]['Time']).julianDate('ET')
    iSeg.time.t0 = iSeg_t0_JD - ideck_t0_JD     # days
    if 'TIME' in cpList[cp_ID]['ControlParams'] :
        iSeg.time.t0.ov=True
    # get dt as (segN1_t0-segN2_t0)/2 --> use half tfly (TODO: user defined option)
    if cp_ID < cpEndID :
        iiSeg_t0_JD = Epoch(cpList[cp_ID + 1]['Time']).julianDate('ET')
        iSeg.time.dt = (iiSeg_t0_JD - iSeg_t0_JD) * bp_tFactor
    else:
        iSeg.time.dt = 0.0
    iSeg.time.tf.use = False   # Not needed, default set {t0,dt} but as remainder

    # --> Set State:
    # func-> def setSegment(segIn,centralBodyId,stateType,stateFrame,stateIn=[]):
    centerID = mcp.passBodyID(cpList[cp_ID]['Center'])
    segFrame = m2cFrame(cpList[cp_ID]['Frame'],centerID)
    paramsSet, segState = getStateParams(cpList[cp_ID])
    setSegment(iSeg, centerID, paramsSet ,segFrame,segState)

    # Set iSeg miscellaneus:
    colorID=str(cp_ID%maxColorsToUse)
    iSeg.plot_data.plot_color = rcpy.colors.CopColor(mcp.colorList[colorID].value)
    iSeg.name=segNameBase+'_('+mcp.colorList[colorID].name +')'

    # --> Set Backwad prop Segment: (skip first cp):
    if cp_ID > cpBeginID:
        iSegBack = deepcopy(iSeg)
        iSegBack = inheritState(iSegBack, iSeg, 't0m')
        iSegBack.time.t0 = {'assume':2,    # ( value(1),inherit(2),DValue(3))
                            'node': mcp.timeNodes['t0'],
                            'inherit_seg': iSeg }
        preSeg_t0_JD = Epoch(cpList[cp_ID - 1]['Time']).julianDate('ET')
        iSegBack.time.dt = (preSeg_t0_JD - iSeg_t0_JD) * bp_tFactorBack
        iSegBack.plot_data.plot_color = rcpy.colors.CopColor(
            mcp.colorDarkList[colorID].value
        )
        iSegBack.name = segNameBase.replace('p++','p--') + '_('+mcp.colorDarkList[colorID].name +')'

        # Add Seg to segement list (segBack first):
        seg_list.append(iSegBack)

    # Add Seg to segement list:
    seg_list.append(iSeg)

    # --> Add maneuvers:
    # find maneuvers between iSeg and iiSeg by epoch:
    iDV = [
        dv for dv in impulseMvrs
        if Epoch(cpList[cp_ID]["Time"])  < Epoch(dv["Time"])  < Epoch(cpList[cp_ID + 1]["Time"])
    ]

    if iDV:
        print('Dv: '+str(iDV[0]['Name'])+ ' From '+str(iDV[0]['Start']['Name']))
        # print(Epoch(iDV[0]['Time']))
        # print(Epoch(iDV[0]['Time']).julianDate('ET'))

    # if len(iDv) > 1 --> Warning!!!

    #dviSeg center?. use same is coming from (not ideal but...)

    # if is eventBase it makes it easier -> find event. ??
    # impulseMvrs
    # find if there is a manuever between segi and seg2.
    # (check break point) --> do breakpoint between iSegDv and iiSeg

    # if iDV --> add iSegDV

# Set Segments 2 (Id=1) to N-1

# for others, create copy, and run time backward. Need only inherent state of seg


# starte with CP 0 -> only forward prop.
 # for next build forward and backeard
 # (Can start without bounds ;) )  -Do time and segment logic.


# Add DV's




#==============================================================================
# Add Segments to Ideck:
#==============================================================================
#seg_list=[iSeg]
myIdeck.segments=seg_list

#-----------------------------------------------------------------
# Save Ideck:
#-----------------------------------------------------------------
myIdeck.save(outputDir + ideckFileName)



#===============================================================================
#===============================================================================

# WARNING: cehck on Monte that longitudeOfNode == raan

#TODO: decide how to degin ParamsSet: by orbitType?. by inputs?

#cpFrame:

#flybyFrame

#flybyCenter

#conicType?:

#hyperbolic?

# ------------------------------------------------------
# Force Variables and Epoch:

raise Exception('exit')

# ------------------------------------------------------
# Global Integrator, Integration Methods and Optimizer:


#------------------------------------------------------
# Graphics:



#Get control Points


#Build Segmets:

    #->

    #->

    #->


#==============================================================================
# End
#==============================================================================


#===============================================================================
#TODO:
#===============================================================================

# 1) cosmic2json_mcp to a function.


#===============================================================================
# NOTE:
#===============================================================================

#Get CP names:
#cpNames=mgr.cp.keys()
#or:
#cpNames= [ xx['Name'] for xx in solJson['ControlPoints']


#Get States:
#mgr.cp['E01'].state()
#or
#
#solJson['ControlPoints'][0]['State']


#Test:
#userSetup=solJson['UserSetup']
