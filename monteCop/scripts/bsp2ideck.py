#!/usr/bin/env python
#===========================================================================
#
# Section 392 Navigation and Mission Design
#
# Copyright (C) 2018, California Institute of Technology.
# U.S. Government Sponsorship under NASA Contract NAS7-03001 is acknowledged.
#
#===========================================================================

'''
MonteCop script: loadBSP.py

Decription:   Tranfer between two consecutive Flybys
              reading Flyby Epochs from Json files.
              and Flyby States from BSP. Each sequence encompass 5 segments
              Reapeat sequence N times

Digram:      .B1_fb ---> DV_CU ---> |BP1|<--Apo.Apo-DV -->|BP2|<--B2_fb.

Example:    bsp2ideck.py cot1_19F20_V4.bsp -sc -650 -c Jupiter -bl Jupiter Europa
            bsp2ideck.py 21F31_V4_mdos.bsp -sc -159 -to 2400 -tl 100 -c Jupiter -bl Jupiter Europa -j copConfig_Euclip.json

  last upfate: Nov.10 2023.
  By Ricardo L. Restrepo
'''

# ==============================================================================
# imports here:

import sys
import os

import monteCop.src.CopPy510.robocoppy as rcpy
from monteCop.src.CopPy510.robocoppy import ColorEnum as copColor
import monteCop.utils.spiceIDs as spiceIDs
import monteCop.utils.copUtils as mcp

from copy import deepcopy
import json
import argparse
import re

from importlib import reload    # for debugging. Reload Module: importlib.reload(module)
from time import strptime
from datetime import datetime, timedelta


# ============================================================================
# Parse User Inputs:
# ============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("inputBSPfile", nargs = '+', metavar="inputBSPfile.bsp",
                    help="BSP input file(s). Multiple files in beta mode. \
                    Only fist file is considered, others are added to the SPICE kernel so user can use them")
parser.add_argument("-o", "--outputName",default="", help = "Optional output file name. Default:... ")
parser.add_argument("-sc","--scID", default= "-999", help="Spice ID or Name (e.g. -150 or 'CASSINI PROBE'). Defautl = '-999'  ")
parser.add_argument('-c', '--center', default= "Earth", help = "Body Center. Default: 'Enceladus' ")
parser.add_argument('-f', '--frame', default="j2000", help = "Visualization Frame = 'j2000' (Some options: eclipj2000, iau_body_fixed) ")
parser.add_argument('-bl','--bodyList', default= ["Earth","Sun","Moon"], nargs='+', help= "Body list to be included" )
parser.add_argument("-to","--tiniOffset", default="", help="t0 offset- days. Default='' (use t0 from bsp)")
parser.add_argument("-tl","--tlInterval", default="", help="Trajectory time line interval from t0_offset. Default='' (use bsp time span)")
parser.add_argument('-dt', '--propSteps', default="10", help = "Step time for plotting (in minutes). Default: 10 (min)" )
parser.add_argument('-co', '--color', default="blue", help = "Segment Color(s). Defualt=blue" )
parser.add_argument('-n', '--numSegs', default="1", help = "Number of segments. Defualt=1" )
parser.add_argument('-j', "--configFile",default="", help = "JSON configuration file (input user params through a 'sysConfig.json') ")
parser.add_argument('-ax', '--auxBody', default= "Moon",
                    help = "auxiliar Body (secondary body), only needed on Rot frames. Default: 'Moon' ")

args = parser.parse_args()

# ============================================================================
# Set User Params:
# ============================================================================
scID = args.scID
visCenter = args.center
visFrameType = args.frame
visAuxBody = args.auxBody
propSteps= float(args.propSteps)*60
numSegs= args.numSegs

#TMP fix to accept several BSPs:
# Need to merge bsp files to accept several bsp files
bspFile_list = args.inputBSPfile
bspFile = bspFile_list[0]

bodies_list = args.bodyList

tOffsetDays  = 2.0/86400   # Cut ini and end of BSP to avoid conflics reading data

# ============================================================================
# Read User Configuration file:
# ============================================================================

Epochs_list_Cal = []

if args.configFile:
    print(" Reading JSON config file")

    # Load your JSON file into a dictionary
    with open(args.configFile) as json_file:
        jsonConfig = json.load(json_file)

    # Read Segments Epoch List, if in Config File:
    try:
        Epochs_list_Cal = jsonConfig['Copernicus']['ControlPoints']['Epochs']
        print(' ... Reading CPs Epochs from config file')
    except:
        pass

# ============================================================================
# Lists:
# ============================================================================
#Color list:
color_list={'1' : copColor.orange,
            '2' : copColor.aqua,
            '3' : copColor.dodgerblue,
            '4' : copColor.yellow,
            '5' : copColor.red,
            '6' : copColor.blue,
            '7' : copColor.blueviolet,
            '8' : copColor.peru,
            '9' : copColor.darkcyan,
            '10': copColor.greenyellow,
            }

#Get color Names
colorName_list=[]
for cn in color_list:
     colorName_list.append(str(color_list[str(cn)]).replace('ColorEnum.',''))

#general lists:
stateNodes={'t0m' : 1,
            't0p' : 2,
            'tfm' : 3,
            'tfp' : 4,
            }

timeNodes={'t0' : 1,
           'dt' : 2,
           'tf' : 3,
          }

# ============================================================================
# Util Functions:
# ============================================================================
def passBodyID(bodyID):
    ''' Check and pass SPICE BodyID (accept str or num)
        = Inputs:
        - bodyID    Body name (e.g. "Earth") or ID (e.g., 399)
        = Outputs:
        - bosyID    SPICE body ID
    '''
    # NOTE: body name is upper-case on SpiceBodyName
    if type(bodyID) == str and bodyID.upper() in spiceIDs.SpiceBodyName.keys():
        return spiceIDs.SpiceBodyName[bodyID.upper()]
    elif bodyID in spiceIDs.SpiceBodyName.values():
        return bodyID
    else:
        #sys.exit('Invalid SPICE Body Name or ID: '+ str(bodyID))
        #print('Invalid SPICE Body Name or ID: '+ str(bodyID))
        return []

def setSegment(segIn,centralBodyId,stateType,stateFrame,stateIn=[]):
    segIn.def_central_body.id=centralBodyId
    segIn.state=deepcopy(stateType)
    segIn.state.frame=stateFrame
    if stateIn:
        for ii in range(6):
            segIn.state[ii].value=stateIn[ii]

def setStateOV(xIn,bounds=[],scaleIn=[],Dpert=[]):
    xIn.ov=True
    if scaleIn: xIn.scale=scaleIn
    if Dpert:   xIn.dpert=Dpert
    if bounds:
        xIn.ge_check=True
        xIn.ge_value=bounds[0]
        xIn.le_check=True
        xIn.le_value=bounds[1]
        # xIn={'le_check':True, 'le_value': bounds[0],
        #      'ge_check':True, 'ge_value': bounds[1],
        #     }

def inheritState(segInhering,segInhereted,nodeIn):
    for xx in range(6):
        segInhering.state[xx].assume      = 2  #inherit ( value(1),inherit(2),DValue(3))
        segInhering.state[xx].node        = stateNodes[nodeIn]  #(t0m,t0p,tfm,tfp)
        segInhering.state[xx].inherit_seg = segInhereted
    return segInhering


def xStateConst(xIn,segTarget,nodeIn,scaleIn=1.0,bounds=[]):
    xIn={'eq_check': True,
         'eq_node' : stateNodes[nodeIn],  #(t0m,t0p,tfm,tfp)
         'eq_seg'  : segTarget,
         'scale'   : scaleIn,
         }
    return xIn

def setStateConst(segIn,segTarget,nodeIn,constFrame,posScaleIn=1E5,velScaleIn=1E3):
    #Extend: add segIn stateNode: (t0m,t0p,tfm,tfp)
    #xFun=segIn.functions.state.t0m
    #xFun=segIn.functions.state.t0p
    #xFun=segIn.functions.state.tf0
    xFun=segIn.functions.state.tfp
    xFun.rx=xStateConst(xFun.rx,segTarget,nodeIn,posScaleIn)
    xFun.ry=xStateConst(xFun.ry,segTarget,nodeIn,posScaleIn)
    xFun.rz=xStateConst(xFun.rz,segTarget,nodeIn,posScaleIn)
    xFun.vx=xStateConst(xFun.vx,segTarget,nodeIn,velScaleIn)
    xFun.vy=xStateConst(xFun.vy,segTarget,nodeIn,velScaleIn)
    xFun.vz=xStateConst(xFun.vz,segTarget,nodeIn,velScaleIn)
    xFun.frame=constFrame


def setDV0(segIn,frameIn,valIn=[1e-4],setOV=False,setObj=False,mnvFrame=[]):
    if len(valIn) == 1: valIn.extend([valIn[0],valIn[0]])
    #Extend: DVf
    segIn.dv0.active = True
    segIn.dv0.mnvrframe.reference_frame=frameIn
    #NOTE: DV0(0):DVx, DV0(1):DVy, DV0(2):DVz, DV0(4):DV_mag, DV0(5):DV_alpha, etc...
    segIn.dv0.engine.isp.use=True           #Use ISP (ISP, CEV,Thrust)
    segIn.dv0.x={'value':valIn[0],'use':True}
    segIn.dv0.y={'value':valIn[1],'use':True}
    segIn.dv0.z={'value':valIn[2],'use':True}
    if setOV:
        segIn.dv0.x.ov=True
        segIn.dv0.y.ov=True
        segIn.dv0.z.ov=True
    if setObj:
        segIn.functions.dv0.mag.obj=True
        segIn.functions.dv0.frame.controls_frame_id=1  #ijk
        if mnvFrame:
            segIn.functions.dv0.frame.reference_frame=mnvFrame
        else:
            segIn.functions.dv0.frame.reference_frame=frameIn


# def inheritTime(segInhering,segInhereted,nodeIn):
#     segInhering.time={'use': True, 'assume':2, 'node':3,
#                       'inherit_seg':segInhereted}

# ============================================================================
# Dictionaries:
# ============================================================================
#Color list:
color_list={'1' : copColor.orange,
            '2' : copColor.aqua,
            '3' : copColor.dodgerblue,
            '4' : copColor.yellow,
            '5' : copColor.red,
            '6' : copColor.blue,
            '7' : copColor.blueviolet,
            '8' : copColor.peru,
            '9' : copColor.darkcyan,
            '10': copColor.greenyellow,
        }

# Can add colors by:
segColor = getattr(copColor,args.color)


# ============================================================================
# Scan BSP:
# ============================================================================
#Read BSP with brief to obtain t0 and tf for scID (need options -t, and -n):

    # NEED TO REDO THE TIMING HANDLING
    # Use a single timing convention, currently mixing all type

# Need to merge bsp files to accept several bsp files:
tmpFile = 'tmp_bspInfo.txt'
os.system('brief ' +  bspFile +  ' -t -n -sec > '+ tmpFile )
tmpData = open(tmpFile, "r")
ll = tmpData.readlines()
tmpData.close()
os.system('rm '+ tmpFile)
for ln in ll:
    #Find line with scID info:
    if len(ln) > 2 and ln.split()[0] == scID:
        li = ln.split()
        t0_cal ={'yy':li[1],'mm':li[2],'dd':li[3],'hh':li[4]}
        tf_cal ={'yy':li[5],'mm':li[6],'dd':li[7],'hh':li[8]}

 #Check bsp info was found
if 't0_cal' in locals():
    print(' Data for object ' + scID + ' Found!')
    print(f" ... t0 (BSP)= {t0_cal['yy']} {t0_cal['mm']} {t0_cal['dd']} {t0_cal['hh']}")
    print(f" ... tf (BSP)= {tf_cal['yy']} {tf_cal['mm']} {tf_cal['dd']} {tf_cal['hh']}")
else:
    print(' ERROR: data for object '+ scID + ' NOT Found!')
    sys.exit()

# Using Calendar Format 2
t0_cal_F2 =t0_cal['dd']+'-'+t0_cal['mm']+'-'+t0_cal['yy']+' '+t0_cal['hh']
tf_cal_F2 =tf_cal['dd']+'-'+tf_cal['mm']+'-'+tf_cal['yy']+' '+tf_cal['hh']

Ideck_Epoch_JD = mcp.date_to_julian(t0_cal_F2)

# Get t0_etsec, tf_etsec: Replace or addapt scripts to read and/or Cal, etsec
os.system('brief ' +  bspFile +  ' -t -n -etsec -sec > '+ tmpFile )
tmpData = open(tmpFile, "r")
ll = tmpData.readlines()
tmpData.close()
os.system('rm '+ tmpFile)
for ln in ll:
    #Find line with scID info:
    if len(ln) > 2 and ln.split()[0] == scID:
        li = ln.split()
        t0_etsec = li[1]
        tf_etsec = li[2]

if args.tiniOffset:
    t0_days = float(args.tiniOffset)
else:
    t0_days = tOffsetDays                                # Add a time offset of 1 sec
if args.tlInterval:
    tf_days = t0_days+  float(args.tlInterval)
else:
    tf_days  =t0_days + (float(tf_etsec) - float(t0_etsec))/86400.0 - 2*tOffsetDays   # Add a time offset of 1 sec  -> Need to make it days

print(f" ... Seg_t0_Days = {t0_days}")
print(f" ... Seg_tf_Days = {tf_days}")
print(f" Center: {visCenter}")
print(f" Frame: {visFrameType}")
print(f" Body_list: {bodies_list}")


# # ============================================================================
# # Initialize trajIdeck:
# # ============================================================================
#
#-----------------------------------------------------------------
# Create Ideck:
#-----------------------------------------------------------------
trajIdeck=rcpy.Ideck()
if args.outputName:
    trajIdeck.filename = args.outputName
else:
    trajIdeck.filename = bspFile.replace('.bsp','_BSP2COP.ideck')

#load *.bsp in tour (and remove template .bsp)
#trajIdeck.spice.append(bspFile)
for fl in bspFile_list:
    trajIdeck.spice.append(fl)

#-----------------------------------------------------------------
# Set Frames
#-----------------------------------------------------------------

# Create J200 frame
j2k = rcpy.Frame()
j2k.frametype_id    = rcpy.CopFrameEnum.j2000    # 1       #-> SPICE: EMO2000
j2k.framecenter_id  = 1                          # main-> 1, aux -> 2, barycenter -> 3, L1 ->4, ...
j2k.mainbody.id     = passBodyID(visCenter)
j2k.auxbody.id      = passBodyID(visAuxBody)

# Set Visualization Frame
visFrame = rcpy.Frame()
visFrame.frametype_id    = getattr(rcpy.CopFrameEnum, visFrameType)   # 1
visFrame.framecenter_id  = 1   #main-> 1, aux -> 2
visFrame.mainbody.id     = passBodyID(visCenter)
visFrame.auxbody.id      = passBodyID(visAuxBody)

# Set spacecraft frame
#scFrame = deepcopy(visFrame)
scFrame = rcpy.Frame()
scFrame.frametype_id    = rcpy.CopFrameEnum.body_fixed  # -40
scFrame.framecenter_id  = 1   #main-> 1, aux -> 2
if passBodyID(scID):
    scFrame.mainbody.id = passBodyID(scID)
else:
    scFrame.mainbody.id = scID

# ----------------------------------------------------------------
# Celestial Bodies list:
bodies_list_ID= [passBodyID(bb) for bb in bodies_list ]


#-----------------------------------------------------------------
# Force Variables/Graphics: frame, Epoch, badies list
#-----------------------------------------------------------------
# Set force model: Epoch, Frame, center, and Bodiy list

# # Use Calendar Epoch
# trajIdeck.force_vars = {
#     'force_date_choice': 1,           #Calendar form:
#     'force_frame':  visFrame,
#     'force_bodies': bodies_list_ID,
#     'force_year':   t0_cal['yy'],
#     'force_month':  strptime(t0_cal['mm'],'%b').tm_mon,  # Convert to number
#     'force_day':    t0_cal['dd'],
#     'force_hour':   t0_cal['hh'].split(':')[0],
#     'force_min':    t0_cal['hh'].split(':')[1],
#     'force_sec':    t0_cal['hh'].split(':')[2],
# }

# Use Julian Date Epoch:
trajIdeck.force_vars = {'force_date_choice': 2,           #Set to Julian Date:
                    'force_frame':  visFrame,
                    'force_bodies': bodies_list_ID,
                    'force_julian': Ideck_Epoch_JD,
                    'force_time_system_cal': 'TDB',    #optional: just to be sure
                    'force_time_system_jd': 'JDTDB',  #optional: just to be sure
                    'force_time_system_et': 'TDB',    #optional: just to be sure
                    }

#Assign frame and bodies to Graphic:
trajIdeck.graphics.vis_frame      = visFrame
trajIdeck.graphics.bodies_to_plot = bodies_list_ID

#-----------------------------------------------------------------
# Global Integrator, and Integration Methods: Set Static
#-----------------------------------------------------------------
trajIdeck.globalvars.global_integrator=1

#TODO: propSteps from user (in sec)
#user Defined Integration method: Statitc
user_IntegMethod_static={'method': rcpy.IntegMethod.static,   #same as:1003
                        'output_time_step_choice': 2,
                        'time_step_out': int(propSteps),
                        'n_steps': 10,
                        }

# --> Mission info:
trajIdeck.globalvars.mission_info='Trajectory generated from the BSP kernel: ' + bspFile + ' and spiceID:  ' + scID

#Graphics:
trajIdeck.graphics={'ogl_enableprintingmode': False,  #use Printer friendly mode color
               'ogl_maxbodytimestep': 0.01, #(days), times step to plot bodies trajecotries
               'grap_axis_length': 1.0E+05,
               }
#-----------------------------------------------------------------
# Create State Object:
#-----------------------------------------------------------------
nullState        = rcpy.State()
nullState.frame  = scFrame

#Defualt state params:R,V. so, not needed to be defined here, but for completeness:
p = nullState.param
p.params_id = [rcpy.Param1Enum.rx,
              rcpy.Param2Enum.ry,
              rcpy.Param3Enum.rz,
              rcpy.Param4Enum.vx,
              rcpy.Param5Enum.vy,
              rcpy.Param6Enum.vz,
              ]
nullState.rx = 0.0; nullState.ry = 0.0; nullState.rz = 0.0
nullState.vx = 0.0; nullState.vy = 0.0; nullState.vz = 0.0


#Define State -> State X (R,V):
#Defualt state params:R,V. so, not needed to be defined here, but for completeness:
stateX        = rcpy.State()
stateX.param.params_id = [
            rcpy.Param1Enum.rx,
            rcpy.Param2Enum.ry,
            rcpy.Param3Enum.rz,
            rcpy.Param4Enum.vx,
            rcpy.Param5Enum.vy,
            rcpy.Param6Enum.vz,
            ]

#-----------------------------------------------------------------
# Epoch List (-> Pass it as JSON )
#-----------------------------------------------------------------

#raise Exception('exit')

if not Epochs_list_Cal:
    Epochs_list_Cal = [mcp.julian_to_date(Ideck_Epoch_JD +t0_days),
                       mcp.julian_to_date(Ideck_Epoch_JD +tf_days)]

# Convert to Julian dates
Epochs_list_JD = [ mcp.date_to_julian(ii) for ii in Epochs_list_Cal]

# Create Segment Epochs in Days from Ideck Epoch:
Segments_Epochs = [ii - Ideck_Epoch_JD for ii in Epochs_list_JD]

#-----------------------------------------------------------------------

# Null segmets: Read states from BSP as static segments
segNull_list=[]
for ii in range(len(Segments_Epochs)-1):
    segName = 'seg'+str(ii+1) #Create segment name
    segN= rcpy.Segment(segName)                                  #Instanciate a segment
    segN.state=deepcopy(nullState)                               #deep Copy reference state
    segN.use_global_integrator=False
    segN.integration=user_IntegMethod_static
    segN.def_central_body.id = scID                           #Clipper
    segN.def_central_body.id = passBodyID(visCenter)
    segN.time.t0  = Segments_Epochs[ii]
    segN.time.dt  = 0.0                                     # Just need state at t0
    segN.plot_data.plot_color   = color_list['1']    #Assing seg color
    segNull_list.append(segN)

# State segmets: Read states null segments by inherent
segState_list=[]
for ii in range(len(Segments_Epochs)-1):
    colorId=(ii+1)%10
    segName = 'seg'+str(ii+1) #Create segment name
    segN= rcpy.Segment(segName)                                  #Instanciate a segment
    segN.state=deepcopy(stateX)                               #deep Copy reference state
    segN.state.frame=j2k
    segN.def_central_body.id = passBodyID(visCenter)
    segN = inheritState(segN, segNull_list[ii], 't0m')
    segN.time.t0  = Segments_Epochs[ii]
    #segN.time.dt  = Segments_Epochs[ii+1]-Segments_Epochs[ii]     #SPEED: Can make dt = 0, then add times
    segN.time.dt  = 0.0                                            #SPEED: Can make dt = 0, then add times
    segN.plot_data.plot_color   = color_list[str(colorId)]        #Assing seg color
    segState_list.append(segN)

#Load segments into Ideck:
trajIdeck.segments = segNull_list + segState_list

# Propagate ideck to load state info into StateSegments from Null segments
print(' Propagating Ideck: ... ')
x = trajIdeck.propagate()
for ii in range(len(trajIdeck.segments)):
    seg = trajIdeck.segments[ii]
    for xx in range(6):              # Change 'inherent' by 'value'
        seg.state[xx].assume = 1     #inherit ( value(1),inherit(2),DValue(3))

# Remove Null segments from ideck:
del trajIdeck.segments[:len(segNull_list)]

# Add Propagation time to Segments
for ii in range(len(Segments_Epochs)-1):
    segN = trajIdeck.segments[ii]
    segN.time.dt  = Segments_Epochs[ii+1]-Segments_Epochs[ii]

# TODO:
# Split segments to be back and forward prop.
# Add constraints to make segments continuous.
# Implement: B1_fb ---> DV_CU ---> |BP1|<--Apo.Apo-DV -->|BP2|<--B2_fb. Scheme


#-----------------------------------------------------------------
# Save Ideck:
#-----------------------------------------------------------------
# #save ideck:
trajIdeck.save(trajIdeck.filename)
print(' Solution Saved to:' + trajIdeck.filename)

#===============================================================================

#raise Exception('exit')

# TODO:
# Create null states from BSP, given epochs, then convert them to states as:
# - Read State from bsp as is, but save given state parametrization
#   and Frame on Config.json
# - Read states at List of epochs given in config.json

# Need then some how to remove inhirent from bsp, and just make segments "value"
