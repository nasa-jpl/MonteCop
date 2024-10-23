#!/usr/bin/env python
# ===========================================================================
# Section 392 Navigation and Mission Design
#
# Copyright (C) 2021, California Institute of Technology.
# U.S. Government Sponsorship under NASA Contract NAS7-03001 is acknowledged.
#
# ===========================================================================


'''
Project: MonteCop (a Monte-Copernicus interface)

bsp2visualCop.py:  ...

Created By  : Ricardo L. Restrepo ricardo.l.restrepo@jpl.nasa.gov
Created Date: 2022/12/01
version ='1.0'

  Examples:
       >> ./bsp2visualCop EuropaClipper_19F23_V1.bsp  -c Jupiter -f j2000 -sc '-159'
       >> ./bsp2visualCop EuropaClipper_19F23_V1.bsp  -c Jupiter -f j2000 -sc 'Europa Clipper'
       >> ./bsp2visualCop EuropaClipper_19F23_V1.bsp  -c Europa -f iau_body_fixed -sc '-159'
       >> ./bsp2visualCop cot1_19F20_V4.bsp  -sc '-650' -c Jupiter  -bl Jupiter Europa -co red
       >> ./bsp2visualCop cot1_19F20_V4.bsp  -sc '-650' -c Jupiter -to 10 -tl  20 -bl Jupiter Europa
       >> ./bsp2visualCop cot1_19F20_V4.bsp  -sc '-650' -c Europa -f iau_body_fixed  -bl Jupiter Europa -o cot1_EuFrame.ideck -dt 100
       >> ./bsp2visualCop ENC_NF5_2034A4.bsp -c Sun -dt 10000 -bl Sun Earth Jupiter Venus Saturn -sc -54856
       >> ./bsp2visualCop Enceladus_Tour0.bsp -sc -303 -c Saturn -bl Saturn Titan Enceladus Rhea Tethys Dione

 Last update: --

'''


# ==============================================================================
# imports here:

import sys
import os

#sys.path.append('/Users/ricgomez/Documents/Copernicus/Copernicus_v5.1/Copernicus_5.1.0/copernicus/utils/CopPy')
#sys.path.append('/Users/ricgomez/Documents/Copernicus/MonteCop/monteCop520/utils')

# import robocoppy as rcpy
# from robocoppy import ColorEnum as copColor
# import spiceIDs as spiceIDs

import monteCop.src.CopPy510.robocoppy as rcpy
from monteCop.src.CopPy510.robocoppy import ColorEnum as copColor
import monteCop.utils.spiceIDs as spiceIDs

from copy import deepcopy
import json
import argparse
import re

from time import strptime

#-----------------------------------------------------------------------------
# Parse User Inputs:
#-----------------------------------------------------------------------------

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

if args.configFile:
    print(" ... Reading JSON config file (test only)")

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
    print(f" ... t0= {t0_cal['yy']} {t0_cal['mm']} {t0_cal['dd']} {t0_cal['hh']}")
    print(f" ... tf= {tf_cal['yy']} {tf_cal['mm']} {tf_cal['dd']} {tf_cal['hh']}")
else:
    print(' ERROR: data for object '+ scID + ' NOT Found!')
    sys.exit()

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
    t0_julDay = float(args.tiniOffset)
else:
    t0_julDay = tOffsetDays                                # Add a time offset of 1 sec
if args.tlInterval:
    tf_julDay = float(args.tlInterval)
else:
    tf_julDay  = (float(tf_etsec) - float(t0_etsec))/86400.0 - 2*tOffsetDays   # Add a time offset of 1 sec  -> Need to make it days

print(f" ... Seg_t0_julDay = {t0_julDay}")
print(f" ... Seg_tf_julDay = {tf_julDay}")
print(f" Center: {visCenter}")
print(f" Frame: {visFrameType}")
print(f" Body_list: {bodies_list}")

#raise Exception('raise')

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
    trajIdeck.filename = bspFile.replace('.bsp','_VisCOP.ideck')

#load *.bsp in tour (and remove template .bsp)
#trajIdeck.spice.append(bspFile)
for fl in bspFile_list:
    trajIdeck.spice.append(fl)

#load *.bsp in tour (and remove template .bsp)
#for ll in ephem_list: trajIdeck.spice.append(ll)

#Import copUtils
# use color_list to color legs (if multiples)
#
# #-----------------------------------------------------------------
# # Set Frames
# #-----------------------------------------------------------------

bodies_list_ID= [passBodyID(bb) for bb in bodies_list ]

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
#scFrame.mainbody.id     = scID
#scFrame.auxbody.id      = passBodyID(visAuxBody)
if passBodyID(scID):
    scFrame.mainbody.id = passBodyID(scID)
else:
    scFrame.mainbody.id = scID

#-----------------------------------------------------------------
# Force Variables/Graphics: frame, Epoch, badies list
#-----------------------------------------------------------------

# Set force model: Epoch, Frame, center, and Bodiy list
trajIdeck.force_vars = {
    'force_date_choice': 1,           #Calendar form:
    'force_frame':  visFrame,
    'force_bodies': bodies_list_ID,
    'force_year':   t0_cal['yy'],
    'force_month':  strptime(t0_cal['mm'],'%b').tm_mon,  # Convert to number
    'force_day':    t0_cal['dd'],
    'force_hour':   t0_cal['hh'].split(':')[0],
    'force_min':    t0_cal['hh'].split(':')[1],
    'force_sec':    t0_cal['hh'].split(':')[2],
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

#-----------------------------------------------------------------
# Create State Object:
#-----------------------------------------------------------------
stateRef        = rcpy.State()
stateRef.frame  = scFrame

#Defualt state params:R,V. so, not needed to be defined here, but for completeness:
p = stateRef.param
p.params_id = [rcpy.Param1Enum.rx,
              rcpy.Param2Enum.ry,
              rcpy.Param3Enum.rz,
              rcpy.Param4Enum.vx,
              rcpy.Param5Enum.vy,
              rcpy.Param6Enum.vz,
              ]
stateRef.rx = 0.0; stateRef.ry = 0.0; stateRef.rz = 0.0
stateRef.vx = 0.0; stateRef.vy = 0.0; stateRef.vz = 0.0

#-----------------------------------------------------------------
# Set segment(s):
#-----------------------------------------------------------------
numSegs = 1
#segName = subphases[nSp][0]+'_('+colorName_list[nSp]+')'
segName = str(scID)
seg= rcpy.Segment(segName)
seg.state=deepcopy(stateRef)                                #deep Copy reference state
seg.use_global_integrator=False
seg.integration=user_IntegMethod_static
seg.def_central_body.id = scID                           #Clipper
seg.time.t0 = t0_julDay
seg.time.dt = tf_julDay
#seg.plot_data.plot_color = color_list[str(6)]           #Assing seg color
seg.plot_data.plot_color = segColor        #Assing seg color
seg_list=[]
seg_list.append(seg)

# TODO :: Make multiple segments (nSegInput)

#-----------------------------------------------------------------
# Save Ideck:
#-----------------------------------------------------------------
#Load segments:
trajIdeck.segments=seg_list

# --> Mission info:
trajIdeck.globalvars.mission_info='Trajectory generated from the BSP kernel: ' + bspFile + ' and spiceID:  ' + scID

#Graphics:
trajIdeck.graphics={'ogl_enableprintingmode': False,  #use Printer friendly mode color
               'ogl_maxbodytimestep': 0.01, #(days), times step to plot bodies trajecotries
               'grap_axis_length': 1.0E+05,
               }

#-----------------------------------------------------------------
# Save Ideck:
#-----------------------------------------------------------------
#ave ideck:
trajIdeck.save(trajIdeck.filename)
print(' Solution Saved to:' + trajIdeck.filename)

#===============================================================================

# TODO:
# Allows input data to be enter through confi.JSON
# Move it to nexus. Then clone in local, use all dependancies from monteCop module
