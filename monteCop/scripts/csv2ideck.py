#!/usr/bin/env mpython_q

# ===========================================================================
# Section 392 Navigation and Mission Design
#
# Copyright (C) 2021, California Institute of Technology.
# U.S. Government Sponsorship under NASA Contract NAS7-03001 is acknowledged.
#
# ===========================================================================

#==============================================================================
#
# Project: MonteCop (a Monte-Copernicus interface)
#
# genIdeck_fromICsFile.py:
#
# Created By  : Ricardo L. Restrepo ricardo.l.restrepo@jpl.nasa.gov
# Created Date: 2022/02/20 ..etc
# version ='0.1'
#
# Last update:
#
#==============================================================================
# ---------------------------------------------------------------------------
""" genIdeck_fromICsFile.py

    Description:
        generate Ideck from Initial Conditions (X,dt) from *.csv file.
        I.C. are asumme on km,km/s,s units

    Inputs:
        ResonanceFile.csv -> File with family of resonances
            XX = X [km]	Y [km]	Z [km]	DX [km/s]	DY [km/s]	DZ [kms]
            Columns 10-15, dt (Period,[s]) -> column 16
        Id resonance: resnace id
        inclination:  inclination of resonant orbit

    Use Example:
        ./genIdeck_fromICsFile.py  resFile -r 2 -i 50.5

"""
#==============================================================================
#imports here:

#import CopPy.robocoppy as rcpy

import monteCop.src.CopPy510.robocoppy as rcpy
from monteCop.src.CopPy510.robocoppy import ColorEnum as copColor
from copy import deepcopy

#from   CopPy.robocoppy import ColorEnum as copColor
#import utils.mcpUtils as mcp
#from   utils.cosmic2json_mcp import cosmic2json

import os
import json
import argparse
import csv


#For debugging:
from importlib import reload

#==============================================================================


#==============================================================================
# Create and Populate parser:
#==============================================================================
parser    = argparse.ArgumentParser(description =
                    'Generate a Copernicus ideck from an File with State I.C.  and dt.')
parser.add_argument("inputFile", metavar="ICsFile.py",
                    help='csv input file with inicila conditions')
parser.add_argument('-r', '--resId', type=int, default=1,
                    help = " Orbit id. Values start at 2")
parser.add_argument('-i', '--resInc', default = None,
                    help = "inclination of Resonant orbit. If selected, resId is ignored")
parser.add_argument('-o', '--outputDirPath',  metavar='OutputPath',
                     default='./', help = 'Optional output directory. Default: .')
# parser.add_argument('-s', action='store_true', help='Save json file solution')

args = parser.parse_args()


#==============================================================================
# Setup Params:
#==============================================================================

# inc_colum = 18
# dt_colum = 15
# state_colum = 9

#SPICE IDs:
Enceladus_ID = 602

#Global t0:
global_t0_JD = 2.468162500000000E+06    # July 1st, 2045

#==============================================================================
#Read file and get Get I.C. from resId or Inc :
#==============================================================================

ICsFile = args.inputFile
baseName = os.path.basename(ICsFile)[:-4]
outputDir = args.outputDirPath
resId = int(args.resId)
#outputIdeck = baseName + '_ID' + str(resId).zfill(3) + '.ideck'

# Read CSV file:
resData = []
with open(ICsFile) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        resData.append(row)

# TODO: Remove header on resData
# resData.pop[0]

# If args.resInc: Select Resonant orbit by Inclination
if args.resInc:
    resInc = float(args.resInc)
    incDiffList = [[ii, abs(float(row[18]) - resInc)] for ii, row in enumerate(resData[1:])]
    id,val = min(incDiffList,key=lambda x: x[1])
    #use 'lambda function directly to find closest inc.'
    #list =[1, 3, 5, 7, 9]
    #a,b = min([ [i,val] for i, val in enumerate(list)],key = lambda x: 1/(x[1]+1))
    resId = id + 1
    # Update resInc with real inc.
    resInc = float(resData[resId][18])
else:
    resInc = float(resData[resId][18])

#else- > find resonant orbit by ID:
Xstate = [float(ii) for ii in resData[resId][9:15] ]
dt = float(resData[resId][15])
print('resId =  {} \nResonant Orb INC (deg) = {:.6}'
      .format(resId+1,float(resInc)))


#==============================================================================
# Functions and Objects set up
#==============================================================================

# Cartesian State
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
# Frames:
#-----------------------------------------------------------------

#Enceladus J2000
j2k_enceladus = rcpy.Frame()
j2k_enceladus.frametype_id    = rcpy.CopFrameEnum.j2000    # 1
j2k_enceladus.framecenter_id  = 1   #main-> 1, aus -> 2, ..
j2k_enceladus.mainbody.id     = 602   #Enceladus
j2k_enceladus.auxbody.id      = 699   # Saturn

#Enceladus body Fixed Frame:
iau_enceladus = rcpy.Frame()
iau_enceladus.frametype_id    = rcpy.CopFrameEnum.iau_body_fixed   # -4
iau_enceladus.framecenter_id  = 1   #main-> 1, aux -> 2, ..
iau_enceladus.mainbody.id     = 602
j2k_enceladus.auxbody.id      = 699   # Saturn


# ----------------------------------------------------------------------------
#Functions:
# ----------------------------------------------------------------------------

def setSegment(segIn,centralBodyId,stateType,stateFrame,stateIn=[]):
    segIn.def_central_body.id=centralBodyId
    segIn.state=deepcopy(stateType)
    segIn.state.frame=stateFrame
    if stateIn:
        for ii in range(6):
            segIn.state[ii].value=stateIn[ii]

#==============================================================================
# Create Ideck Sat-Enc system, and add X0,dt:
#==============================================================================

#-----------------------------------------------------------------
# Create Ideck:
#-----------------------------------------------------------------
newIdeck = rcpy.Ideck()
#newIdeck.filename = outputIdeck

#-----------------------------------------------------------------
# Force Variables:
#-----------------------------------------------------------------
force_bodies_list= [699,602]

newIdeck.force_vars = {
    'force_date_choice': 2,           #Set to Julian Date:
    'force_julian': global_t0_JD,
    'force_frame':  j2k_enceladus,
    'force_time_system_cal': 'TDB',    #optional: just to be sure
    'force_time_system_jd': 'JDTDB',  #optional: just to be sure
    'force_time_system_et': 'TDB',    #optional: just to be sure
    'force_bodies': force_bodies_list,
}

#-----------------------------------------------------------------
#ADD pck.sat427.tpc to SPICE kernels:
#-----------------------------------------------------------------
newIdeck.spice.append('/Users/ricgomez/.copernicus/support_files/gravity/pck.sat427.tpc')

#-----------------------------------------------------------------
# Grapichs:
#-----------------------------------------------------------------
newIdeck.graphics.vis_frame      = iau_enceladus
newIdeck.graphics.bodies_to_plot = [699,602]
newIdeck.graphics={
    'ogl_enableprintingmode': False,   #use Printer friendly mode color
    'ogl_maxbodytimestep': 0.01,       #(days), times step to plot bodies trajecotries
    'grap_axis_length': 1.0E+03,
}

#-----------------------------------------------------------------
#Create Segment and  Define state:
#-----------------------------------------------------------------
# Create Segment:
segName =  'Res_' + str(resId)
newSeg=rcpy.Segment(segName)

# # Create State:
# stateX0 = rcpy.State()
# # set Frame:
# stateX0.frame = iau_enceladus

# set State:
setSegment(newSeg,Enceladus_ID,stateX,iau_enceladus,Xstate)

# Set time:
newSeg.time.t0 = 0.0      #S1.time.t0 = {'ov': False, 'value': 0.0}
newSeg.time.dt = dt/86400.0  #secs
newSeg.time.dt.ov  = True

# Add Enceladus Spherical Harmonics
newSeg.gravfilename = 'Enceladus_Harmonics.csv'
newSeg.central_body_grav_model = 2
newSeg.grav_degree = 3
newSeg.grav_order = 3
newSeg.gravframe = -40

# Segment Grpahics:
newSeg.plot_data.plot_color=copColor.blue  # BLUE

#-----------------------------------------------------------------
# Set Ideck:
#-----------------------------------------------------------------
#Load segments:
seg_list=[newSeg]
newIdeck.segments=seg_list

#-----------------------------------------------------------------
# Save Ideck:
#-----------------------------------------------------------------

# #Save new files with resId:
# outputIdeck = baseName + '_ID' + str(resId +1).zfill(3) + '.ideck'

#Save new files with resInc:
outputIdeck = baseName + '_Inc' + str(round(resInc)) + 'Deg.ideck'


newIdeck.save(outputIdeck)
print('{} Saved'.format(outputIdeck))
# ===============================================
