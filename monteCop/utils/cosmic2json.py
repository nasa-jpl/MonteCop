#!/usr/bin/env mpython
#===========================================================================
#
# Section 392 Navigation and Mission Design
#
# Copyright (C) 2018, California Institute of Technology.
# U.S. Government Sponsorship under NASA Contract NAS7-03001 is acknowledged.
#
#===========================================================================

""" Convert a Cosmic input.py file into a .json file.
WARNING: This script does NOT function if incomplete sets of controls are
used when a subset of the controls is fixed. The script changes the state
parameters, rendering Monte's interpretation of fixed controls inadequate.
This can be remedied by using a template that includes correct apoStateParams
and flybyStateParams.
"""
from __future__ import print_function
from builtins import str
from builtins import range

#===========================================================================
# Place all imports after here.
#
from collections import OrderedDict
import argparse
import sys
import re
import json

import Monte as M
from mpy.opt.cosmic import Manager
from mpy.units import s, km, deg, year, hour, day, sec

# Place all imports before here.
#===========================================================================
#Globals:
global frame
global flybyStateParams
global apoStateParams

frame= None

# Regex for finding CosmicEvent start point and delta
# Need 2 different ones because the format of the cosmicEvent.name() changed,
# so regex1 will work with M<140, regex2 with M>144
regex1 = re.compile("^.*('.*').*('.*')")
regex2 = re.compile("^.*('.*').*\s([0-9].*)")

# State Params Types:

#Default:
stateParams = [
   "Cartesian.x",
   "Cartesian.y",
   "Cartesian.z",
   "Cartesian.dx",
   "Cartesian.dy",
   "Cartesian.dz"
]

# cartStateParams = []
#
# flybyStateParams = []
#
# hyperStateParams = []
#
# apoStateParams = []
#
# periStateParams = []

# Read state ParamsSet From configFileJSON, if default
# (e.g. flybyStateParams = 'Default') use the defaults here
# for future, if flybyStateParams = 'Inputs' use the sateParam input, so
# in Monte-Cop, convert Monte params to Copernicus params as they are inputs
# something like this need to be done:
#        ==:
#        Conic.vInfinity  --> Param1Enum.vin
#        Conic.inboundRA  --> Param2Enum.ra_vinf
#        Conic.inboundDec --> Param3Enum.dec_vinf
#           ...
#        Conic.longitudeOfNode  --> Param4Enum.raan
#        Conic.argumentOfPeriapsis --> Param5Enum.aop
#
# If there is some incompatibility or Monte State do not exist on COPernicus
# use Cartesian as the default paramSet

#===========================================================================
# Utils functions:

# Disable print on screen
def _blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore
def _enablePrint():
    sys.stdout = sys.__stdout__

#===========================================================================
def formatControlBounds( controls) :

   bounds = []
   for control in controls:
      bounds.append([str(control.minBound()),
                     str(control.maxBound())])
   return bounds

#===========================================================================
def controlPointJSON( cp ):
   """ Create a dictionary with the control point information. This dictionary
   can be easily translated into JSON format.

   = INPUT VARIABLES
   - cp           controlPoint to format

   = RETURN VALUE
   - An OrderedDict containing the controlPoint information
   """

   global frame
   global flybyStateParams
   global apoStateParams

   cpJson = OrderedDict()
   cpJson[ 'Type' ]   = 'ControlPoint'
   cpJson[ 'Time' ]   = cp.time().format( 'full' )
   cpJson[ 'Name' ]   = cp.name()
   cpJson[ 'Center' ] = cp.state().center()
   cpJson[ 'Mass' ]   = str(cp.mass())

   state = cp.state()
   if frame is not None:
      state = state.relativeTo(cp.state().center(), frame)
   cpJson[ 'Frame' ]  = state.frame()


   #cpJson[ 'OrbitType']    = XXX
   # To difne orbitType, check eccentricity wrt center:
   #  -> inside SOI? use OEH1, else, OEH2
   # for now only use   OEH

   # Define the param sets according to the orbit type: e.g.
   # hyperbolic( usualy flybys): openOrbit -> flybyStateParams
   # Elliptic and circular: closeOrbit -> apoStateParams or stateParams
   #if openOrbit, check if flyby or just hyperbolic-> check if inside S.O.I.

   #ecc=Conic.eccentricity(state)
   #if Conic.eccentricity(state) > 0:
   #    cpJson[ 'StateParams' ] = apoStateParams
   #else:
   #    cpJson[ 'StateParams' ] = flybyStateParams

   #Define frames accordingly

   # For Bounds: if control is not a CP param skip an warn user ,
   # BUT a coulple of eception can be added:
   # tarnsform bounds (EXCEPTINONS),
   # For flybyStateParams the only changes are:
   #    -> timeFromPeriapsis to trueAnomaly
   #    -> periapsisAltitude to periapsisRange

   if any(x in cpJson[ 'Name' ] for x in ['Apo', 'Peri']):
      cpJson[ 'StateParams' ] = apoStateParams
      cpJson[ 'State' ]  = [
         str( makeCoordinateFunc( param )( state ) ) for param in apoStateParams
      ]
      cpJson['StateParamsType'] = 'apoStateParams'
   else:
      cpJson[ 'StateParams' ] = flybyStateParams
      cpJson[ 'State' ]  = [
         str( makeCoordinateFunc( param )( state ) ) for param in flybyStateParams
      ]
      cpJson['StateParamsType'] = 'flybyStateParams'
   cpJson['ControlParams'] = [par.rsplit('/',1)[-1] for par in cp.controls().params().params()]
   cpJson['ControlBounds'] = formatControlBounds(cp.controls())

   return cpJson

#===========================================================================
def breakPointJSON( bp ):
   """ Create a dictionary with the break point information. This dictionary
   can be easily translated into JSON format.

   = INPUT VARIABLES
   - bp           breakPoint to format

   = RETURN VALUE
   - An OrderedDict containing the breakPoint information
   """
   global frame
   global flybyStateParams
   global apoStateParams

   bpJson = OrderedDict()
   bpJson[ 'Type' ]   = 'BreakPoint'
   bpJson[ 'Name' ]   = bp.name()

   # Handle breakpoint mode
   if bp.mode() == M.BreakPoint.Mode.MATCH_ELEM:
      bpJson[ 'Mode' ]   = 'MATCH_ELEM'
   elif bp.mode() == M.BreakPoint.Mode.MATCH_MAG:
      bpJson[ 'Mode' ]   = 'MATCH_MAG'
   elif bp.mode() == M.BreakPoint.Mode.FREE_DV_MAG:
      bpJson[ 'Mode' ]   = 'FREE_DV_MAG'
   elif bp.mode() == M.BreakPoint.Mode.FREE_DV_ELEM:
      bpJson[ 'Mode' ]   = 'FREE_DV_ELEM'
   else:
      print('MODE NOT RECOGNIZED!! ERROR')

   # Handle time
   if bp.isEventBased():
      eventName = bp.event().name()
      res = regex.match(eventName)
      bpJson[ 'Start' ] = OrderedDict()
      bpJson[ 'Start' ][ 'Name' ] = res.group(1)
      bpJson[ 'Start' ][ 'Delta' ] = res.group(2)
   else:
      bpJson[ 'Time' ] = bp.time().format( 'full' )

   if 'Unassigned' not in bp.center():
      bpJson[ 'Center' ] = bp.center()

   bpJson[ 'Frame' ]  = bp.frame()
   bpJson[ 'PosTol']  = str(bp.posTol())
   bpJson[ 'VelTol']  = str(bp.velTol())
   bpJson[ 'MassTol'] = str(bp.massTol())

   # Get the states and difference
   sMinus = M.State()
   sPlus  = M.State()
   bp.state(sMinus, sPlus)
   bpJson[ 'State'] = OrderedDict()
   bpJson[ 'State' ][ 'Minus' ] = [str(x *km) for x in sMinus.pos()] + \
                                  [str(x *km/s) for x in sMinus.vel()]
   bpJson[ 'State' ][ 'Plus' ]  = [str(x *km) for x in sPlus.pos()] + \
                                  [str(x *km/s) for x in sPlus.vel()]
   bpJson[ 'State' ][ 'Diff' ]  = [str(x *km) for x in sPlus.pos() - sMinus.pos()] + \
                                  [str(x *km/s) for x in sPlus.vel() - sMinus.vel()]

   return bpJson

#===========================================================================
def makeCoordinateFunc( funcStr ):
   """ Return the Monte method (function) given its name. This is used for state
   coordinates. Examples: Conic.inclination, Cartesian.x.

   = INPUT VARIABLES
   - funcStr      name of the module and method

   = RETURN VALUE
   - state coordinate method
   """
   mod, method = funcStr.split( '.' )
   return getattr( getattr( M, mod ), method )

#===========================================================================
def impulseBurnJSON( burn, burn2 ):
   """ Create a dictionary with the impulse burn information. This dictionary
   can be easily translated into JSON format.

   = INPUT VARIABLES
   - burn         burn to format (from ImpulseBurnMgr)
   - burn2        burn to format (from OptCosmic. Represents the same burn,
                  but has different info stored)

   = RETURN VALUE
   - An OrderedDict containing the breakPoint information
   """
   # Check that we have the same burn twice
   if burn.name() not in burn2.name():
      print('ERROR')
      print(burn)
      print(burn2)
      exit()

   dvJson = OrderedDict()
   dvJson[ 'Type' ] = 'ImpulseBurn'
   dvJson[ 'Name' ] = burn.name()

   # # Handle time
   dvJson[ 'Time' ] = burn.time().format( 'full' )
   dvJson[ 'isEventBased' ] = str(burn.isEventBased())
   if burn.isEventBased():
      eventName = burn.event().name()
      res = regex1.match(eventName)
      if not res:
         res = regex2.match(eventName)
      dvJson[ 'Start' ] = OrderedDict()
      dvJson[ 'Start' ][ 'Name' ] = res.group(1).strip('\'')
      dvJson[ 'Start' ][ 'Delta' ] = res.group(2).strip('\'')

   dvJson[ 'Frame' ]    = burn.frame()
   dvJson[ 'DeltaVel' ] = [ str( dv_i * km/sec ) for dv_i in burn.dvel() ]
   dvJson['ControlParams'] = [p.rsplit('/',1)[-1] for p in burn2.controls().params().params()]
   dvJson['ControlBounds'] = formatControlBounds(burn2.controls())

   return dvJson

# ======================================================================
# Main
# ======================================================================

def cosmic2json(inputFile,jsonTemplate,saveToFile = False,jsonFileOut=None):
    """ Main function of cosmic2json module

    Convert a Cosmic input.py file into a json solution. json2cosmic use a json template with user
    definitions and default params (like frames, ParamsSet, Copernicus setup)

    = INPUT VARIABLES
    - inputFile         cosmic input file
    - jsonTemplate      json template with predefined setup information
    - saveToFile        save json solution to file (default: False)
    - jsonFileOut       name of json file to be generated. if None, cosmic name is used

    = RETURN VALUE
    - solJson           json solution
    - mgr               cosmic Manager
    """

    global frame
    global flybyStateParams
    global apoStateParams

    with open(jsonTemplate, 'r') as temp:
       template = json.load(temp, object_pairs_hook=OrderedDict)

    # Create OrderedDict. If template file inputted, start filling the dict
    solJson = OrderedDict()
    # Start filling the dict
    for key in template:
      solJson[key] = template[key]
    # If template contains stateParams, replace default ones
    if 'StateParams' in template['Defaults']:
      stateParams = template['Defaults']['StateParams']
    if 'apoStateParams' in template['Defaults']:
      apoStateParams = template['Defaults']['apoStateParams']
    else:
      apoStateParams = stateParams
    if 'flybyStateParams' in template['Defaults']:
      flybyStateParams = template['Defaults']['flybyStateParams']
    else:
      flybyStateParams = stateParams

    # If the breakpoint center is not defined in cosmic file
    if 'BPcenter' in template['Defaults']:
      BPcenter = template['Defaults']['BPcenter']

    else:
       print('No template file inputted. The timeline only will be written to output')

    #rlr:
    #frame = None
    #Others
    if 'Frame' in template['Defaults']:
        frame = template['Defaults']['Frame']


    #===========================================================================
    # Start Cosmic:
    #===========================================================================
    #raise Exception('exit')
    print('Generating trajectory ... ')
    boa=M.BoaLoad()
    mgr=Manager(boa)
    mgr.loadInput(inputFile)
    mgr.tl.createTraj(mgr.boa, mgr.problem, True, False, False)
    print('Trajectory created. ')


#    raise

    # Read and add Control Points
    controlPoints = []
    for i in range(mgr.cosmic.timeline().numControl()):
       cp = mgr.cosmic.timeline().controlPoint(i)
       controlPoints.append(controlPointJSON(cp))
    solJson['ControlPoints'] = controlPoints

    # Read and add Breakpoints
    breakPoints = []
    for i in range(mgr.cosmic.timeline().numBreak()):
       bp = mgr.cosmic.timeline().breakPoint(i)
       breakPoints.append(breakPointJSON(bp))
    solJson['breakPoints'] = breakPoints

    # Read and add Burns
    burnList  = M.ImpulseBurnMgrBoa.read(mgr.boa, mgr.cosmic.timeline().body())
    burnList2 = mgr.cosmic.burns()
    burns = []
    for i in range(len(burnList)):
       burns.append(impulseBurnJSON(burnList[i], burnList2.find(burnList[i].name())))
    solJson['ImpulseBurns'] = burns

    # Write to json file
    if saveToFile:
        if not jsonFileOut: jsonFileOut=inputFile.replace('.py','.json')
        with open( jsonFileOut, 'w+' ) as outfile:
           json.dump( solJson, outfile, indent = 4, separators=(',', ': ') )

    return solJson, mgr
