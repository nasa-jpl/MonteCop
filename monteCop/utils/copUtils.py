# ==============================================================================
#  Module MonteCop Utils (mcpUtils):
#
#   Decription:  Utils for Monte-copernicus interface...
#
#   last upfate: July 24 2020.
#   By Ricardo L. Restrepo
#
# ==============================================================================

# ==============================================================================
# imports here:

import sys

# import robocoppy as rcpy
import monteCop.src.CopPy510.robocoppy as rcpy
from monteCop.src.CopPy510.robocoppy import ColorEnum as copColor
from monteCop.utils.spiceIDs import *

# from MonteCop_Dev.mcp.CopPy.robocoppy.frames import *
# from MonteCop_Dev.mcp.CopPy.robocoppy.serialization import S11n


#For debugging:
from importlib import reload

# from copy import deepcopy

# ==============================================================================


# ==============================================================================
# Handle Monte and SPICE IDs
# ==============================================================================

def passBodyID(bodyID):
    ''' Check and pass SPICE BodyID (accept str or num)
        = Inputs:
        - bodyID    Body name (e.g. "Earth") or ID (e.g., 399)
        = Outputs:
        - bosyID    SPICE body ID
    '''
    # NOTE: body name is upper-case on SpiceBodyName
    if type(bodyID) == str and bodyID.upper() in SpiceBodyName.keys():
        return SpiceBodyName[bodyID.upper()]
    elif bodyID in SpiceBodyName.values():
        return bodyID
    else:
        #sys.exit('Invalid SPICE Body Name or ID: '+ str(bodyID))
        print('Invalid SPICE Body Name or ID: '+ str(bodyID))

def passFrameID(frameName):
    ''' Convert Monte frame name into Copernicus frame ID (SPICE)
        = Inputs:
        - frameName  Monte frame name
        = Outputs:
        - frameID   Copernicus frame ID
    '''
    if frameName in m2cFrameID.keys():
        return m2cFrameID[frameName]
    else:
        #sys.exit('Invalid Frame Name: '+ str(frameName))
        print('Invalid Frame Name: '+ str(frameName))

# ==============================================================================
# Plots Utils
# ==============================================================================

#---------------------------------------------
# Lists:
#---------------------------------------------
#Color list:
colorListA={
        '0' : copColor.seagreen,
        '1' : copColor.orange,
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

# A subset of Copernicus color list (CopPY.colors) which have color and dark-colors
# has been selected as defaults. dark-color are going to be used for backward prop. seg_list
# later, for user defined colors, generate the dark using the CopColor Class and mod rbg
colorList={
        '0' : copColor.blue,
        '1' : copColor.lightcyan,
        '2' : copColor.green,
        '3' : copColor.fuchsia,
        '4' : copColor.red,
        '5' : copColor.violet,
        '6' : copColor.orange,
        '7' : copColor.salmon,
        '8' : copColor.seagreen ,
        '9': copColor.gray
        }

colorDarkList={
        '0' : copColor.darkblue,
        '1' : copColor.darkcyan,
        '2' : copColor.darkgreen,
        '3' : copColor.darkmagenta,
        '4' : copColor.darkred,
        '5' : copColor.darkviolet,
        '6' : copColor.darkorange,
        '7' : copColor.darksalmon,
        '8' : copColor.darkseagreen,
        '9': copColor.darkgray
        }

# example:
# cl['1'].value
# cl['1'].name

# genColorLight()
# chekc if light color exist, otherwise, substract 10 (-10) to each color comp-
# component, e.g. r,b,g   (if colorNumber > than 10, otherwise c = int(colorNumber/2) )


# ==============================================================================
# Cop Lists:
# ==============================================================================

# #Get color Names
# colorNames=[]
# for cn in colorList:
#      colorNames.append(str(colorList[str(cn)]).replace('ColorEnum.',''))

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

#ASSUME values: (value(1) inherit(2), DValue(3))

#myColor=rcpy.colors.CopColor([128, 128, 128])

# ==============================================================================
# State Params:
# ==============================================================================

#m2c_paramSet = []

#Defualt state params:R,V. so, not need to be defined here, but added for completeness:
stateX        = rcpy.State()
stateX.param.params_id = [
            rcpy.Param1Enum.rx,
            rcpy.Param2Enum.ry,
            rcpy.Param3Enum.rz,
            rcpy.Param4Enum.vx,
            rcpy.Param5Enum.vy,
            rcpy.Param6Enum.vz,
            ]

#------> State OE (Classic Orbtila Elements v1):
stateOE1        = rcpy.State()
stateOE1.param.params_id = [
            rcpy.Param1Enum.sma,
            rcpy.Param2Enum.ecc,
            rcpy.Param3Enum.inc,
            rcpy.Param4Enum.raan,
            rcpy.Param5Enum.aop,
            rcpy.Param6Enum.ta,
            ]
stateOE1.param.angle_unit = rcpy.AngleUnits.deg  # (need ot be defined)

#------> State OE (Classic Orbtila Elements):
stateOE2        = rcpy.State()
stateOE2.param.params_id = [
            rcpy.Param1Enum.period,
            rcpy.Param2Enum.rp,
            rcpy.Param3Enum.inc,
            rcpy.Param4Enum.raan,
            rcpy.Param5Enum.aop,
            rcpy.Param6Enum.ta,
            ]
stateOE2.param.angle_unit = rcpy.AngleUnits.deg  # (need ot be defined)

#------> State OEH (Hyperbolic Orbtila Elements):
stateOEH        = rcpy.State()
stateOEH.param.params_id = [
            rcpy.Param1Enum.vin,        #vInfinity
            rcpy.Param2Enum.ra_vinf,     #inboundRA
            rcpy.Param3Enum.dec_vinf,   #inboundDec
            rcpy.Param4Enum.btheta,     #bPlaneTheta
            rcpy.Param5Enum.rp,         #periapsisRange
            rcpy.Param6Enum.ta,         #TRUE-ANOMALY (careful, NOT timeFromPeriapsis)
            ]
stateOEH.param.angle_unit = rcpy.AngleUnits.deg  # (need ot be defined)

# ==============================================================================
# Frames:
# ==============================================================================
j2k = rcpy.Frame()
j2k.frametype_id = rcpy.CopFrameEnum.j2000  # 1           #-> SPICE: EME2000
j2k.framecenter_id = 1  # main-> 1, aux -> 2, barycenter -> 3, L1 ->4, ...
# j2k.mainbody.id     = 399
# j2k.auxbody.id      = 301

eclipj2k = rcpy.Frame()
eclipj2k.frametype_id = rcpy.CopFrameEnum.eclipj2000  # 17     #-> SPICE: EM02000
eclipj2k.framecenter_id = 1  # main-> 1, aus -> 2, ..
# eclipj2k.mainbody.id     = 399
# eclipj2k.auxbody.id      = 301

# iau_body_fixed:
iau_bodyFixed = rcpy.Frame()
iau_bodyFixed.frametype_id = rcpy.CopFrameEnum.iau_body_fixed  # -4
iau_bodyFixed.framecenter_id = 1  # main-> 1, aux -> 2, ..
# iau_bodyFixed.mainbody.id     = 399
# iau_bodyFixed.auxbody.id      = 301

# iau_europa:
iau_europa = rcpy.Frame()
iau_europa.frametype_id = rcpy.CopFrameEnum.iau_body_fixed  # -4
iau_europa.framecenter_id = 1  # main-> 1, aux -> 2, ..
iau_europa.mainbody.id = 502
# or: use IauEuropaFixed

# create Copernicus frame from monte name:
def m2cFrame(monteFrame, mainBodyId="Earth", auxBody="Moon"):
    newFrame.frametype_id = passFrameID(monteFrame)
    newFrame.mainbody.id  = passBodyID(mainBodyId)
    newFrame.auxbody.id   = passBodyID(auxBody)
    newFrame.framecenter_id = 1
    return newFrame


# #Not working:(
# class m2cFrameC(S11n):
#     def __init__(self,monteFrame, mainBodyId="Earth", auxBody="Moon"):
#         self.frametype_id = rcpy.CopFrameEnum.j2000
#         self.framecenter_id  = rcpy.FrameCenter.main    #1
#         #self.mainbody    = rcpy.Body(numSpiceBodyID(mainBodyId))
#         #self.auxbody     = rcpy.Body(numSpiceBodyID(auxBody))
#         self.mainbody    = rcpy.Body(rcpy.SpiceBodyEnum(399))
#         self.auxbody     = rcpy.Body(rcpy.SpiceBodyEnum(301))
#         self.scale = 100000.0
#         self.et = 0.0
#         self.zaxis_id = rcpy.ZAxisEnum.ang_mom
#         super().__init__(*args, **kwargs)

# ==============================================================================

# ==============================================================================
# Util functions:
# ==============================================================================

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






# ==============================================================================
#  ...
# ==============================================================================

# Explicittly definte Params set:

# flybyStateParams:

# hypeStateParams: (not flyby but hyperbolic state)

# stateParams

# apoStateParams

# periStateParams


# ==============================================================================
# XXX:
# ==============================================================================


# List of features to add:

# - getLightColorFrom()
# - getDarkColorFrom()
# - Add SPICE ID's dictionary
