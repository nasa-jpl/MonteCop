#===========================================================================
#
# Section 392 Navigation and Mission Design
#
# Copyright (C) 2018, California Institute of Technology.
# U.S. Government Sponsorship under NASA Contract NAS7-03001 is acknowledged.
#
#===========================================================================

#===========================================================================
# cosmicBatch()  template
#
#   NOTE: To be used with cosmicBatch() v1.0
#
# __version__ = '1.0'
# __author__  = 'Ricardo L. Restrepo (392M)'
#===========================================================================
import Monte as M
import datetime
import os.path
import sys

import atnLib.utilities.customController as customController
from atnLib.utilities.translatScaler import translatScaler
import rlrLibs.cosmicPlus.cosmicPlus as cosmicPlus

#set Problem name:
inputFileName = __file__
problemName=os.path.basename(inputFileName)[:-3]

#Gloabal Params:
scName = 'encnf5'
startTime = datetime.datetime.now()
stastOut='stat.out'

#--------------------------------------------------------------------------
# Setup problem tolerances:
#--------------------------------------------------------------------------
PosTolIn     = 1.0E-1
VelTolIn     = 1.0E-4
MajorFeasTol = 1e-5
MajorOptTol  = 1e-3

# PosTolIn     = 1.0E-3
# VelTolIn     = 1.0E-6
# MajorFeasTol = 1e-8
# MajorOptTol  = 1e-3

#--------------------------------------------------------------------------
# Setup SNOPT params for cosmicBatch() run:
#--------------------------------------------------------------------------
MaxIterLimit   = 300                 #Major iterations limit

#set SNOPT params (these values are overwritten on cosmicBatch.py calls):
MajorStepLimit = 0.1
LinesearchTol  = 0.99

paramList={
        'MajStL': MajorStepLimit,
        'LinSeT': LinesearchTol,
        'MaxIte': MaxIterLimit,
        'MajFeT': MajorFeasTol,
        'MinOpT': MajorOptTol,
        'PosTol': PosTolIn,
        'VelTol': VelTolIn,
        }

runInBatch=False                   #DO NOT MODIFY (overwritten by cosmicBatch.py calls)
if runInBatch:
    baseName='_'.join(problemName.split('_')[:-2])
    stastOut='stat_'+baseName+'.out'

snoptOut = -1
for i in range(7,100,2):
   if not os.path.isfile('fort.{0}'.format(i)):
      snoptOut = i
      break
# Overwrite over the default...
if snoptOut < 0:
   snoptOut = 30

if runInBatch:
    snoptOut=-1

opt = NewSnOpt(
   Name          = 'SNOPT',
   MaxIterations = MaxIterLimit,
   Options       = {
      'Major iterations limit'      : MaxIterLimit,
      'Major feasibility tolerance' : MajorFeasTol,
      'Major Optimality Tolerance'  : MajorOptTol,
      'Scale Option'                : 0,
      #'Print file'                  : snoptOut,
      #'Summary file'                : snoptOut + 1,
      'Minor Print Level'           : 0,
      'Major Print Level'           : 10100,
      'Print frequency'             : 1,
      'Summary frequency'           : 100,
      'Major Step Limit'            : MajorStepLimit,
      'Hessian Flush'               : 999999,
      'Linesearch tolerance'        : LinesearchTol,
      'Verify level'                : 0,
      'System information yes'      : '',
      'Scale Print'                 : '',
      'Function Precision'          : 1e-10,
   }
)

#--------------------------------------------------------------------------
# <-- SNOPT Param setup|
#--------------------------------------------------------------------------

DefaultData(
   Inputs = [
      # Default Data - the initial default are not used since lockfile
      # contains necessary data
      "time",
      "body",
      "frame/IAU 2000",
      "frame/inertial",
      # Ephemeris
      "/nav/common/import/ephem/de430.boa",
      "/nav/common/import/ephem/sat375l.boa",
      #"/home/bdanders/encnf5/ephem/sat427.bsp",
      # Lockfile. Always load last
      #"/home/bdanders/encnf5/lockfile/lockfile_0.6.0.boa",
      #"/home/bdanders/encnf5/lockfile/lockfile_0.6.1.boa",
      "/group/encnf5/lockfile/lockfile_0.6.2.boa",
   ]
)

Defaults(
   # In order to reuse lockfile force model and integrator, make sure Body
   # is 'Europa Clipper'
   Body  = scName,
   Frame = 'EMO2000',
   Mass  = 4500.000 *kg,
   Isp   = 315.0 *sec,

   Propagator = 'DIVA',
   PosTol     = PosTolIn *km,
   VelTol     = VelTolIn *km/sec,
   MassTol    = 0 *kg,
)

Unit.setFormat('deg')
Epoch.setFormat('full')

EditErrorControl(
   Entries = {
      ErrorId.FORCE_MODEL_DUPLICATE : 'IGNORE',
      ErrorId.COSMIC_DUPLICATE_BURN: 'IGNORE',
      ErrorId.GRAV_BELOW_SURFACE: 'IGNORE', # ignore warnings of sub-surface integration, which happens during optimization of very low orbits
      },
   )


# # Turn off below surface warning
# errControl = M.ErrorControl( boa )
# errControl.setAction( M.ErrorId.GRAV_BELOW_SURFACE, M.ErrorAction.IGNORE )

# NewDivaPropagator(
#     Name     = "DIVA", # Has to be the same name as the one used in defaults above
#     StateTol = 1e-10,
#     MassTol  = 1e-6,
#     MinStep  = 1e-6 *sec,
#     MaxStep  = 100 *day,
#     Forces   = [
#        "Gravity",
#        "Impulse Burn",
#        #"Finite Burn"
#        ]
#     )

# Add forces
EditIntegState(
   AddForces = [
      "Gravity",
      "Impulse Burn",
      "Finite Burn"
   ],
   DelForces = [
      "Solar Pressure",
   ],
)

#===========================================================================
# Constraints inputs
#===========================================================================
# Create dictionary of user defined tolerances for translatScaler
userTol = {}

#===========================================================================
# Optimizer inputs
#===========================================================================
NewDblseOpt(
   Name             = 'DBLSE',
   MaxIterations    = 100,
   DblseTol         = 1.0E-9,
   CostRelTol       = 1.0E-2,
   CostAbsTol       = 1.0E-2,
   ConstraintRelTol = 1.0E-4,
)

NewIpOpt(
   Name          = "IPOPT",
   MaxIterations = 200,
   ConstraintTol = 1e-4,
   NlpTol        = 1e-4,
)

# Add custom Scaler
myScaler = PyOptScaler(translatScaler(userTol, constraintScaleFactor=1.0E-5))
#myScaler = PyOptScaler(translatScaler(userTol, constraintScaleFactor=MajorFeasTol))

#--- Create Cosmic
NewCosmic(
   # Name         = problemName,
   MassModel    = False,
   Optimizer    = 'SNOPT',
   Cost         = 'MIN DV',
   # ScalingUnits = [ 'LU', 'TU', 'MU' ],
   Scalers      = [
      myScaler,
   ],
   # Scalers      =[
   #     NewOptInitValueScaler(
   #        ControlScaling = True,
   #        ConstraintScaling = True,
   #        CostScaling = False,
   #    ),
   # ],
   CostFrame       = "EMO2000",
   BreakPointFrame = "EMO2000",
)

#--- Create Cosmic Output
NewCosmicOutput(
   TrajFile = '',
   Elements = {
      'Traj Info' : cosmicPlus.outputCostAndCons,
      #'largeMvrBrief' : output.largeMvrBrief,
   },
   Groups = {
      'Sol Begin' : [ 'tlBrief' ],
      'Iter Init' : [ 'iterNum' ],
      'Iter End' : [ 'Traj Info' ],
   },
)

#--- Add custom controller
customController = M.PyOptController(
   customController.CustomController(
      M.OptCosmicBoa.read(boa, M.OptCosmicBoa.getAll(boa)[0]),
      problemName   = problemName.split('_CV_')[0],
      saveChkPt     = 100,
      chkPtTemplate = inputFileName,
      startTime     = startTime
   )
)

DebugFlag.setLevel('OptProblem', 1)

#===========================================================================
# Frames:
#===========================================================================
# NewInertialFrame(
#    Frame = 'SE_ROT_INERTIAL',
#    RefFrame = 'SE_ROT_FRAME',
# )
#
# NewBodyPosDirFrame(
#    Frame = 'SE_ROT_FRAME',
#    BaseFrame = CoordName.EME2000,
#    Body = 'Enceladus',
#    Center = 'Saturn',
#    )
#
# NewBodyVelDirFrame(
#    Frame = 'VUW',
#    BaseFrame = CoordName.EME2000,
#    Body = 'encnf5',
#    Center = 'Sun',
#    )
#
# NewBodyPosDirFrame(
#    Frame = 'VUW_JUPITER',
#    BaseFrame = CoordName.EME2000,
#    Body = 'encnf5',
#    Center = 'Jupiter',
#    )
#
# NewBodyPosDirFrame(
#    Frame     = 'SatEnc Rotating Frame',
#    BaseFrame = 'EME2000',
#    Body      = 'Enceladus',
#    Center    = 'Saturn Barycenter',
#    )

# Build mean equator and equinox frames for the Jupiter bodies
#
NewMeanEquatorEquinoxFrame(
   Frame = "Saturn Mean Equator And Equinox of Epoch",
   Body = "Saturn",
   )


fb = [ "Mimas", "Enceladus", "Tethys", "Dione", "Rhea", "Titan", "Iapetus" ]
for body in fb:
   NewMeanEquatorEquinoxFrame(
      Frame = "%s Mean Equator And Equinox of Epoch" % body,
      PoleBody = body,
      OrbitBody = "Saturn",
      )

# add GM from Ryan Park's recommended file (pck.sat427.tpc)
#EditGm(
#   Body = BodyName.Enceladus,
#   Gm = 7.210497553340731 *km**3/sec**2
#)

# M.BodyVelDirFrame(boa, 'DVFrame', 'EME2000', M.TimeInterval(), sc, secondary)

#-----------------------------------------------------------------------------
# Add gravity models:

#basicGrav.add( boa, sc, [primary,secondary])
#basicGrav.add( boa, sc, [primary,secondary], [], {'Enceladus':'Enceladus'})
#basicGrav.add( boa, sc, [primary,secondary], [], {'Enceladus':'Enceladus365'})
#basicGrav.add( boa, sc, [primary,secondary], [], {'Saturn':'Saturn'})
#basicGrav.add( boa, sc, [primary,secondary], [], {'Saturn':'Saturn365'})
#basicGrav.add( boa, sc, [primary,secondary], [], {'Enceladus':'Enceladus','Saturn':'Saturn'})
#basicGrav.add( boa, sc, [primary,secondary], [], {'Enceladus':'Enceladus365','Saturn':'Saturn365'})
#basicGrav.add( boa, sc, [primary,secondary,'Titan','Jupiter Barycenter','Mercury','Venus','Earth Barycenter','Mars Barycenter','Neptune Barycenter','Uranus Barycenter'], [], {'Enceladus':'Enceladus365','Saturn':'Saturn365'})

#===========================================================================
# Cosmic Timeline
#===========================================================================

Timeline = []

# ======================================================================
# Run Cosmic
# ======================================================================
# These commands will run automatically in non-interactive mode
def runCosmic():
   cosmicPlus.runCosmic(Cosmic, problemName,
                  unitNum=snoptOut,
                  startTime=startTime,
                  saveOutput=True,
                  globalStatsFile=stastOut,
                  runTimes=2,  #if max iterations, or minLocal/noConvergence, run again
                  paramsIn=paramList,
                  iterInit=True,
                  )
