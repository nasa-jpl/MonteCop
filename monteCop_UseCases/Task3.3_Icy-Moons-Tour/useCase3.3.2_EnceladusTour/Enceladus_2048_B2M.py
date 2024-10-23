# ===========================================================================
# Section 392 Navigation and Mission Design
#
# Copyright (C) 2021, California Institute of Technology.
# U.S. Government Sponsorship under NASA Contract NAS7-03001 is acknowledged.
#
# ===========================================================================


#===========================================================================
# cosmicBatch()  template
#
#   NOTE: To be used with cosmicBatch() v1.0
#
# __version__ = '1.0'
# __author__  = 'Ricardo L. Restrepo (392M)'
#===========================================================================

import datetime
import os.path
import sys

import Monte as M

import atnLib.utilities.customController as customController
from atnLib.utilities.cosmicConstraints import LS_constraint, Shadow_constraint
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
      #"/home/bdanders/encnf5/ephem/sat427.2040_2050.boa",
      # EOP File
      # "$euclip_eop/simulated.eop",
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
NewInertialFrame(
   Frame = 'SE_ROT_INERTIAL',
   RefFrame = 'SE_ROT_FRAME',
)

NewBodyPosDirFrame(
   Frame = 'SE_ROT_FRAME',
   BaseFrame = CoordName.EME2000,
   Body = 'Enceladus',
   Center = 'Saturn',
   )

NewBodyVelDirFrame(
   Frame = 'VUW',
   BaseFrame = CoordName.EME2000,
   Body = 'encnf5',
   Center = 'Sun',
   )

NewBodyPosDirFrame(
   Frame = 'VUW_JUPITER',
   BaseFrame = CoordName.EME2000,
   Body = 'encnf5',
   Center = 'Jupiter',
   )

NewBodyPosDirFrame(
   Frame     = 'SatEnc Rotating Frame',
   BaseFrame = 'EME2000',
   Body      = 'Enceladus',
   Center    = 'Saturn Barycenter',
   )

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

#===========================================================================
# Cosmic Timeline
#===========================================================================

Timeline = [
   ControlPoint(
      Name = 'ETF-01',
      Active = True,
      Propagator = 'DIVA',
      Time = '21-JUN-2047 00:00:00.000000000 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  4.854651306153028e+01 *km ),
         Conic.bPlaneTheta(  1.635981695717843e+02 *deg ),
         Conic.vInfinity(  7.674743464981433e-01 *km/sec ),
         Conic.inboundDec( -2.309878130320300e+01 *deg ),
         Conic.inboundRA(  4.261888756538856e+01 *deg ),
         Conic.trueAnomaly( -1.078082752008818e-14 *deg ),
         ],
      Controls = [
         [ -8.640000000000000e+04 *sec, 'Cosmic/Cosmic/ETF-01/TIME',  8.640000000000000e+04 *sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-01/Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.000000000000000e+01 *km, 'Cosmic/Cosmic/ETF-01/Conic.periapsisAltitude',  2.000000000000000e+02 *km ],
         [  2.000000000000000e-01 *km/sec, 'Cosmic/Cosmic/ETF-01/Conic.vInfinity',  1.000000000000000e+00 *km/sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-01/Conic.inboundDec',  3.600000000000000e+02 *deg ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-01/Conic.inboundRA',  3.600000000000000e+02 *deg ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo01',
      Start = '22-JUN-2047 06:59:54.649700000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri01',
      Start = '22-JUN-2047 13:13:27.834000000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo02',
      Start = '29-JUN-2047 05:23:54.483800000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri02',
      Start = '30-JUN-2047 14:54:19.915000000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   BreakPoint(
      Name = 'ETF-01->ETF-02',
      Mode = 'MATCH_ELEM',
      Time = '04-JUL-2047 17:12:57.888750000 ET',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo03',
      Start = '05-JUL-2047 22:38:55.677800000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri03',
      Start = '08-JUL-2047 16:34:37.085300000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo04',
      Start = '12-JUL-2047 22:19:16.094200000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DV00',
      Start = '16-JUL-2047 17:53:20.000000000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  3.948379243916926e-03 *km/sec,
      # DeltaVelRA  =  2.291197394791373e+02 *deg,
      # DeltaVelDec =  2.440882024771363e+01 *deg,
      DeltaVel = [
         -2.353166710116206e-03 *km/sec,
         -2.718461859387844e-03 *km/sec,
          1.631646471282444e-03 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -4.000000000000000e-02 *km/sec, 'DX',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DY',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DZ',  4.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri04',
      Start = '16-JUL-2047 18:14:37.302400000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'ETF-02',
      Active = True,
      Propagator = 'DIVA',
      Time = '18-JUL-2047 10:25:55.777500000 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  4.923069984183336e+01 *km ),
         Conic.bPlaneTheta(  1.648081581238715e+02 *deg ),
         Conic.vInfinity(  7.662180298343662e-01 *km/sec ),
         Conic.inboundDec( -2.402631082030650e+01 *deg ),
         Conic.inboundRA(  4.619585979260837e+01 *deg ),
         Conic.trueAnomaly(  7.091252994054027e-06 *deg ),
         ],
      Controls = [
         [ -8.640000000000000e+04 *sec, 'Cosmic/Cosmic/ETF-02/TIME',  8.640000000000000e+04 *sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-02/Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.000000000000000e+01 *km, 'Cosmic/Cosmic/ETF-02/Conic.periapsisAltitude',  2.000000000000000e+02 *km ],
         [  2.000000000000000e-01 *km/sec, 'Cosmic/Cosmic/ETF-02/Conic.vInfinity',  1.000000000000000e+00 *km/sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-02/Conic.inboundDec',  3.600000000000000e+02 *deg ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-02/Conic.inboundRA',  3.600000000000000e+02 *deg ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo05',
      Start = '21-JUL-2047 03:59:36.799600000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   BreakPoint(
      Name = 'ETF-02->ETF-03',
      Mode = 'MATCH_ELEM',
      Time = '23-JUL-2047 03:38:13.496250000 ET',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri05',
      Start = '26-JUL-2047 08:29:29.579000000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DV01',
      Start = '27-JUL-2047 03:28:38.000000000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  8.515429146424126e-03 *km/sec,
      # DeltaVelRA  =  2.459439456223511e+02 *deg,
      # DeltaVelDec =  2.671002934452197e+01 *deg,
      DeltaVel = [
         -3.100749543824222e-03 *km/sec,
         -6.946100755655493e-03 *km/sec,
          3.827475683335635e-03 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -4.000000000000000e-02 *km/sec, 'DX',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DY',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DZ',  4.000000000000000e-02 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'ETF-03',
      Active = True,
      Propagator = 'DIVA',
      Time = '27-JUL-2047 20:50:31.215000000 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  4.956027765583785e+01 *km ),
         Conic.bPlaneTheta(  4.006255105388450e-01 *deg ),
         Conic.vInfinity(  7.384021849065285e-01 *km/sec ),
         Conic.inboundDec( -2.803347304979622e+01 *deg ),
         Conic.inboundRA(  8.066218069019823e+01 *deg ),
         Conic.trueAnomaly( -4.244185921666981e-06 *deg ),
         ],
      Controls = [
         [ -8.640000000000000e+04 *sec, 'Cosmic/Cosmic/ETF-03/TIME',  8.640000000000000e+04 *sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-03/Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.000000000000000e+01 *km, 'Cosmic/Cosmic/ETF-03/Conic.periapsisAltitude',  2.000000000000000e+02 *km ],
         [  2.000000000000000e-01 *km/sec, 'Cosmic/Cosmic/ETF-03/Conic.vInfinity',  1.000000000000000e+00 *km/sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-03/Conic.inboundDec',  3.600000000000000e+02 *deg ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-03/Conic.inboundRA',  3.600000000000000e+02 *deg ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo06',
      Start = '29-JUL-2047 07:16:57.720400000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri06',
      Start = '04-AUG-2047 21:14:34.725300000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo07',
      Start = '05-AUG-2047 04:07:55.998300000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   BreakPoint(
      Name = 'ETF-03->ETF-04',
      Mode = 'MATCH_ELEM',
      Time = '07-AUG-2047 05:51:32.929100000 ET',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo08',
      Start = '12-AUG-2047 01:46:53.666300000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri07',
      Start = '12-AUG-2047 19:16:12.235600000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'ETF-04',
      Active = True,
      Propagator = 'DIVA',
      Time = '17-AUG-2047 14:52:34.643200000 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  4.960760777243501e+01 *km ),
         Conic.bPlaneTheta(  1.689951897374885e+02 *deg ),
         Conic.vInfinity(  7.047582027090139e-01 *km/sec ),
         Conic.inboundDec( -2.594364180477014e+01 *deg ),
         Conic.inboundRA(  5.539972514144985e+01 *deg ),
         Conic.trueAnomaly( -2.532841591086933e-06 *deg ),
         ],
      Controls = [
         [ -8.640000000000000e+04 *sec, 'Cosmic/Cosmic/ETF-04/TIME',  8.640000000000000e+04 *sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-04/Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.000000000000000e+01 *km, 'Cosmic/Cosmic/ETF-04/Conic.periapsisAltitude',  2.000000000000000e+02 *km ],
         [  2.000000000000000e-01 *km/sec, 'Cosmic/Cosmic/ETF-04/Conic.vInfinity',  1.000000000000000e+00 *km/sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-04/Conic.inboundDec',  3.600000000000000e+02 *deg ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-04/Conic.inboundRA',  3.600000000000000e+02 *deg ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo09',
      Start = '20-AUG-2047 07:54:14.117400000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri08',
      Start = '22-AUG-2047 05:50:43.869300000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   BreakPoint(
      Name = 'ETF-04->ETF-05',
      Mode = 'MATCH_ELEM',
      Time = '23-AUG-2047 02:33:40.761450000 ET',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo10',
      Start = '27-AUG-2047 04:20:55.786300000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'ETF-05',
      Active = True,
      Propagator = 'DIVA',
      Time = '28-AUG-2047 14:14:46.879700000 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  5.039850752047949e+01 *km ),
         Conic.bPlaneTheta(  1.679419127285657e+02 *deg ),
         Conic.vInfinity(  7.062757017453376e-01 *km/sec ),
         Conic.inboundDec( -2.557224434380545e+01 *deg ),
         Conic.inboundRA(  5.332774757257473e+01 *deg ),
         Conic.trueAnomaly(  1.336447075840307e-07 *deg ),
         ],
      Controls = [
         [ -8.640000000000000e+04 *sec, 'Cosmic/Cosmic/ETF-05/TIME',  8.640000000000000e+04 *sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-05/Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.000000000000000e+01 *km, 'Cosmic/Cosmic/ETF-05/Conic.periapsisAltitude',  2.000000000000000e+02 *km ],
         [  2.000000000000000e-01 *km/sec, 'Cosmic/Cosmic/ETF-05/Conic.vInfinity',  1.000000000000000e+00 *km/sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-05/Conic.inboundDec',  3.600000000000000e+02 *deg ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-05/Conic.inboundRA',  3.600000000000000e+02 *deg ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri09',
      Start = '31-AUG-2047 14:26:17.276900000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo11',
      Start = '04-SEP-2047 07:08:47.434400000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri10',
      Start = '08-SEP-2047 08:37:09.115500000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DV02',
      Start = '09-SEP-2047 02:53:56.000000000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.365282953619952e-02 *km/sec,
      # DeltaVelRA  =  2.537124012063888e+02 *deg,
      # DeltaVelDec =  2.852338541161988e+01 *deg,
      DeltaVel = [
         -3.364296023968816e-03 *km/sec,
         -1.151424587408378e-02 *km/sec,
          6.519463824485427e-03 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -4.000000000000000e-02 *km/sec, 'DX',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DY',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DZ',  4.000000000000000e-02 *km/sec ],
         ],
      ),
   BreakPoint(
      Name = 'ETF-05->ETF-06',
      Mode = 'MATCH_ELEM',
      Time = '09-SEP-2047 03:37:16.359400000 ET',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo12',
      Start = '11-SEP-2047 04:55:17.301500000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri11',
      Start = '16-SEP-2047 03:24:30.885400000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo13',
      Start = '18-SEP-2047 00:24:46.473300000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'ETF-06',
      Active = True,
      Propagator = 'DIVA',
      Time = '20-SEP-2047 16:59:45.839100000 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  1.938073783668421e+02 *km ),
         Conic.bPlaneTheta(  7.382095334008017e+00 *deg ),
         Conic.vInfinity(  6.594444211613025e-01 *km/sec ),
         Conic.inboundDec( -2.696485764540171e+01 *deg ),
         Conic.inboundRA(  9.686582070526300e+01 *deg ),
         Conic.trueAnomaly( -3.704541763473408e-06 *deg ),
         ],
      Controls = [
         [ -8.640000000000000e+04 *sec, 'Cosmic/Cosmic/ETF-06/TIME',  8.640000000000000e+04 *sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-06/Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.000000000000000e+01 *km, 'Cosmic/Cosmic/ETF-06/Conic.periapsisAltitude',  2.000000000000000e+02 *km ],
         [  2.000000000000000e-01 *km/sec, 'Cosmic/Cosmic/ETF-06/Conic.vInfinity',  1.000000000000000e+00 *km/sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-06/Conic.inboundDec',  3.600000000000000e+02 *deg ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-06/Conic.inboundRA',  3.600000000000000e+02 *deg ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri12',
      Start = '25-SEP-2047 10:44:44.812900000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo14',
      Start = '26-SEP-2047 08:20:27.708100000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   BreakPoint(
      Name = 'ETF-06->ETF-07',
      Mode = 'MATCH_ELEM',
      Time = '26-SEP-2047 21:33:38.147250000 ET',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DV03',
      Start = '27-SEP-2047 17:45:25.000000000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.437195774862945e-02 *km/sec,
      # DeltaVelRA  =  2.625668033104578e+02 *deg,
      # DeltaVelDec =  2.745071814919245e+01 *deg,
      DeltaVel = [
         -1.649959029291104e-03 *km/sec,
         -1.264660746800494e-02 *km/sec,
          6.625264091359107e-03 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -4.000000000000000e-02 *km/sec, 'DX',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DY',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DZ',  4.000000000000000e-02 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'ETF-07',
      Active = True,
      Propagator = 'DIVA',
      Time = '03-OCT-2047 02:07:30.455400000 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  4.953123292533633e+01 *km ),
         Conic.bPlaneTheta(  1.006765553999817e+01 *deg ),
         Conic.vInfinity(  5.620869258052134e-01 *km/sec ),
         Conic.inboundDec( -2.626656155820896e+01 *deg ),
         Conic.inboundRA(  1.020176681159576e+02 *deg ),
         Conic.trueAnomaly( -5.717760590056391e-07 *deg ),
         ],
      Controls = [
         [ -8.640000000000000e+04 *sec, 'Cosmic/Cosmic/ETF-07/TIME',  8.640000000000000e+04 *sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-07/Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.000000000000000e+01 *km, 'Cosmic/Cosmic/ETF-07/Conic.periapsisAltitude',  2.000000000000000e+02 *km ],
         [  2.000000000000000e-01 *km/sec, 'Cosmic/Cosmic/ETF-07/Conic.vInfinity',  1.000000000000000e+00 *km/sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-07/Conic.inboundDec',  3.600000000000000e+02 *deg ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-07/Conic.inboundRA',  3.600000000000000e+02 *deg ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo15',
      Start = '04-OCT-2047 11:42:43.718200000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri13',
      Start = '04-OCT-2047 17:15:13.353500000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   BreakPoint(
      Name = 'ETF-07->ETF-08',
      Mode = 'MATCH_ELEM',
      Time = '09-OCT-2047 22:44:41.122950000 ET',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DV04',
      Start = '10-OCT-2047 03:26:05.000000000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  2.930963836332256e-03 *km/sec,
      # DeltaVelRA  =  2.987877551179990e+02 *deg,
      # DeltaVelDec =  1.997130779753551e+01 *deg,
      DeltaVel = [
          1.326574166689021e-03 *km/sec,
         -2.414251482456109e-03 *km/sec,
          1.001069313094405e-03 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -4.000000000000000e-02 *km/sec, 'DX',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DY',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DZ',  4.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo16',
      Start = '11-OCT-2047 08:36:57.260200000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri14',
      Start = '12-OCT-2047 08:08:08.986100000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'ETF-08',
      Active = True,
      Propagator = 'DIVA',
      Time = '16-OCT-2047 19:21:51.790500000 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  1.851276468677661e+02 *km ),
         Conic.bPlaneTheta(  1.645657385396549e+01 *deg ),
         Conic.vInfinity(  5.571208270211176e-01 *km/sec ),
         Conic.inboundDec( -2.334774827085500e+01 *deg ),
         Conic.inboundRA(  1.157545907725322e+02 *deg ),
         Conic.trueAnomaly( -2.690852681843373e-06 *deg ),
         ],
      Controls = [
         [ -8.640000000000000e+04 *sec, 'Cosmic/Cosmic/ETF-08/TIME',  8.640000000000000e+04 *sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-08/Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.000000000000000e+01 *km, 'Cosmic/Cosmic/ETF-08/Conic.periapsisAltitude',  2.000000000000000e+02 *km ],
         [  2.000000000000000e-01 *km/sec, 'Cosmic/Cosmic/ETF-08/Conic.vInfinity',  1.000000000000000e+00 *km/sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-08/Conic.inboundDec',  3.600000000000000e+02 *deg ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-08/Conic.inboundRA',  3.600000000000000e+02 *deg ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo17',
      Start = '19-OCT-2047 16:02:16.002700000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri15',
      Start = '21-OCT-2047 10:48:07.706900000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   BreakPoint(
      Name = 'ETF-08->ETF-09',
      Mode = 'MATCH_ELEM',
      Time = '24-OCT-2047 09:06:45.070950000 ET',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DV05',
      Start = '25-OCT-2047 04:45:36.000000000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  2.096521785848085e-02 *km/sec,
      # DeltaVelRA  =  2.841902641506967e+02 *deg,
      # DeltaVelDec =  2.432929566152896e+01 *deg,
      DeltaVel = [
          4.683046969911953e-03 *km/sec,
         -1.852045418112347e-02 *km/sec,
          8.637256963647033e-03 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -4.000000000000000e-02 *km/sec, 'DX',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DY',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DZ',  4.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo18',
      Start = '26-OCT-2047 10:50:17.155700000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri16',
      Start = '29-OCT-2047 00:00:00.589800000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'ETF-09',
      Active = True,
      Propagator = 'DIVA',
      Time = '31-OCT-2047 22:51:38.351400000 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  5.147255188496728e+01 *km ),
         Conic.bPlaneTheta(  1.939926945410127e+01 *deg ),
         Conic.vInfinity(  4.294345596189092e-01 *km/sec ),
         Conic.inboundDec( -2.175479750658329e+01 *deg ),
         Conic.inboundRA(  1.209272799084547e+02 *deg ),
         Conic.trueAnomaly(  4.243574282733001e-06 *deg ),
         ],
      Controls = [
         [ -8.640000000000000e+04 *sec, 'Cosmic/Cosmic/ETF-09/TIME',  8.640000000000000e+04 *sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-09/Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.000000000000000e+01 *km, 'Cosmic/Cosmic/ETF-09/Conic.periapsisAltitude',  2.000000000000000e+02 *km ],
         [  2.000000000000000e-01 *km/sec, 'Cosmic/Cosmic/ETF-09/Conic.vInfinity',  1.000000000000000e+00 *km/sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-09/Conic.inboundDec',  3.600000000000000e+02 *deg ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-09/Conic.inboundRA',  3.600000000000000e+02 *deg ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo19',
      Start = '03-NOV-2047 18:44:37.930600000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri17',
      Start = '07-NOV-2047 00:09:18.084000000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   BreakPoint(
      Name = 'ETF-09->ETF-10',
      Mode = 'MATCH_ELEM',
      Time = '09-NOV-2047 22:46:59.969950000 ET',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo20',
      Start = '10-NOV-2047 14:54:24.358200000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri18',
      Start = '14-NOV-2047 10:39:05.046600000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DV06',
      Start = '17-NOV-2047 09:42:31.000000000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.768704498043583e-02 *km/sec,
      # DeltaVelRA  =  2.979668652388782e+02 *deg,
      # DeltaVelDec =  2.091520396969117e+01 *deg,
      DeltaVel = [
          7.748003382727298e-03 *km/sec,
         -1.459222692179907e-02 *km/sec,
          6.314025434216219e-03 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -4.000000000000000e-02 *km/sec, 'DX',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DY',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DZ',  4.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo21',
      Start = '17-NOV-2047 14:39:34.583000000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'ETF-10',
      Active = True,
      Propagator = 'DIVA',
      Time = '18-NOV-2047 22:42:21.588500000 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  4.960524141503497e+01 *km ),
         Conic.bPlaneTheta( -1.628947043652629e+02 *deg ),
         Conic.vInfinity(  3.718909710569697e-01 *km/sec ),
         Conic.inboundDec( -2.518195253775352e+01 *deg ),
         Conic.inboundRA(  1.067178216306318e+02 *deg ),
         Conic.trueAnomaly(  1.069588794724196e-06 *deg ),
         ],
      Controls = [
         [ -8.640000000000000e+04 *sec, 'Cosmic/Cosmic/ETF-10/TIME',  8.640000000000000e+04 *sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-10/Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.000000000000000e+01 *km, 'Cosmic/Cosmic/ETF-10/Conic.periapsisAltitude',  2.000000000000000e+02 *km ],
         [  2.000000000000000e-01 *km/sec, 'Cosmic/Cosmic/ETF-10/Conic.vInfinity',  1.000000000000000e+00 *km/sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-10/Conic.inboundDec',  3.600000000000000e+02 *deg ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-10/Conic.inboundRA',  3.600000000000000e+02 *deg ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri19',
      Start = '23-NOV-2047 05:36:55.661200000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo22',
      Start = '25-NOV-2047 18:26:29.239900000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DV07',
      Start = '29-NOV-2047 18:42:20.000000000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  2.409271857049694e-02 *km/sec,
      # DeltaVelRA  =  2.941780251139279e+02 *deg,
      # DeltaVelDec =  2.575525199689197e+01 *deg,
      DeltaVel = [
          8.887455140278888e-03 *km/sec,
         -1.979578694150064e-02 *km/sec,
          1.046895642427703e-02 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -4.000000000000000e-02 *km/sec, 'DX',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DY',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DZ',  4.000000000000000e-02 *km/sec ],
         ],
      ),
   BreakPoint(
      Name = 'ETF-10->ETF-11',
      Mode = 'MATCH_ELEM',
      Time = '29-NOV-2047 20:15:28.419750000 ET',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri20',
      Start = '30-NOV-2047 12:48:42.588500000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo23',
      Start = '02-DEC-2047 11:30:52.139800000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri21',
      Start = '07-DEC-2047 20:53:55.039200000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo24',
      Start = '09-DEC-2047 12:02:36.442100000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'ETF-11',
      Active = True,
      Propagator = 'DIVA',
      Time = '10-DEC-2047 17:48:35.251000000 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  4.967509453923496e+01 *km ),
         Conic.bPlaneTheta(  7.772425781397184e+00 *deg ),
         Conic.vInfinity(  2.647540613888633e-01 *km/sec ),
         Conic.inboundDec( -2.092557433015346e+01 *deg ),
         Conic.inboundRA(  1.283172064374212e+02 *deg ),
         Conic.trueAnomaly(  2.394990728148501e-06 *deg ),
         ],
      Controls = [
         [ -8.640000000000000e+04 *sec, 'Cosmic/Cosmic/ETF-11/TIME',  8.640000000000000e+04 *sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-11/Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.000000000000000e+01 *km, 'Cosmic/Cosmic/ETF-11/Conic.periapsisAltitude',  2.000000000000000e+02 *km ],
         [  2.000000000000000e-01 *km/sec, 'Cosmic/Cosmic/ETF-11/Conic.vInfinity',  1.000000000000000e+00 *km/sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-11/Conic.inboundDec',  3.600000000000000e+02 *deg ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-11/Conic.inboundRA',  3.600000000000000e+02 *deg ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri22',
      Start = '16-DEC-2047 14:18:48.380300000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo25',
      Start = '17-DEC-2047 17:51:56.595000000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri23',
      Start = '23-DEC-2047 17:53:10.620200000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo26',
      Start = '24-DEC-2047 15:41:25.706300000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   BreakPoint(
      Name = 'ETF-11->ETF-Last',
      Mode = 'MATCH_ELEM',
      Time = '27-DEC-2047 06:09:16.312610000 ET',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DV08',
      Start = '27-DEC-2047 06:54:04.000000000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.214356160804756e-02 *km/sec,
      # DeltaVelRA  =  3.189959069477645e+02 *deg,
      # DeltaVelDec =  1.073069600752149e+01 *deg,
      DeltaVel = [
          9.004038831712080e-03 *km/sec,
         -7.828220895017779e-03 *km/sec,
          2.261046409137712e-03 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -4.000000000000000e-02 *km/sec, 'DX',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DY',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DZ',  4.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri24',
      Start = '30-DEC-2047 21:46:25.699600000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo27',
      Start = '31-DEC-2047 08:38:11.090000000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri25',
      Start = '07-JAN-2048 01:53:13.797600000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo28',
      Start = '07-JAN-2048 07:05:19.338700000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.732050807568877e-05 *km/sec,
      # DeltaVelRA  =  4.500000000000000e+01 *deg,
      # DeltaVelDec =  3.526438968275467e+01 *deg,
      DeltaVel = [
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
          1.000000000000000e-05 *km/sec,
         ],
      Controls = [
         [ -1.800000000000000e+04 *sec, 'TIME',  1.800000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'ETF-Last',
      Active = True,
      Propagator = 'DIVA',
      Time = '12-JAN-2048 18:29:57.374220000 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  1.072637741975505e+02 *km ),
         Conic.bPlaneTheta(  1.167483611222127e+02 *deg ),
         Conic.vInfinity(  1.655778434328415e-01 *km/sec ),
         Conic.inboundDec( -1.892171314948090e+01 *deg ),
         Conic.inboundRA(  1.508085090090735e+02 *deg ),
         Conic.trueAnomaly( -5.164288644538340e-08 *deg ),
         ],
      Controls = [
         [ -8.640000000000000e+04 *sec, 'Cosmic/Cosmic/ETF-Last/TIME',  8.640000000000000e+04 *sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-Last/Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.000000000000000e+01 *km, 'Cosmic/Cosmic/ETF-Last/Conic.periapsisAltitude',  2.000000000000000e+02 *km ],
         [  2.000000000000000e-01 *km/sec, 'Cosmic/Cosmic/ETF-Last/Conic.vInfinity',  1.000000000000000e+00 *km/sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-Last/Conic.inboundDec',  3.600000000000000e+02 *deg ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-Last/Conic.inboundRA',  3.600000000000000e+02 *deg ],
         ],
      ),
   ]


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
