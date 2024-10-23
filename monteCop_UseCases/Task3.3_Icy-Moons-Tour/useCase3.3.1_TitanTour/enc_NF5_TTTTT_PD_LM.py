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
# PosTolIn     = 1.0E-1
# VelTolIn     = 1.0E-4
# MajorFeasTol = 1e-5
# MajorOptTol  = 1e-3

PosTolIn     = 1.0E-3
VelTolIn     = 1.0E-6
MajorFeasTol = 1e-8
MajorOptTol  = 1e-3

#--------------------------------------------------------------------------
# Setup SNOPT params for cosmicBatch() run:
#--------------------------------------------------------------------------
MaxIterLimit   = 300                 #Major iterations limit

#set SNOPT params (these values are overwritten on cosmicBatch.py calls):
MajorStepLimit = 0.5
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
      "/nav/common/import/ephem/de440.boa",
      "/group/encnf5/import/ephem/sat441.bsp",
      "/group/encnf5/import/ephem/jup380.boa",
      #"/group/encnf5/lockfile/lockfile_0.6.2.boa",
      #"/group/encnf5/lockfile/locksmith/v1.0.0/lockfile_v1.0.0.boa"
      "/group/encnf5/lockfile/locksmith/v1.2.0/lockfile_v1.2.0.boa"
      #"/group/encnf5/lockfile/locksmith/v1.2.1/lockfile_v1.2.1.boa"
   ]
)

Defaults(
   # In order to reuse lockfile force model and integrator, make sure Body
   # is 'Europa Clipper'
   Body  = scName,
   Frame = 'EMO2000',
   #Mass  = 4500.000 *kg,
   #Isp   = 315.0 *sec,

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

Timeline = [
   ControlPoint(
      Name = 'TTF-01',
      Active = True,
      Propagator = 'DIVA',
      Time = '19-OCT-2044 19:23:06.345200000 ET',
      Body = 'encnf5',
      Center = 'Titan',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  1.732109009161299e+03 *km ),
         Conic.bPlaneTheta(  1.749342385845971e+02 *deg ),
         Conic.vInfinity(  2.283735044876229e+00 *km/sec ),
         Conic.inboundDec(  2.450324193119108e+01 *deg ),
         Conic.inboundRA(  2.699870187833350e+02 *deg ),
         Conic.trueAnomaly( -3.288236804784989e-07 *deg ),
         ],
      Controls = [
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo01',
      Start = '11-NOV-2044 07:38:41.734100000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.626441988559945e-03 *km/sec,
      # DeltaVelRA  =  2.612100139992911e+02 *deg,
      # DeltaVelDec = -3.187927431923592e+00 *deg,
      DeltaVel = [
         -2.481570479300439e-04 *km/sec,
         -1.604852246169097e-03 *km/sec,
         -9.044827073792280e-05 *km/sec,
         ],
      Controls = [
         [ -3.600000000000000e+04 *sec, 'TIME',  0.000000000000000e+00 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   BreakPoint(
      Name = 'TTF-01->TTF-02',
      Mode = 'MATCH_ELEM',
      Time = '12-NOV-2044 16:26:34.807300000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri01',
      Start = '04-DEC-2044 19:23:27.035800000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  3.573912611075197e-14 *km/sec,
      # DeltaVelRA  =  2.459071623360718e+01 *deg,
      # DeltaVelDec =  1.183709739133726e+01 *deg,
      DeltaVel = [
          3.180664307645658e-14 *km/sec,
          1.455598463313765e-14 *km/sec,
          7.331159712282263e-15 *km/sec,
         ],
      Controls = [
         [  0.000000000000000e+00 *sec, 'TIME',  3.600000000000000e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'TTF-02',
      Active = True,
      Propagator = 'DIVA',
      Time = '06-DEC-2044 14:10:00.118563615 ET',
      Body = 'encnf5',
      Center = 'Titan',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  3.063014412247131e+03 *km ),
         Conic.bPlaneTheta( -1.699182390276931e+02 *deg ),
         Conic.vInfinity(  2.294131463094499e+00 *km/sec ),
         Conic.inboundDec(  2.292803098363881e+01 *deg ),
         Conic.inboundRA( -1.264457324023122e+02 *deg ),
         Conic.trueAnomaly( -5.214787178610423e-08 *deg ),
         ],
      Controls = [
         [ -8.879684916361500e+04 *sec, 'Cosmic/Cosmic/TTF-02/TIME',  8.400315083638500e+04 *sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/TTF-02/Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  1.500000000000000e+03 *km, 'Cosmic/Cosmic/TTF-02/Conic.periapsisAltitude',  5.000000000000000e+03 *km ],
         [  2.000000000000000e+00 *km/sec, 'Cosmic/Cosmic/TTF-02/Conic.vInfinity',  3.000000000000000e+00 *km/sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/TTF-02/Conic.inboundDec',  3.600000000000000e+02 *deg ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/TTF-02/Conic.inboundRA',  3.600000000000000e+02 *deg ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo02',
      Start = '14-DEC-2044 15:44:08.337900000 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.524160245942347e-13 *km/sec,
      # DeltaVelRA  =  2.392794621336713e+02 *deg,
      # DeltaVelDec =  3.990804113414500e+00 *deg,
      DeltaVel = [
         -7.767309720220821e-14 *km/sec,
         -1.307096527353893e-13 *km/sec,
          1.060760132737960e-14 *km/sec,
         ],
      Controls = [
         [ -3.600000000000000e+04 *sec, 'TIME',  0.000000000000000e+00 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri02',
      Start = '25-DEC-2044 02:14:44.441108125 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.090206168982551e-13 *km/sec,
      # DeltaVelRA  =  2.264314716865740e+02 *deg,
      # DeltaVelDec = -6.075000356862639e+00 *deg,
      DeltaVel = [
         -7.471740438116156e-14 *km/sec,
         -7.854735172499288e-14 *km/sec,
         -1.153767430673692e-14 *km/sec,
         ],
      Controls = [
         [ -1.875205990812500e+04 *sec, 'TIME',  1.724794009187500e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   BreakPoint(
      Name = 'TTF-02->TTF-03',
      Mode = 'MATCH_ELEM',
      Time = '07-JAN-2045 10:11:45.145150000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DV01',
      Start = '26-JAN-2045 03:15:34.441814175 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  6.987015237508353e-14 *km/sec,
      # DeltaVelRA  =  2.509844370515237e+02 *deg,
      # DeltaVelDec =  6.988110110103534e+01 *deg,
      DeltaVel = [
         -7.830615588337675e-15 *km/sec,
         -2.272170783990143e-14 *km/sec,
          6.560673475538342e-14 *km/sec,
         ],
      Controls = [
         [ -3.301411481417500e+04 *sec, 'TIME',  2.985885185825000e+03 *sec ],
         [ -4.000000000000000e-02 *km/sec, 'DX',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DY',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DZ',  4.000000000000000e-02 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'TTF-03',
      Active = True,
      Propagator = 'DIVA',
      Time = '08-FEB-2045 07:29:10.772622485 ET',
      Body = 'encnf5',
      Center = 'Titan',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  1.694643475349066e+03 *km ),
         Conic.bPlaneTheta( -1.602182708579391e+02 *deg ),
         Conic.vInfinity(  2.300817708203693e+00 *km/sec ),
         Conic.inboundDec(  1.585161515959451e+01 *deg ),
         Conic.inboundRA( -1.542615803453485e+02 *deg ),
         Conic.trueAnomaly( -1.427675711635282e-06 *deg ),
         ],
      Controls = [
         [ -8.854375172248500e+04 *sec, 'Cosmic/Cosmic/TTF-03/TIME',  8.425624827751500e+04 *sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/TTF-03/Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  1.500000000000000e+03 *km, 'Cosmic/Cosmic/TTF-03/Conic.periapsisAltitude',  5.000000000000000e+03 *km ],
         [  2.000000000000000e+00 *km/sec, 'Cosmic/Cosmic/TTF-03/Conic.vInfinity',  3.000000000000000e+00 *km/sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/TTF-03/Conic.inboundDec',  3.600000000000000e+02 *deg ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/TTF-03/Conic.inboundRA',  3.600000000000000e+02 *deg ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo05',
      Start = '10-FEB-2045 16:47:59.131389543 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.021466535479363e-13 *km/sec,
      # DeltaVelRA  =  1.194607104681123e+02 *deg,
      # DeltaVelDec = -9.886382285160419e+00 *deg,
      DeltaVel = [
         -4.949240996855553e-14 *km/sec,
          8.761767195957030e-14 *km/sec,
         -1.753806565535925e-14 *km/sec,
         ],
      Controls = [
         [ -1.680632918954300e+04 *sec, 'TIME',  1.919367081045700e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri05',
      Start = '16-FEB-2045 09:59:57.418055068 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.962351646337587e-03 *km/sec,
      # DeltaVelRA  =  3.292333215644402e+02 *deg,
      # DeltaVelDec =  1.003799757856650e+01 *deg,
      DeltaVel = [
          1.660354263315112e-03 *km/sec,
         -9.884614866549216e-04 *km/sec,
          3.420403391050706e-04 *km/sec,
         ],
      Controls = [
         [ -1.552635145506800e+04 *sec, 'TIME',  2.047364854493200e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   BreakPoint(
      Name = 'TTF-03->TTF-04',
      Mode = 'MATCH_ELEM',
      Time = '05-MAR-2045 19:24:15.439500000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri08',
      Start = '23-MAR-2045 03:04:15.621920361 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  4.479573448923230e-14 *km/sec,
      # DeltaVelRA  =  2.361389338670085e+02 *deg,
      # DeltaVelDec =  1.855097182416417e+00 *deg,
      DeltaVel = [
         -2.494624928812587e-14 *km/sec,
         -3.717848327234770e-14 *km/sec,
          1.450122909355276e-15 *km/sec,
         ],
      Controls = [
         [ -3.595172932036100e+04 *sec, 'TIME',  4.827067963900000e+01 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo09',
      Start = '28-MAR-2045 11:51:52.619262498 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.459599715792402e-13 *km/sec,
      # DeltaVelRA  =  1.075272316922313e+02 *deg,
      # DeltaVelDec = -1.646465348426840e+01 *deg,
      DeltaVel = [
         -4.215469254499004e-14 *km/sec,
          1.334763655365315e-13 *km/sec,
         -4.136852720375944e-14 *km/sec,
         ],
      Controls = [
         [ -4.360068762498000e+03 *sec, 'TIME',  3.163993123750200e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'TTF-04',
      Active = True,
      Propagator = 'DIVA',
      Time = '31-MAR-2045 08:31:52.533563165 ET',
      Body = 'encnf5',
      Center = 'Titan',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  1.500000000000000e+03 *km ),
         Conic.bPlaneTheta(  3.274616816608582e+01 *deg ),
         Conic.vInfinity(  2.335879850980339e+00 *km/sec ),
         Conic.inboundDec( -1.717521100000940e+01 *deg ),
         Conic.inboundRA(  1.327250061450487e+02 *deg ),
         Conic.trueAnomaly(  7.431439729293590e-07 *deg ),
         ],
      Controls = [
         [ -8.860867546316500e+04 *sec, 'Cosmic/Cosmic/TTF-04/TIME',  8.419132453683500e+04 *sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/TTF-04/Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  1.500000000000000e+03 *km, 'Cosmic/Cosmic/TTF-04/Conic.periapsisAltitude',  5.000000000000000e+03 *km ],
         [  2.000000000000000e+00 *km/sec, 'Cosmic/Cosmic/TTF-04/Conic.vInfinity',  3.000000000000000e+00 *km/sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/TTF-04/Conic.inboundDec',  3.600000000000000e+02 *deg ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/TTF-04/Conic.inboundRA',  3.600000000000000e+02 *deg ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri09',
      Start = '03-APR-2045 09:35:05.541522939 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.453253464257066e-13 *km/sec,
      # DeltaVelRA  =  5.887170385996627e+00 *deg,
      # DeltaVelDec =  7.103109378795369e+01 *deg,
      DeltaVel = [
          4.698958210505250e-14 *km/sec,
          4.845266963400272e-15 *km/sec,
          1.374334707245565e-13 *km/sec,
         ],
      Controls = [
         [ -1.374095942293900e+04 *sec, 'TIME',  2.225904057706100e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo10',
      Start = '07-APR-2045 08:58:14.397326268 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  9.488267815691180e-14 *km/sec,
      # DeltaVelRA  =  1.285509872059332e+02 *deg,
      # DeltaVelDec =  7.445960559045264e+01 *deg,
      DeltaVel = [
         -1.584248504443155e-14 *km/sec,
          1.988039398412846e-14 *km/sec,
          9.141393885442688e-14 *km/sec,
         ],
      Controls = [
         [ -1.224646802626800e+04 *sec, 'TIME',  2.375353197373200e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DV03',
      Start = '07-APR-2045 11:54:11.836811478 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.543679113143803e-13 *km/sec,
      # DeltaVelRA  =  2.527725798111264e+02 *deg,
      # DeltaVelDec = -5.493049075856048e+01 *deg,
      DeltaVel = [
         -2.626841174991363e-14 *km/sec,
         -8.471598351688606e-14 *km/sec,
         -1.263432812794889e-13 *km/sec,
         ],
      Controls = [
         [ -2.271583681147800e+04 *sec, 'TIME',  1.328416318852200e+04 *sec ],
         [ -4.000000000000000e-02 *km/sec, 'DX',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DY',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DZ',  4.000000000000000e-02 *km/sec ],
         ],
      ),
   BreakPoint(
      Name = 'TTF-04->TTF-05',
      Mode = 'MATCH_ELEM',
      Time = '08-APR-2045 07:00:13.334900000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri10',
      Start = '11-APR-2045 08:28:06.094489013 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  2.331559720159270e-13 *km/sec,
      # DeltaVelRA  =  9.925858443332604e+01 *deg,
      # DeltaVelDec =  4.849805512291208e+01 *deg,
      DeltaVel = [
         -2.485753290035851e-14 *km/sec,
          1.524869681296706e-13 *km/sec,
          1.746182547599768e-13 *km/sec,
         ],
      Controls = [
         [ -1.224019968901300e+04 *sec, 'TIME',  2.375980031098700e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo11',
      Start = '15-APR-2045 07:48:53.547790694 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  7.296166321493116e-14 *km/sec,
      # DeltaVelRA  =  2.126116228850955e+02 *deg,
      # DeltaVelDec =  2.030049945027962e+01 *deg,
      DeltaVel = [
         -5.764129807893147e-14 *km/sec,
         -3.687962550043775e-14 *km/sec,
          2.531359867368271e-14 *km/sec,
         ],
      Controls = [
         [ -1.226839439069400e+04 *sec, 'TIME',  2.373160560930600e+04 *sec ],
         [ -1.000000000000000e-02 *km/sec, 'DX',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DY',  1.000000000000000e-02 *km/sec ],
         [ -1.000000000000000e-02 *km/sec, 'DZ',  1.000000000000000e-02 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'TTF-05',
      Active = True,
      Propagator = 'DIVA',
      Time = '16-APR-2045 06:14:46.709882950 ET',
      Body = 'encnf5',
      Center = 'Titan',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  1.542405181631449e+03 *km ),
         Conic.bPlaneTheta(  2.505177967558319e+01 *deg ),
         Conic.vInfinity(  2.340377525837590e+00 *km/sec ),
         Conic.inboundDec(  2.339208076380495e+00 *deg ),
         Conic.inboundRA(  1.606716593539945e+02 *deg ),
         Conic.trueAnomaly( -2.209313619938414e-07 *deg ),
         ],
      Controls = [
         [ -8.696389818295000e+04 *sec, 'Cosmic/Cosmic/TTF-05/TIME',  8.583610181705000e+04 *sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/TTF-05/Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  1.500000000000000e+03 *km, 'Cosmic/Cosmic/TTF-05/Conic.periapsisAltitude',  5.000000000000000e+03 *km ],
         [  2.000000000000000e+00 *km/sec, 'Cosmic/Cosmic/TTF-05/Conic.vInfinity',  3.000000000000000e+00 *km/sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/TTF-05/Conic.inboundDec',  3.600000000000000e+02 *deg ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/TTF-05/Conic.inboundRA',  3.600000000000000e+02 *deg ],
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
