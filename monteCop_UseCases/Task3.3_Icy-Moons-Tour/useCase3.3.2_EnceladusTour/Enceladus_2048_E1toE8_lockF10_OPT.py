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
      #"/group/encnf5/lockfile/lockfile_0.6.2.boa",
      "/group/encnf5/lockfile/locksmith/v1.0.0/lockfile_v1.0.0.boa",
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
#     StateTol = 1e-11,
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
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo01',
      Start = '22-JUN-2047 01:49:08.276143559 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  2.079372564945525e-09 *km/sec,
      # DeltaVelRA  =  1.480912492489416e+02 *deg,
      # DeltaVelDec =  5.849110115205955e+00 *deg,
      DeltaVel = [
         -1.755970680874360e-09 *km/sec,
          1.093367616250112e-09 *km/sec,
          2.119067896347647e-10 *km/sec,
         ],
      Controls = [
         [ -7.135362644355900e+04 *sec, 'TIME',  3.664637355644100e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri01',
      Start = '22-JUN-2047 14:22:39.780473943 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  4.387408006968845e-06 *km/sec,
      # DeltaVelRA  =  2.447880627047095e+02 *deg,
      # DeltaVelDec =  2.484426140616608e+01 *deg,
      DeltaVel = [
         -1.695934280746075e-06 *km/sec,
         -3.602095612469541e-06 *km/sec,
          1.843383610004251e-06 *km/sec,
         ],
      Controls = [
         [ -2.215194647394300e+04 *sec, 'TIME',  1.384805352605700e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo02',
      Start = '29-JUN-2047 00:21:56.944434675 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  7.171325468961806e-06 *km/sec,
      # DeltaVelRA  =  2.400753963681130e+02 *deg,
      # DeltaVelDec =  2.044500492648375e+01 *deg,
      DeltaVel = [
         -3.352133662626046e-06 *km/sec,
         -5.823750625507724e-06 *km/sec,
          2.505002503305723e-06 *km/sec,
         ],
      Controls = [
         [ -7.188246063467499e+04 *sec, 'TIME',  3.611753936532500e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri02',
      Start = '30-JUN-2047 16:03:10.827774641 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  2.273339809547358e-05 *km/sec,
      # DeltaVelRA  =  2.475517335872465e+02 *deg,
      # DeltaVelDec =  2.677152477545692e+01 *deg,
      DeltaVel = [
         -7.750238285504801e-06 *km/sec,
         -1.875861913497686e-05 *km/sec,
          1.023989275500760e-05 *km/sec,
         ],
      Controls = [
         [ -2.213091277464100e+04 *sec, 'TIME',  1.386908722535900e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   BreakPoint(
      Name = 'ETF-01->ETF-02',
      Mode = 'MATCH_ELEM',
      Time = '04-JUL-2047 17:12:57.888750000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo03',
      Start = '05-JUL-2047 17:38:15.467341166 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.255510580204827e-09 *km/sec,
      # DeltaVelRA  =  2.288016630524703e+00 *deg,
      # DeltaVelDec = -8.403252727206995e+00 *deg,
      DeltaVel = [
          1.241041277746859e-09 *km/sec,
          4.958538736188396e-11 *km/sec,
         -1.834792993633742e-10 *km/sec,
         ],
      Controls = [
         [ -7.195978954116600e+04 *sec, 'TIME',  3.604021045883400e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri03',
      Start = '08-JUL-2047 15:53:00.059581405 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.071760486865566e-04 *km/sec,
      # DeltaVelRA  =  2.359786772401826e+02 *deg,
      # DeltaVelDec =  1.743029732320631e+01 *deg,
      DeltaVel = [
         -5.721167280198072e-05 *km/sec,
         -8.475174077284557e-05 *km/sec,
          3.210408608315529e-05 *km/sec,
         ],
      Controls = [
         [ -1.550297428140500e+04 *sec, 'TIME',  2.049702571859500e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo04',
      Start = '12-JUL-2047 19:45:02.007735563 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  9.166027173418408e-09 *km/sec,
      # DeltaVelRA  =  8.863763821559547e+01 *deg,
      # DeltaVelDec = -7.946398290332789e+01 *deg,
      DeltaVel = [
          3.984864004301158e-11 *km/sec,
          1.675567107343342e-09 *km/sec,
         -9.011489394019196e-09 *km/sec,
         ],
      Controls = [
         [ -8.745913535563001e+03 *sec, 'TIME',  2.725408646443700e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri04',
      Start = '16-JUL-2047 17:37:21.271717512 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.310396639647292e-03 *km/sec,
      # DeltaVelRA  =  2.394156366864482e+02 *deg,
      # DeltaVelDec =  2.136084264273471e+01 *deg,
      DeltaVel = [
         -6.209367026807418e-04 *km/sec,
         -1.050600904343969e-03 *km/sec,
          4.772993863965428e-04 *km/sec,
         ],
      Controls = [
         [ -1.576396931751200e+04 *sec, 'TIME',  2.023603068248800e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DV00',
      Start = '16-JUL-2047 17:37:23.902956621 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.368678255457480e-03 *km/sec,
      # DeltaVelRA  =  2.349054315666430e+02 *deg,
      # DeltaVelDec =  2.029965352367862e+01 *deg,
      DeltaVel = [
         -7.380180384348392e-04 *km/sec,
         -1.050305084743602e-03 *km/sec,
          4.748355198051832e-04 *km/sec,
         ],
      Controls = [
         [ -1.704390295662100e+04 *sec, 'TIME',  1.895609704337900e+04 *sec ],
         [ -4.000000000000000e-02 *km/sec, 'DX',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DY',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DZ',  4.000000000000000e-02 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'ETF-02',
      Active = True,
      Propagator = 'DIVA',
      Time = '18-JUL-2047 10:25:40.496105982 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  4.000000000000000e+01 *km ),
         Conic.bPlaneTheta(  1.643790510125787e+02 *deg ),
         Conic.vInfinity(  7.652900276765036e-01 *km/sec ),
         Conic.inboundDec( -2.408787598592023e+01 *deg ),
         Conic.inboundRA(  4.636217101556048e+01 *deg ),
         Conic.trueAnomaly(  7.091252994054027e-06 *deg ),
         ],
      Controls = [
         [ -8.638471860598199e+04 *sec, 'Cosmic/Cosmic/ETF-02/TIME',  8.641528139401801e+04 *sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-02/Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.000000000000000e+01 *km, 'Cosmic/Cosmic/ETF-02/Conic.periapsisAltitude',  2.000000000000000e+02 *km ],
         [  2.000000000000000e-01 *km/sec, 'Cosmic/Cosmic/ETF-02/Conic.vInfinity',  8.000000000000000e-01 *km/sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-02/Conic.inboundDec',  3.600000000000000e+02 *deg ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-02/Conic.inboundRA',  3.600000000000000e+02 *deg ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo05',
      Start = '20-JUL-2047 22:52:15.163125608 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.166137248948827e-09 *km/sec,
      # DeltaVelRA  =  1.431694594855355e+02 *deg,
      # DeltaVelDec =  1.363069778146938e+01 *deg,
      DeltaVel = [
         -9.071011493800973e-10 *km/sec,
          6.793527722091679e-10 *km/sec,
          2.748152089513294e-10 *km/sec,
         ],
      Controls = [
         [ -7.155836352560800e+04 *sec, 'TIME',  3.644163647439200e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   BreakPoint(
      Name = 'ETF-02->ETF-03',
      Mode = 'MATCH_ELEM',
      Time = '23-JUL-2047 03:38:13.496250000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri05',
      Start = '26-JUL-2047 10:06:10.339880555 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  6.547224476129427e-10 *km/sec,
      # DeltaVelRA  =  2.387734775007992e+02 *deg,
      # DeltaVelDec = -1.119210161544333e+01 *deg,
      DeltaVel = [
         -3.329679387110689e-10 *km/sec,
         -5.492214866252241e-10 *km/sec,
         -1.270810521988353e-10 *km/sec,
         ],
      Controls = [
         [ -2.380076088055500e+04 *sec, 'TIME',  1.219923911944500e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DV01',
      Start = '27-JUL-2047 03:25:12.818013160 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  9.318317331557396e-03 *km/sec,
      # DeltaVelRA  =  2.433094292586773e+02 *deg,
      # DeltaVelDec =  2.009113000445361e+01 *deg,
      DeltaVel = [
         -3.930826926418759e-03 *km/sec,
         -7.818784717718668e-03 *km/sec,
          3.200975336380026e-03 *km/sec,
         ],
      Controls = [
         [ -1.779481801316000e+04 *sec, 'TIME',  1.820518198684000e+04 *sec ],
         [ -4.000000000000000e-02 *km/sec, 'DX',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DY',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DZ',  4.000000000000000e-02 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'ETF-03',
      Active = True,
      Propagator = 'DIVA',
      Time = '27-JUL-2047 20:54:28.704448659 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  4.000042363449073e+01 *km ),
         Conic.bPlaneTheta(  2.816052969315793e+00 *deg ),
         Conic.vInfinity(  7.315662562160318e-01 *km/sec ),
         Conic.inboundDec( -2.792741131464248e+01 *deg ),
         Conic.inboundRA(  8.054578818574910e+01 *deg ),
         Conic.trueAnomaly( -4.244185921666981e-06 *deg ),
         ],
      Controls = [
         [ -8.663748944865901e+04 *sec, 'Cosmic/Cosmic/ETF-03/TIME',  8.616251055134099e+04 *sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-03/Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.000000000000000e+01 *km, 'Cosmic/Cosmic/ETF-03/Conic.periapsisAltitude',  2.000000000000000e+02 *km ],
         [  2.000000000000000e-01 *km/sec, 'Cosmic/Cosmic/ETF-03/Conic.vInfinity',  8.000000000000000e-01 *km/sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-03/Conic.inboundDec',  3.600000000000000e+02 *deg ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-03/Conic.inboundRA',  3.600000000000000e+02 *deg ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo06',
      Start = '29-JUL-2047 02:08:22.855524844 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  2.580693546753974e-09 *km/sec,
      # DeltaVelRA  =  3.485389457200965e+02 *deg,
      # DeltaVelDec =  6.636466974786028e+01 *deg,
      DeltaVel = [
          1.014005553432375e-09 *km/sec,
         -2.055840795019033e-10 *km/sec,
          2.364213845260462e-09 *km/sec,
         ],
      Controls = [
         [ -7.148513512484400e+04 *sec, 'TIME',  3.651486487515600e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri06',
      Start = '04-AUG-2047 22:47:57.973573599 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  2.768372092044757e-09 *km/sec,
      # DeltaVelRA  =  3.145827040054404e+02 *deg,
      # DeltaVelDec =  7.125999587612985e+01 *deg,
      DeltaVel = [
          6.243084245099293e-10 *km/sec,
         -6.334692188608102e-10 *km/sec,
          2.621610150242417e-09 *km/sec,
         ],
      Controls = [
         [ -2.360324827359900e+04 *sec, 'TIME',  1.239675172640100e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo07',
      Start = '04-AUG-2047 23:03:33.829109053 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  9.669758249028980e-11 *km/sec,
      # DeltaVelRA  =  1.889623390045192e+02 *deg,
      # DeltaVelDec = -2.807317071962352e+01 *deg,
      DeltaVel = [
         -8.427916692705921e-11 *km/sec,
         -1.329172750724744e-11 *km/sec,
         -4.550576294734957e-11 *km/sec,
         ],
      Controls = [
         [ -7.173783080905300e+04 *sec, 'TIME',  3.626216919094700e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   BreakPoint(
      Name = 'ETF-03->ETF-04',
      Mode = 'MATCH_ELEM',
      Time = '07-AUG-2047 05:51:32.929100000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo08',
      Start = '11-AUG-2047 23:33:34.488055165 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.848030401497753e-10 *km/sec,
      # DeltaVelRA  =  2.138950754908812e+02 *deg,
      # DeltaVelDec =  4.352782162252906e+01 *deg,
      DeltaVel = [
         -1.112194390468981e-10 *km/sec,
         -7.472248665095578e-11 *km/sec,
          1.272750958176443e-10 *km/sec,
         ],
      Controls = [
         [ -1.000082175516500e+04 *sec, 'TIME',  2.599917824483500e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri07',
      Start = '12-AUG-2047 19:31:36.609837815 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  6.039222001166158e-05 *km/sec,
      # DeltaVelRA  =  2.540382769858369e+02 *deg,
      # DeltaVelDec =  2.974482652982379e+01 *deg,
      DeltaVel = [
         -1.441941268696706e-05 *km/sec,
         -5.041355462475643e-05 *km/sec,
          2.996288180059739e-05 *km/sec,
         ],
      Controls = [
         [ -1.892437423781500e+04 *sec, 'TIME',  1.707562576218500e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'ETF-04',
      Active = True,
      Propagator = 'DIVA',
      Time = '17-AUG-2047 14:51:33.405993394 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  4.000000000000000e+01 *km ),
         Conic.bPlaneTheta(  1.795156641265789e+02 *deg ),
         Conic.vInfinity(  6.988101451695857e-01 *km/sec ),
         Conic.inboundDec( -2.582802376901141e+01 *deg ),
         Conic.inboundRA(  5.605343992956716e+01 *deg ),
         Conic.trueAnomaly( -2.532841591086933e-06 *deg ),
         ],
      Controls = [
         [ -8.633876279339399e+04 *sec, 'Cosmic/Cosmic/ETF-04/TIME',  8.646123720660601e+04 *sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-04/Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.000000000000000e+01 *km, 'Cosmic/Cosmic/ETF-04/Conic.periapsisAltitude',  2.000000000000000e+02 *km ],
         [  2.000000000000000e-01 *km/sec, 'Cosmic/Cosmic/ETF-04/Conic.vInfinity',  8.000000000000000e-01 *km/sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-04/Conic.inboundDec',  3.600000000000000e+02 *deg ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-04/Conic.inboundRA',  3.600000000000000e+02 *deg ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo09',
      Start = '20-AUG-2047 02:45:45.890878459 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.084950543010652e-09 *km/sec,
      # DeltaVelRA  =  2.843261747133387e+02 *deg,
      # DeltaVelDec = -3.162608445563124e+01 *deg,
      DeltaVel = [
          2.285922264276768e-10 *km/sec,
         -8.950943835403626e-10 *km/sec,
         -5.689194313349172e-10 *km/sec,
         ],
      Controls = [
         [ -7.149177347845900e+04 *sec, 'TIME',  3.650822652154100e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri08',
      Start = '22-AUG-2047 06:31:18.897274991 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  4.694570982878985e-04 *km/sec,
      # DeltaVelRA  =  2.537001595861027e+02 *deg,
      # DeltaVelDec =  2.856393208119297e+01 *deg,
      DeltaVel = [
         -1.157224756265380e-04 *km/sec,
         -3.957440957671549e-04 *km/sec,
          2.244657800939731e-04 *km/sec,
         ],
      Controls = [
         [ -2.043502797499100e+04 *sec, 'TIME',  1.556497202500900e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   BreakPoint(
      Name = 'ETF-04->ETF-05',
      Mode = 'MATCH_ELEM',
      Time = '23-AUG-2047 02:33:40.761450000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo10',
      Start = '26-AUG-2047 22:23:52.464484557 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.355751703175909e-03 *km/sec,
      # DeltaVelRA  =  2.506910251932607e+02 *deg,
      # DeltaVelDec =  3.106334226510597e+01 *deg,
      DeltaVel = [
         -3.840090752827150e-04 *km/sec,
         -1.096007397806389e-03 *km/sec,
          6.995480646235543e-04 *km/sec,
         ],
      Controls = [
         [ -6.857667818455699e+04 *sec, 'TIME',  3.942332181544300e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'ETF-05',
      Active = True,
      Propagator = 'DIVA',
      Time = '28-AUG-2047 14:13:23.488127860 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  4.000000000000000e+01 *km ),
         Conic.bPlaneTheta(  1.746333237148080e+02 *deg ),
         Conic.vInfinity(  6.955906579287227e-01 *km/sec ),
         Conic.inboundDec( -2.646188548360593e+01 *deg ),
         Conic.inboundRA(  5.379388356364024e+01 *deg ),
         Conic.trueAnomaly(  1.336447075840307e-07 *deg ),
         ],
      Controls = [
         [ -8.631660842786000e+04 *sec, 'Cosmic/Cosmic/ETF-05/TIME',  8.648339157214000e+04 *sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-05/Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.000000000000000e+01 *km, 'Cosmic/Cosmic/ETF-05/Conic.periapsisAltitude',  2.000000000000000e+02 *km ],
         [  2.000000000000000e-01 *km/sec, 'Cosmic/Cosmic/ETF-05/Conic.vInfinity',  8.000000000000000e-01 *km/sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-05/Conic.inboundDec',  3.600000000000000e+02 *deg ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-05/Conic.inboundRA',  3.600000000000000e+02 *deg ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri09',
      Start = '31-AUG-2047 13:47:54.552747406 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.320634130502532e-08 *km/sec,
      # DeltaVelRA  =  2.818732646698740e+02 *deg,
      # DeltaVelDec =  2.429468063941222e+01 *deg,
      DeltaVel = [
          2.476543858078300e-09 *km/sec,
         -1.177928091845231e-08 *km/sec,
          5.433481593602546e-09 *km/sec,
         ],
      Controls = [
         [ -1.569727584740600e+04 *sec, 'TIME',  2.030272415259400e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo11',
      Start = '04-SEP-2047 10:26:22.379464085 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  4.354863713150638e-10 *km/sec,
      # DeltaVelRA  =  2.087055327861350e+02 *deg,
      # DeltaVelDec =  7.751983996182416e+01 *deg,
      DeltaVel = [
         -8.254322406457439e-11 *km/sec,
         -4.520145568050253e-11 *km/sec,
          4.251962184292690e-10 *km/sec,
         ],
      Controls = [
         [ -2.985494506408500e+04 *sec, 'TIME',  6.145054935915000e+03 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri10',
      Start = '08-SEP-2047 07:26:59.865346142 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  6.446146843602368e-10 *km/sec,
      # DeltaVelRA  =  9.823263917728931e+01 *deg,
      # DeltaVelDec =  9.289933547738926e+00 *deg,
      DeltaVel = [
         -9.109349522491723e-11 *km/sec,
          6.296042356323412e-10 *km/sec,
          1.040604290514456e-10 *km/sec,
         ],
      Controls = [
         [ -1.379074984614200e+04 *sec, 'TIME',  2.220925015385800e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DV02',
      Start = '09-SEP-2047 02:41:02.693291442 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  2.055154824433115e-02 *km/sec,
      # DeltaVelRA  =  2.516074778776871e+02 *deg,
      # DeltaVelDec =  2.820934547990970e+01 *deg,
      DeltaVel = [
         -5.714339868295592e-03 *km/sec,
         -1.718542768515442e-02 *km/sec,
          9.714603974790037e-03 *km/sec,
         ],
      Controls = [
         [ -1.722669329144200e+04 *sec, 'TIME',  1.877330670855800e+04 *sec ],
         [ -4.000000000000000e-02 *km/sec, 'DX',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DY',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DZ',  4.000000000000000e-02 *km/sec ],
         ],
      ),
   BreakPoint(
      Name = 'ETF-05->ETF-06',
      Mode = 'MATCH_ELEM',
      Time = '09-SEP-2047 03:37:16.359400000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo12',
      Start = '11-SEP-2047 09:57:18.693762434 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  5.699732766258416e-10 *km/sec,
      # DeltaVelRA  =  1.768024844467806e+02 *deg,
      # DeltaVelDec =  2.181415992790226e+00 *deg,
      DeltaVel = [
         -5.686735249427651e-10 *km/sec,
          3.176905031007503e-11 *km/sec,
          2.169528840933733e-11 *km/sec,
         ],
      Controls = [
         [ -3.612139226243400e+04 *sec, 'TIME',  7.187860773756600e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri11',
      Start = '16-SEP-2047 02:37:44.667688013 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  7.111587003670183e-10 *km/sec,
      # DeltaVelRA  =  2.024014939595779e+02 *deg,
      # DeltaVelDec = -5.595949708327535e+00 *deg,
      DeltaVel = [
         -6.543584773268839e-10 *km/sec,
         -2.697270626912402e-10 *km/sec,
         -6.934689546118295e-11 *km/sec,
         ],
      Controls = [
         [ -1.519378228801300e+04 *sec, 'TIME',  2.080621771198700e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo13',
      Start = '17-SEP-2047 19:26:48.638104111 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  9.312203318183270e-10 *km/sec,
      # DeltaVelRA  =  4.616380326672672e+01 *deg,
      # DeltaVelDec =  8.694975170076552e+01 *deg,
      DeltaVel = [
          3.431956019561540e-11 *km/sec,
          3.574287812685406e-11 *km/sec,
          9.299010274448693e-10 *km/sec,
         ],
      Controls = [
         [ -7.212216480411100e+04 *sec, 'TIME',  3.587783519588900e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'ETF-06',
      Active = True,
      Propagator = 'DIVA',
      Time = '20-SEP-2047 17:35:21.106607769 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  6.153290294019016e+01 *km ),
         Conic.bPlaneTheta( -3.055234341125120e+01 *deg ),
         Conic.vInfinity(  5.938093955051477e-01 *km/sec ),
         Conic.inboundDec( -2.858001048873222e+01 *deg ),
         Conic.inboundRA(  9.437472789344294e+01 *deg ),
         Conic.trueAnomaly( -3.704541763473407e-06 *deg ),
         ],
      Controls = [
         [ -8.853526750776899e+04 *sec, 'Cosmic/Cosmic/ETF-06/TIME',  8.426473249223101e+04 *sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-06/Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.000000000000000e+01 *km, 'Cosmic/Cosmic/ETF-06/Conic.periapsisAltitude',  2.000000000000000e+02 *km ],
         [  2.000000000000000e-01 *km/sec, 'Cosmic/Cosmic/ETF-06/Conic.vInfinity',  8.000000000000000e-01 *km/sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-06/Conic.inboundDec',  3.600000000000000e+02 *deg ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-06/Conic.inboundRA',  3.600000000000000e+02 *deg ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri12',
      Start = '25-SEP-2047 10:51:44.570293407 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.367867629275197e-09 *km/sec,
      # DeltaVelRA  =  2.012641102305900e+02 *deg,
      # DeltaVelDec = -6.376723710427962e+01 *deg,
      DeltaVel = [
         -5.634596688431256e-10 *km/sec,
         -2.192772446729161e-10 *km/sec,
         -1.226985143665072e-09 *km/sec,
         ],
      Controls = [
         [ -1.841975739340700e+04 *sec, 'TIME',  1.758024260659300e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo14',
      Start = '26-SEP-2047 05:18:19.856516240 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.290579684547180e-02 *km/sec,
      # DeltaVelRA  =  2.642061427779134e+02 *deg,
      # DeltaVelDec =  3.604495893769786e+01 *deg,
      DeltaVel = [
         -1.053414808218280e-03 *km/sec,
         -1.038174621607950e-02 *km/sec,
          7.594027585104978e-03 *km/sec,
         ],
      Controls = [
         [ -7.072148416240000e+03 *sec, 'TIME',  2.892785158376000e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   BreakPoint(
      Name = 'ETF-06->ETF-07',
      Mode = 'MATCH_ELEM',
      Time = '26-SEP-2047 21:33:38.147250000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DV03',
      Start = '27-SEP-2047 19:38:03.586521052 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.852163739870590e-04 *km/sec,
      # DeltaVelRA  =  2.738643937966798e+02 *deg,
      # DeltaVelDec =  3.715587174004913e+01 *deg,
      DeltaVel = [
          9.948659833910016e-06 *km/sec,
         -1.472809579763771e-04 *km/sec,
          1.118679971125965e-04 *km/sec,
         ],
      Controls = [
         [ -2.475858652105200e+04 *sec, 'TIME',  1.124141347894800e+04 *sec ],
         [ -4.000000000000000e-02 *km/sec, 'DX',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DY',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DZ',  4.000000000000000e-02 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'ETF-07',
      Active = True,
      Propagator = 'DIVA',
      Time = '03-OCT-2047 02:42:55.722202448 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  7.621757180637350e+01 *km ),
         Conic.bPlaneTheta(  2.815796453273849e+01 *deg ),
         Conic.vInfinity(  5.127149809565379e-01 *km/sec ),
         Conic.inboundDec( -3.313497183217454e+01 *deg ),
         Conic.inboundRA(  1.011526524618190e+02 *deg ),
         Conic.trueAnomaly( -5.717760590056391e-07 *deg ),
         ],
      Controls = [
         [ -8.852526680244799e+04 *sec, 'Cosmic/Cosmic/ETF-07/TIME',  8.427473319755201e+04 *sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-07/Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.000000000000000e+01 *km, 'Cosmic/Cosmic/ETF-07/Conic.periapsisAltitude',  2.000000000000000e+02 *km ],
         [  2.000000000000000e-01 *km/sec, 'Cosmic/Cosmic/ETF-07/Conic.vInfinity',  8.000000000000000e-01 *km/sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-07/Conic.inboundDec',  3.600000000000000e+02 *deg ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-07/Conic.inboundRA',  3.600000000000000e+02 *deg ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri13',
      Start = '04-OCT-2047 16:51:53.694782221 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  2.895221063388400e-09 *km/sec,
      # DeltaVelRA  =  3.122051593065601e+02 *deg,
      # DeltaVelDec =  5.507375816046149e+01 *deg,
      DeltaVel = [
          1.113538600288676e-09 *km/sec,
         -1.227839769651624e-09 *km/sec,
          2.373762054549032e-09 *km/sec,
         ],
      Controls = [
         [ -1.660034128222100e+04 *sec, 'TIME',  1.939965871777900e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo15',
      Start = '04-OCT-2047 17:01:49.992547386 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  7.519040557750525e-10 *km/sec,
      # DeltaVelRA  =  2.437078303632005e+02 *deg,
      # DeltaVelDec = -3.833208440999536e+01 *deg,
      DeltaVel = [
         -2.612579834725081e-10 *km/sec,
         -5.287972889301284e-10 *km/sec,
         -4.663447248366794e-10 *km/sec,
         ],
      Controls = [
         [ -3.714627434738600e+04 *sec, 'TIME',  7.085372565261400e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   BreakPoint(
      Name = 'ETF-07->ETF-08',
      Mode = 'MATCH_ELEM',
      Time = '09-OCT-2047 22:44:41.122950000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DV04',
      Start = '10-OCT-2047 01:00:43.541369621 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.029960375324610e-02 *km/sec,
      # DeltaVelRA  =  2.678350894228342e+02 *deg,
      # DeltaVelDec =  2.531881733963164e+01 *deg,
      DeltaVel = [
         -3.517022427166046e-04 *km/sec,
         -9.303600627489705e-03 *km/sec,
          4.404674604385948e-03 *km/sec,
         ],
      Controls = [
         [ -9.278541369621000e+03 *sec, 'TIME',  2.672145863037900e+04 *sec ],
         [ -4.000000000000000e-02 *km/sec, 'DX',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DY',  4.000000000000000e-02 *km/sec ],
         [ -4.000000000000000e-02 *km/sec, 'DZ',  4.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVApo16',
      Start = '10-OCT-2047 19:39:41.694614683 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.256685865835733e-02 *km/sec,
      # DeltaVelRA  =  2.708597115739079e+02 *deg,
      # DeltaVelDec =  3.040171342238947e+01 *deg,
      DeltaVel = [
          1.626293557524587e-04 *km/sec,
         -1.083767706116986e-02 *km/sec,
          6.359578928637550e-03 *km/sec,
         ],
      Controls = [
         [ -4.336443441468300e+04 *sec, 'TIME',  6.463556558531700e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   OptImpulseBurn(
      Body = 'encnf5',
      Name = 'DVPeri14',
      Start = '12-OCT-2047 08:33:59.912825571 ET',
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  2.457169940246198e-02 *km/sec,
      # DeltaVelRA  =  2.738797651904713e+02 *deg,
      # DeltaVelDec =  3.533214489426378e+01 *deg,
      DeltaVel = [
          1.356365770394802e-03 *km/sec,
         -1.999997779477012e-02 *km/sec,
          1.421019252616071e-02 *km/sec,
         ],
      Controls = [
         [ -1.955092672557100e+04 *sec, 'TIME',  1.644907327442900e+04 *sec ],
         [ -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec ],
         [ -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'ETF-08',
      Active = True,
      Propagator = 'DIVA',
      Time = '16-OCT-2047 21:42:11.241747309 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  1.200000000000000e+02 *km ),
         Conic.bPlaneTheta(  9.000000000000000e+01 *deg ),
         Conic.vInfinity(  3.950000000000000e-01 *km/sec ),
         Conic.inboundDec( -3.042037365115197e+01 *deg ),
         Conic.inboundRA(  9.681749059725192e+01 *deg ),
         Conic.trueAnomaly( -2.690852681843373e-06 *deg ),
         ],
      Controls = [
         [ -1.812194512473090e+05 *sec, 'Cosmic/Cosmic/ETF-08/TIME',  1.643805487526910e+05 *sec ],
         [  9.000000000000000e+01 *deg, 'Cosmic/Cosmic/ETF-08/Conic.bPlaneTheta',  1.000000000000000e+02 *deg ],
         [  1.200000000000000e+02 *km, 'Cosmic/Cosmic/ETF-08/Conic.periapsisAltitude',  1.300000000000000e+02 *km ],
         [  2.000000000000000e-01 *km/sec, 'Cosmic/Cosmic/ETF-08/Conic.vInfinity',  3.950000000000000e-01 *km/sec ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-08/Conic.inboundDec',  3.600000000000000e+02 *deg ],
         [ -3.600000000000000e+02 *deg, 'Cosmic/Cosmic/ETF-08/Conic.inboundRA',  3.600000000000000e+02 *deg ],
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
