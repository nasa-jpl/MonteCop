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

runInBatch=True                   #DO NOT MODIFY (overwritten by cosmicBatch.py calls)
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

Timeline = [
   Begin( '16-JUN-2045 23:05:53.890555221 ET' ),
   ControlPoint(
      Name = 'CP01',
      Active = True,
      Propagator = 'DIVA',
      Time = '16-JUN-2045 23:05:53.890555221 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.942586832311219e+02 *km ),
         Conic.eccentricity(  3.368218509521310e-02 ),
         Conic.inclination(  1.012143609805231e+02 *deg ),
         Conic.longitudeOfNode( -3.285652729380276e+01 *deg ),
         Conic.argumentOfPeriapsis(  2.528660122297729e+02 *deg ),
         Conic.trueAnomaly(  2.464035156555066e+02 *deg ),
         ],
      Controls = [
         [ -3.538905552210000e+02 *sec, 'Cosmic/Cosmic/CP01/TIME',  6.846109444779000e+03 *sec ],
         [  3.197673356597354e+02 *km, 'Cosmic/Cosmic/CP01/Conic.semiMajorAxis',  4.326263953043478e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP01/Conic.eccentricity', None ],
         [  8.618716692913719e+01 *deg, 'Cosmic/Cosmic/CP01/Conic.inclination',  1.053398706911677e+02 *deg ],
         [ None, 'Cosmic/Cosmic/CP01/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP01/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP01/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP01->CP02',
      Mode = 'MATCH_ELEM',
      Time = '17-JUN-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP02',
      Active = True,
      Propagator = 'DIVA',
      Time = '17-JUN-2045 23:07:24.575026889 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.855482603741045e+02 *km ),
         Conic.eccentricity(  3.057596750610494e-02 ),
         Conic.inclination(  9.867634684855528e+01 *deg ),
         Conic.longitudeOfNode(  7.137601862383184e+01 *deg ),
         Conic.argumentOfPeriapsis(  2.170794872125627e+02 *deg ),
         Conic.trueAnomaly( -2.235253468258400e+01 *deg ),
         ],
      Controls = [
         [ -4.445750268890000e+02 *sec, 'Cosmic/Cosmic/CP02/TIME',  6.755424973111000e+03 *sec ],
         [  3.237552836954470e+02 *km, 'Cosmic/Cosmic/CP02/Conic.semiMajorAxis',  4.380218544114870e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP02/Conic.eccentricity', None ],
         [  8.693123327706984e+01 *deg, 'Cosmic/Cosmic/CP02/Conic.inclination',  1.062492851164187e+02 *deg ],
         [ None, 'Cosmic/Cosmic/CP02/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP02/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP02/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP02->CP03',
      Mode = 'MATCH_ELEM',
      Time = '18-JUN-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP03',
      Active = True,
      Propagator = 'DIVA',
      Time = '18-JUN-2045 23:00:40.605647620 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.871383551135907e+02 *km ),
         Conic.eccentricity(  3.854016963517256e-02 ),
         Conic.inclination(  9.964268703189406e+01 *deg ),
         Conic.longitudeOfNode(  1.734559603844818e+02 *deg ),
         Conic.argumentOfPeriapsis(  1.760874369465258e+02 *deg ),
         Conic.trueAnomaly(  5.574457330310647e+01 *deg ),
         ],
      Controls = [
         [ -4.060564762000000e+01 *sec, 'Cosmic/Cosmic/CP03/TIME',  7.159394352380000e+03 *sec ],
         [  3.385645448922926e+02 *km, 'Cosmic/Cosmic/CP03/Conic.semiMajorAxis',  4.580579136778076e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP03/Conic.eccentricity', None ],
         [  9.194993924715762e+01 *deg, 'Cosmic/Cosmic/CP03/Conic.inclination',  1.123832590798593e+02 *deg ],
         [ None, 'Cosmic/Cosmic/CP03/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP03/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP03/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP03->CP04',
      Mode = 'MATCH_ELEM',
      Time = '19-JUN-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP04',
      Active = True,
      Propagator = 'DIVA',
      Time = '19-JUN-2045 23:05:47.304719875 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.746460168073114e+02 *km ),
         Conic.eccentricity(  1.900697387116974e-02 ),
         Conic.inclination(  9.060132800720889e+01 *deg ),
         Conic.longitudeOfNode( -9.479222537914349e+01 *deg ),
         Conic.argumentOfPeriapsis(  2.873630570972637e+02 *deg ),
         Conic.trueAnomaly( -1.183380571104932e+00 *deg ),
         ],
      Controls = [
         [ -3.473047198750000e+02 *sec, 'Cosmic/Cosmic/CP04/TIME',  6.852695280125000e+03 *sec ],
         [  3.208001014977074e+02 *km, 'Cosmic/Cosmic/CP04/Conic.semiMajorAxis',  4.340236667321923e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP04/Conic.eccentricity', None ],
         [  8.098286918727727e+01 *deg, 'Cosmic/Cosmic/CP04/Conic.inclination',  9.897906234000557e+01 *deg ],
         [ None, 'Cosmic/Cosmic/CP04/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP04/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP04/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP04->CP05',
      Mode = 'MATCH_ELEM',
      Time = '20-JUN-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP05',
      Active = True,
      Propagator = 'DIVA',
      Time = '20-JUN-2045 23:12:40.068682302 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  4.073392202932741e+02 *km ),
         Conic.eccentricity(  9.286210635162387e-02 ),
         Conic.inclination(  1.036112668496676e+02 *deg ),
         Conic.longitudeOfNode(  2.353071231230814e+00 *deg ),
         Conic.argumentOfPeriapsis(  3.557251584944241e+02 *deg ),
         Conic.trueAnomaly( -1.427680983441845e+01 *deg ),
         ],
      Controls = [
         [ -7.600686823020000e+02 *sec, 'Cosmic/Cosmic/CP05/TIME',  6.439931317698000e+03 *sec ],
         [  3.413733959936587e+02 *km, 'Cosmic/Cosmic/CP05/Conic.semiMajorAxis',  4.618581239914207e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP05/Conic.eccentricity', None ],
         [  9.264725672697743e+01 *deg, 'Cosmic/Cosmic/CP05/Conic.inclination',  1.132355359996391e+02 *deg ],
         [ None, 'Cosmic/Cosmic/CP05/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP05/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP05/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP05->CP06',
      Mode = 'MATCH_ELEM',
      Time = '21-JUN-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP06',
      Active = True,
      Propagator = 'DIVA',
      Time = '21-JUN-2045 23:18:13.307010810 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.810518261443618e+02 *km ),
         Conic.eccentricity(  3.304354697454496e-02 ),
         Conic.inclination(  9.567719486572118e+01 *deg ),
         Conic.longitudeOfNode(  1.060470734057157e+02 *deg ),
         Conic.argumentOfPeriapsis(  3.767655208163106e+01 *deg ),
         Conic.trueAnomaly( -3.501285769745293e+00 *deg ),
         ],
      Controls = [
         [ -1.093307010810000e+03 *sec, 'Cosmic/Cosmic/CP06/TIME',  6.106692989190000e+03 *sec ],
         [  3.255455565488222e+02 *km, 'Cosmic/Cosmic/CP06/Conic.semiMajorAxis',  4.404439882719359e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP06/Conic.eccentricity', None ],
         [  8.756721692319934e+01 *deg, 'Cosmic/Cosmic/CP06/Conic.inclination',  1.070265984616881e+02 *deg ],
         [ None, 'Cosmic/Cosmic/CP06/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP06/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP06/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP06->CP07',
      Mode = 'MATCH_ELEM',
      Time = '22-JUN-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP07',
      Active = True,
      Propagator = 'DIVA',
      Time = '22-JUN-2045 23:06:26.082653899 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.773064052281417e+02 *km ),
         Conic.eccentricity(  2.529248822186036e-02 ),
         Conic.inclination(  9.684970698225185e+01 *deg ),
         Conic.longitudeOfNode( -1.542518126109579e+02 *deg ),
         Conic.argumentOfPeriapsis(  3.116294792024663e+02 *deg ),
         Conic.trueAnomaly(  1.209584485611298e+02 *deg ),
         ],
      Controls = [
         [ -3.860826538990000e+02 *sec, 'Cosmic/Cosmic/CP07/TIME',  6.813917346101000e+03 *sec ],
         [  3.238238070231065e+02 *km, 'Cosmic/Cosmic/CP07/Conic.semiMajorAxis',  4.381145624430264e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP07/Conic.eccentricity', None ],
         [  8.710544158689844e+01 *deg, 'Cosmic/Cosmic/CP07/Conic.inclination',  1.064622063839870e+02 *deg ],
         [ None, 'Cosmic/Cosmic/CP07/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP07/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP07/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP07->CP08',
      Mode = 'MATCH_ELEM',
      Time = '23-JUN-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP08',
      Active = True,
      Propagator = 'DIVA',
      Time = '23-JUN-2045 23:04:20.954958951 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.756625344110731e+02 *km ),
         Conic.eccentricity(  7.088172590873492e-03 ),
         Conic.inclination(  9.311174332420546e+01 *deg ),
         Conic.longitudeOfNode( -6.060627674621874e+01 *deg ),
         Conic.argumentOfPeriapsis(  1.637688795404580e+02 *deg ),
         Conic.trueAnomaly( -5.617729117293970e+01 *deg ),
         ],
      Controls = [
         [ -2.609549589510000e+02 *sec, 'Cosmic/Cosmic/CP08/TIME',  6.939045041049000e+03 *sec ],
         [  3.234023007423810e+02 *km, 'Cosmic/Cosmic/CP08/Conic.semiMajorAxis',  4.375442892396920e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP08/Conic.eccentricity', None ],
         [  8.390501355697792e+01 *deg, 'Cosmic/Cosmic/CP08/Conic.inclination',  1.025505721251952e+02 *deg ],
         [ None, 'Cosmic/Cosmic/CP08/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP08/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP08/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP08->CP09',
      Mode = 'MATCH_ELEM',
      Time = '24-JUN-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP09',
      Active = True,
      Propagator = 'DIVA',
      Time = '24-JUN-2045 23:23:45.826367377 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  4.002136707838361e+02 *km ),
         Conic.eccentricity(  7.712471797569620e-02 ),
         Conic.inclination(  1.019784969079959e+02 *deg ),
         Conic.longitudeOfNode(  3.788598914137624e+01 *deg ),
         Conic.argumentOfPeriapsis(  1.787398333202362e+02 *deg ),
         Conic.trueAnomaly(  1.058587181745952e+01 *deg ),
         ],
      Controls = [
         [ -1.425826367377000e+03 *sec, 'Cosmic/Cosmic/CP09/TIME',  5.774173632623000e+03 *sec ],
         [  3.384827803945301e+02 *km, 'Cosmic/Cosmic/CP09/Conic.semiMajorAxis',  4.579472911220113e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP09/Conic.eccentricity', None ],
         [  9.141192917787066e+01 *deg, 'Cosmic/Cosmic/CP09/Conic.inclination',  1.117256912173975e+02 *deg ],
         [ None, 'Cosmic/Cosmic/CP09/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP09/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP09/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP09->CP10',
      Mode = 'MATCH_ELEM',
      Time = '25-JUN-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP10',
      Active = True,
      Propagator = 'DIVA',
      Time = '25-JUN-2045 23:26:34.540376527 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.822710789931491e+02 *km ),
         Conic.eccentricity(  3.199074362079189e-02 ),
         Conic.inclination(  9.631346609450016e+01 *deg ),
         Conic.longitudeOfNode(  1.389708583849721e+02 *deg ),
         Conic.argumentOfPeriapsis(  1.930139521065696e+02 *deg ),
         Conic.trueAnomaly(  3.794314975264071e+01 *deg ),
         ],
      Controls = [
         [ -1.594540376527000e+03 *sec, 'Cosmic/Cosmic/CP10/TIME',  5.605459623473000e+03 *sec ],
         [  3.325788077946532e+02 *km, 'Cosmic/Cosmic/CP10/Conic.semiMajorAxis',  4.499595634868838e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP10/Conic.eccentricity', None ],
         [  8.926687573708443e+01 *deg, 'Cosmic/Cosmic/CP10/Conic.inclination',  1.091039592342143e+02 *deg ],
         [ None, 'Cosmic/Cosmic/CP10/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP10/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP10/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP10->CP11',
      Mode = 'MATCH_ELEM',
      Time = '26-JUN-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP11',
      Active = True,
      Propagator = 'DIVA',
      Time = '26-JUN-2045 23:11:37.703660804 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.740539754439043e+02 *km ),
         Conic.eccentricity(  1.000000000000000e-05 ),
         Conic.inclination(  9.223572228973720e+01 *deg ),
         Conic.longitudeOfNode( -1.225176169532189e+02 *deg ),
         Conic.argumentOfPeriapsis(  2.453554841424897e+02 *deg ),
         Conic.trueAnomaly(  2.084273681511839e+01 *deg ),
         ],
      Controls = [
         [ -6.977036608040000e+02 *sec, 'Cosmic/Cosmic/CP11/TIME',  6.502296339196000e+03 *sec ],
         [  3.222691909833110e+02 *km, 'Cosmic/Cosmic/CP11/Conic.semiMajorAxis',  4.360112583891854e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP11/Conic.eccentricity', None ],
         [  8.350583684917571e+01 *deg, 'Cosmic/Cosmic/CP11/Conic.inclination',  1.020626894823259e+02 *deg ],
         [ None, 'Cosmic/Cosmic/CP11/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP11/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP11/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP11->CP12',
      Mode = 'MATCH_ELEM',
      Time = '27-JUN-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP12',
      Active = True,
      Propagator = 'DIVA',
      Time = '27-JUN-2045 23:17:05.876993844 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.882281597523331e+02 *km ),
         Conic.eccentricity(  4.492005335300961e-02 ),
         Conic.inclination(  9.944806361577976e+01 *deg ),
         Conic.longitudeOfNode( -2.988505908549115e+01 *deg ),
         Conic.argumentOfPeriapsis(  5.683494017942817e+00 *deg ),
         Conic.trueAnomaly( -5.603669629126096e+01 *deg ),
         ],
      Controls = [
         [ -1.025876993844000e+03 *sec, 'Cosmic/Cosmic/CP12/TIME',  6.174123006156000e+03 *sec ],
         [  3.283766725447578e+02 *km, 'Cosmic/Cosmic/CP12/Conic.semiMajorAxis',  4.442743216782016e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP12/Conic.eccentricity', None ],
         [  8.805697765448532e+01 *deg, 'Cosmic/Cosmic/CP12/Conic.inclination',  1.076251949110376e+02 *deg ],
         [ None, 'Cosmic/Cosmic/CP12/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP12/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP12/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP12->CP13',
      Mode = 'MATCH_ELEM',
      Time = '28-JUN-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP13',
      Active = True,
      Propagator = 'DIVA',
      Time = '28-JUN-2045 23:23:18.218404462 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.851949959822089e+02 *km ),
         Conic.eccentricity(  4.007079248179419e-02 ),
         Conic.inclination(  9.841584877223625e+01 *deg ),
         Conic.longitudeOfNode(  7.304437486549931e+01 *deg ),
         Conic.argumentOfPeriapsis(  3.640556840557826e+02 *deg ),
         Conic.trueAnomaly(  9.091883461053097e+00 *deg ),
         ],
      Controls = [
         [ -1.398218404462000e+03 *sec, 'Cosmic/Cosmic/CP13/TIME',  5.801781595538000e+03 *sec ],
         [  3.277013786595105e+02 *km, 'Cosmic/Cosmic/CP13/Conic.semiMajorAxis',  4.433606887746317e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP13/Conic.eccentricity', None ],
         [  8.897509045124687e+01 *deg, 'Cosmic/Cosmic/CP13/Conic.inclination',  1.087473327737462e+02 *deg ],
         [ None, 'Cosmic/Cosmic/CP13/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP13/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP13/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP13->CP14',
      Mode = 'MATCH_ELEM',
      Time = '29-JUN-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP14',
      Active = True,
      Propagator = 'DIVA',
      Time = '29-JUN-2045 23:19:55.995084311 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.867719732527812e+02 *km ),
         Conic.eccentricity(  4.323377929142732e-02 ),
         Conic.inclination(  9.958344134524580e+01 *deg ),
         Conic.longitudeOfNode(  1.746693379729472e+02 *deg ),
         Conic.argumentOfPeriapsis( -8.262374811825808e+00 *deg ),
         Conic.trueAnomaly(  6.115247999111143e+01 *deg ),
         ],
      Controls = [
         [ -1.195995084311000e+03 *sec, 'Cosmic/Cosmic/CP14/TIME',  6.004004915689000e+03 *sec ],
         [  3.357420986537114e+02 *km, 'Cosmic/Cosmic/CP14/Conic.semiMajorAxis',  4.542393099432566e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP14/Conic.eccentricity', None ],
         [  9.065006938931516e+01 *deg, 'Cosmic/Cosmic/CP14/Conic.inclination',  1.107945292536074e+02 *deg ],
         [ None, 'Cosmic/Cosmic/CP14/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP14/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP14/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP14->CP15',
      Mode = 'MATCH_ELEM',
      Time = '30-JUN-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP15',
      Active = True,
      Propagator = 'DIVA',
      Time = '30-JUN-2045 23:25:16.383050781 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.747224489330362e+02 *km ),
         Conic.eccentricity(  1.781054024946509e-02 ),
         Conic.inclination(  9.058080420745652e+01 *deg ),
         Conic.longitudeOfNode( -9.356106630813473e+01 *deg ),
         Conic.argumentOfPeriapsis(  9.562135793413105e+01 *deg ),
         Conic.trueAnomaly(  1.022412530950818e+01 *deg ),
         ],
      Controls = [
         [ -1.516383050781000e+03 *sec, 'Cosmic/Cosmic/CP15/TIME',  5.683616949219000e+03 *sec ],
         [  3.221444120566615e+02 *km, 'Cosmic/Cosmic/CP15/Conic.semiMajorAxis',  4.358424398413655e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP15/Conic.eccentricity', None ],
         [  8.102172183402739e+01 *deg, 'Cosmic/Cosmic/CP15/Conic.inclination',  9.902654890825569e+01 *deg ],
         [ None, 'Cosmic/Cosmic/CP15/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP15/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP15/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP15->CP16',
      Mode = 'MATCH_ELEM',
      Time = '01-JUL-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP16',
      Active = True,
      Propagator = 'DIVA',
      Time = '01-JUL-2045 23:22:59.987362975 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  4.003492774552055e+02 *km ),
         Conic.eccentricity(  7.377856021586114e-02 ),
         Conic.inclination(  1.023098613371765e+02 *deg ),
         Conic.longitudeOfNode(  4.030960152495059e+00 *deg ),
         Conic.argumentOfPeriapsis(  1.761216782739830e+02 *deg ),
         Conic.trueAnomaly( -2.790680052520652e+01 *deg ),
         ],
      Controls = [
         [ -1.379987362975000e+03 *sec, 'Cosmic/Cosmic/CP16/TIME',  5.820012637025000e+03 *sec ],
         [  3.367962920737132e+02 *km, 'Cosmic/Cosmic/CP16/Conic.semiMajorAxis',  4.556655716291413e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP16/Conic.eccentricity', None ],
         [  9.125380561439783e+01 *deg, 'Cosmic/Cosmic/CP16/Conic.inclination',  1.115324290842640e+02 *deg ],
         [ None, 'Cosmic/Cosmic/CP16/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP16/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP16/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP16->CP17',
      Mode = 'MATCH_ELEM',
      Time = '02-JUL-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP17',
      Active = True,
      Propagator = 'DIVA',
      Time = '02-JUL-2045 23:34:16.107177482 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.825530125666999e+02 *km ),
         Conic.eccentricity(  3.597540623229321e-02 ),
         Conic.inclination(  9.662207310534792e+01 *deg ),
         Conic.longitudeOfNode(  1.074873615028808e+02 *deg ),
         Conic.argumentOfPeriapsis(  2.128470370080220e+02 *deg ),
         Conic.trueAnomaly( -5.496347073300691e+00 *deg ),
         ],
      Controls = [
         [ -2.056107177482000e+03 *sec, 'Cosmic/Cosmic/CP17/TIME',  5.143892822518000e+03 *sec ],
         [  3.243054405620119e+02 *km, 'Cosmic/Cosmic/CP17/Conic.semiMajorAxis',  4.387661842897809e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP17/Conic.eccentricity', None ],
         [  8.797885479515001e+01 *deg, 'Cosmic/Cosmic/CP17/Conic.inclination',  1.075297114162945e+02 *deg ],
         [ None, 'Cosmic/Cosmic/CP17/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP17/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP17/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP17->CP18',
      Mode = 'MATCH_ELEM',
      Time = '03-JUL-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP18',
      Active = True,
      Propagator = 'DIVA',
      Time = '03-JUL-2045 23:23:57.581867975 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.792754561981763e+02 *km ),
         Conic.eccentricity(  2.655663653502589e-02 ),
         Conic.inclination(  9.722928976562744e+01 *deg ),
         Conic.longitudeOfNode( -1.521914271392588e+02 *deg ),
         Conic.argumentOfPeriapsis(  1.449969211780183e+02 *deg ),
         Conic.trueAnomaly(  1.027982551737514e+02 *deg ),
         ],
      Controls = [
         [ -1.437581867975000e+03 *sec, 'Cosmic/Cosmic/CP18/TIME',  5.762418132025000e+03 *sec ],
         [  3.245366518924512e+02 *km, 'Cosmic/Cosmic/CP18/Conic.semiMajorAxis',  4.390789996191986e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP18/Conic.eccentricity', None ],
         [  8.774043247509245e+01 *deg, 'Cosmic/Cosmic/CP18/Conic.inclination',  1.072383063584463e+02 *deg ],
         [ None, 'Cosmic/Cosmic/CP18/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP18/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP18/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP18->CP19',
      Mode = 'MATCH_ELEM',
      Time = '04-JUL-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP19',
      Active = True,
      Propagator = 'DIVA',
      Time = '04-JUL-2045 23:28:48.186116395 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.763353915996310e+02 *km ),
         Conic.eccentricity(  8.873827607936056e-03 ),
         Conic.inclination(  9.344710427435130e+01 *deg ),
         Conic.longitudeOfNode( -6.060357598872034e+01 *deg ),
         Conic.argumentOfPeriapsis(  3.376866150329296e+02 *deg ),
         Conic.trueAnomaly( -4.712350617960661e+01 *deg ),
         ],
      Controls = [
         [ -1.728186116395000e+03 *sec, 'Cosmic/Cosmic/CP19/TIME',  5.471813883605000e+03 *sec ],
         [  3.228305892581812e+02 *km, 'Cosmic/Cosmic/CP19/Conic.semiMajorAxis',  4.367707972316568e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP19/Conic.eccentricity', None ],
         [  8.361474884758911e+01 *deg, 'Cosmic/Cosmic/CP19/Conic.inclination',  1.021958041470534e+02 *deg ],
         [ None, 'Cosmic/Cosmic/CP19/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP19/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP19/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP19->CP20',
      Mode = 'MATCH_ELEM',
      Time = '05-JUL-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP20',
      Active = True,
      Propagator = 'DIVA',
      Time = '05-JUL-2045 23:45:18.513333865 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  4.001200268307048e+02 *km ),
         Conic.eccentricity(  7.578513028736762e-02 ),
         Conic.inclination(  1.019383367754399e+02 *deg ),
         Conic.longitudeOfNode(  3.850075212863845e+01 *deg ),
         Conic.argumentOfPeriapsis(  3.586621373947089e+02 *deg ),
         Conic.trueAnomaly(  9.116908661339602e+00 *deg ),
         ],
      Controls = [
         [ -2.718513333865000e+03 *sec, 'Cosmic/Cosmic/CP20/TIME',  4.481486666135000e+03 *sec ],
         [  3.405904168058466e+02 *km, 'Cosmic/Cosmic/CP20/Conic.semiMajorAxis',  4.607987992079101e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP20/Conic.eccentricity', None ],
         [  9.220126833331287e+01 *deg, 'Cosmic/Cosmic/CP20/Conic.inclination',  1.126904390740491e+02 *deg ],
         [ None, 'Cosmic/Cosmic/CP20/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP20/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP20/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP20->CP21',
      Mode = 'MATCH_ELEM',
      Time = '06-JUL-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP21',
      Active = True,
      Propagator = 'DIVA',
      Time = '06-JUL-2045 23:46:00.755356997 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.844209442950950e+02 *km ),
         Conic.eccentricity(  3.620891446607521e-02 ),
         Conic.inclination(  9.709938041884003e+01 *deg ),
         Conic.longitudeOfNode(  1.403245179294436e+02 *deg ),
         Conic.argumentOfPeriapsis(  1.573631379665908e+01 *deg ),
         Conic.trueAnomaly(  3.050988778849241e+01 *deg ),
         ],
      Controls = [
         [ -2.760755356997000e+03 *sec, 'Cosmic/Cosmic/CP21/TIME',  4.439244643003000e+03 *sec ],
         [  3.280145438038667e+02 *km, 'Cosmic/Cosmic/CP21/Conic.semiMajorAxis',  4.437843827934667e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP21/Conic.eccentricity', None ],
         [  8.729722751010860e+01 *deg, 'Cosmic/Cosmic/CP21/Conic.inclination',  1.066966114012438e+02 *deg ],
         [ None, 'Cosmic/Cosmic/CP21/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP21/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP21/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP21->CP22',
      Mode = 'MATCH_ELEM',
      Time = '07-JUL-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP22',
      Active = True,
      Propagator = 'DIVA',
      Time = '07-JUL-2045 23:51:36.267809423 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.753315317406874e+02 *km ),
         Conic.eccentricity(  1.039502835280058e-02 ),
         Conic.inclination(  9.288643648073433e+01 *deg ),
         Conic.longitudeOfNode( -1.279797303048334e+02 *deg ),
         Conic.argumentOfPeriapsis(  1.620346502840061e+02 *deg ),
         Conic.trueAnomaly( -5.472725668701469e+01 *deg ),
         ],
      Controls = [
         [ -3.096267809423000e+03 *sec, 'Cosmic/Cosmic/CP22/TIME',  4.103732190577000e+03 *sec ],
         [  3.214053839740200e+02 *km, 'Cosmic/Cosmic/CP22/Conic.semiMajorAxis',  4.348425783177917e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP22/Conic.eccentricity', None ],
         [  8.398641107744410e+01 *deg, 'Cosmic/Cosmic/CP22/Conic.inclination',  1.026500579835428e+02 *deg ],
         [ None, 'Cosmic/Cosmic/CP22/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP22/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP22/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP22->CP23',
      Mode = 'MATCH_ELEM',
      Time = '08-JUL-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP23',
      Active = True,
      Propagator = 'DIVA',
      Time = '08-JUL-2045 23:52:53.196519015 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.967106248792879e+02 *km ),
         Conic.eccentricity(  6.509186967853918e-02 ),
         Conic.inclination(  1.013015189065963e+02 *deg ),
         Conic.longitudeOfNode( -3.109710959173560e+01 *deg ),
         Conic.argumentOfPeriapsis(  1.765409426917757e+02 *deg ),
         Conic.trueAnomaly( -3.096475486018146e+01 *deg ),
         ],
      Controls = [
         [ -3.173196519015000e+03 *sec, 'Cosmic/Cosmic/CP23/TIME',  4.026803480985000e+03 *sec ],
         [  3.342761552328238e+02 *km, 'Cosmic/Cosmic/CP23/Conic.semiMajorAxis',  4.522559747267616e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP23/Conic.eccentricity', None ],
         [  9.011444845439897e+01 *deg, 'Cosmic/Cosmic/CP23/Conic.inclination',  1.101398814442654e+02 *deg ],
         [ None, 'Cosmic/Cosmic/CP23/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP23/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP23/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP23->CP24',
      Mode = 'MATCH_ELEM',
      Time = '09-JUL-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP24',
      Active = True,
      Propagator = 'DIVA',
      Time = '10-JUL-2045 00:08:31.238109987 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.813468045177559e+02 *km ),
         Conic.eccentricity(  2.623729670675605e-02 ),
         Conic.inclination(  9.586968535956927e+01 *deg ),
         Conic.longitudeOfNode(  7.202593421200645e+01 *deg ),
         Conic.argumentOfPeriapsis(  2.089813631438251e+02 *deg ),
         Conic.trueAnomaly(  1.237902062736203e+01 *deg ),
         ],
      Controls = [
         [ -4.111238109987000e+03 *sec, 'Cosmic/Cosmic/CP24/TIME',  3.088761890013000e+03 *sec ],
         [  3.275993492111934e+02 *km, 'Cosmic/Cosmic/CP24/Conic.semiMajorAxis',  4.432226489327911e+02 *km ],
         [  1.000000000000000e-05, 'Cosmic/Cosmic/CP24/Conic.eccentricity', None ],
         [  8.844384133937027e+01 *deg, 'Cosmic/Cosmic/CP24/Conic.inclination',  1.080980283036748e+02 *deg ],
         [ None, 'Cosmic/Cosmic/CP24/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP24/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP24/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP24->CP25',
      Mode = 'MATCH_ELEM',
      Time = '10-JUL-2045 12:00:00.000000000 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-03 *km,
      VelTol =  1.000000000000000e-06 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP25',
      Active = True,
      Propagator = 'DIVA',
      Time = '10-JUL-2045 23:56:47.103736421 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Cartesian.x(  1.310347074294615e+02 *km ),
         Cartesian.y( -6.783710445626366e+01 *km ),
         Cartesian.z( -3.553686382535691e+02 *km ),
         Cartesian.dx( -1.248359941535009e-01 *km/sec ),
         Cartesian.dy(  1.147667843457802e-02 *km/sec ),
         Cartesian.dz( -5.195031253923422e-02 *km/sec ),
         ],
      Controls = [
         [  1.302627735635875e+02 *km, 'Cosmic/Cosmic/CP25/X',  1.482627735635875e+02 *km ],
         [ -8.699746462514287e+01 *km, 'Cosmic/Cosmic/CP25/Y', -6.699746462514287e+01 *km ],
         [ -3.603169810278590e+02 *km, 'Cosmic/Cosmic/CP25/Z', -3.403169810278590e+02 *km ],
         [ -1.436041459403040e-01 *km/sec, 'Cosmic/Cosmic/CP25/DX', -1.036041459403040e-01 *km/sec ],
         [  1.147274737082190e-02 *km/sec, 'Cosmic/Cosmic/CP25/DY',  4.447274737082191e-02 *km/sec ],
         [ -5.650363465222682e-02 *km/sec, 'Cosmic/Cosmic/CP25/DZ', -3.550363465222682e-02 *km/sec ],
         ],
      ),
   End( '10-JUL-2045 23:56:47.103736421 ET' ),
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
