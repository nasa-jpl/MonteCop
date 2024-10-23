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

PosTolIn     = 1.0E-6
VelTolIn     = 1.0E-9
MajorFeasTol = 1e-8
MajorOptTol  = 1e-3

#--------------------------------------------------------------------------
# Setup SNOPT params for cosmicBatch() run:
#--------------------------------------------------------------------------
MaxIterLimit   = 300                 #Major iterations limit

#set SNOPT params (these values are overwritten on cosmicBatch.py calls):
MajorStepLimit = 0.9
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

Timeline = [
   Begin( '16-JUN-2045 23:05:53.890555221 ET' ),
   ControlPoint(
      Name = 'CP01',
      Active = True,
      Propagator = 'DIVA',
      Time = '16-JUN-2045 23:03:21.200725430 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.944113343132354e+02 *km ),
         Conic.eccentricity(  3.457825549134075e-02 ),
         Conic.inclination(  1.012321263062204e+02 *deg ),
         Conic.longitudeOfNode( -3.146968383377891e+01 *deg ),
         Conic.argumentOfPeriapsis(  2.525900494285806e+02 *deg ),
         Conic.trueAnomaly(  2.464040211404262e+02 *deg ),
         ],
      Controls = [
         [ -2.012007254300000e+02 *sec, 'Cosmic/Cosmic/CP01/TIME',  6.998799274570000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP02',
      Active = True,
      Propagator = 'DIVA',
      Time = '17-JUN-2045 23:05:05.334962557 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.852350900835588e+02 *km ),
         Conic.eccentricity(  2.966576594910706e-02 ),
         Conic.inclination(  9.857467469882933e+01 *deg ),
         Conic.longitudeOfNode(  7.260124723211128e+01 *deg ),
         Conic.argumentOfPeriapsis(  2.172345154539130e+02 *deg ),
         Conic.trueAnomaly( -2.344585512960203e+01 *deg ),
         ],
      Controls = [
         [ -3.053349625570000e+02 *sec, 'Cosmic/Cosmic/CP02/TIME',  6.894665037443000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP03',
      Active = True,
      Propagator = 'DIVA',
      Time = '18-JUN-2045 23:00:00.000000000 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.867421736763657e+02 *km ),
         Conic.eccentricity(  3.762268678574027e-02 ),
         Conic.inclination(  9.953433503951430e+01 *deg ),
         Conic.longitudeOfNode(  1.744335942888573e+02 *deg ),
         Conic.argumentOfPeriapsis(  1.743277596359993e+02 *deg ),
         Conic.trueAnomaly(  5.830402584731463e+01 *deg ),
         ],
      Controls = [
         [  0.000000000000000e+00 *sec, 'Cosmic/Cosmic/CP03/TIME',  7.200000000000000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP04',
      Active = True,
      Propagator = 'DIVA',
      Time = '19-JUN-2045 23:05:45.183779052 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.746741903384587e+02 *km ),
         Conic.eccentricity(  1.882649621291443e-02 ),
         Conic.inclination(  9.059294450585527e+01 *deg ),
         Conic.longitudeOfNode( -9.404144335767884e+01 *deg ),
         Conic.argumentOfPeriapsis(  2.873600394233990e+02 *deg ),
         Conic.trueAnomaly( -4.964877269450681e-01 *deg ),
         ],
      Controls = [
         [ -3.451837790520000e+02 *sec, 'Cosmic/Cosmic/CP04/TIME',  6.854816220948000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP05',
      Active = True,
      Propagator = 'DIVA',
      Time = '20-JUN-2045 23:12:24.015753466 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  4.072187437648613e+02 *km ),
         Conic.eccentricity(  9.242544927866743e-02 ),
         Conic.inclination(  1.035318039817752e+02 *deg ),
         Conic.longitudeOfNode(  3.193927289991667e+00 *deg ),
         Conic.argumentOfPeriapsis(  3.556865250059313e+02 *deg ),
         Conic.trueAnomaly( -1.417640283470790e+01 *deg ),
         ],
      Controls = [
         [ -7.440157534660000e+02 *sec, 'Cosmic/Cosmic/CP05/TIME',  6.455984246534000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP06',
      Active = True,
      Propagator = 'DIVA',
      Time = '21-JUN-2045 23:18:35.281215740 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.810880239863536e+02 *km ),
         Conic.eccentricity(  3.299264031847432e-02 ),
         Conic.inclination(  9.563229312689698e+01 *deg ),
         Conic.longitudeOfNode(  1.067539834374459e+02 *deg ),
         Conic.argumentOfPeriapsis(  3.763315244002453e+01 *deg ),
         Conic.trueAnomaly( -3.364314124893465e+00 *deg ),
         ],
      Controls = [
         [ -1.115281215740000e+03 *sec, 'Cosmic/Cosmic/CP06/TIME',  6.084718784260000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP07',
      Active = True,
      Propagator = 'DIVA',
      Time = '22-JUN-2045 23:06:50.382546676 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.772759656655941e+02 *km ),
         Conic.eccentricity(  2.517427881490700e-02 ),
         Conic.inclination(  9.675091502250356e+01 *deg ),
         Conic.longitudeOfNode( -1.535813288802330e+02 *deg ),
         Conic.argumentOfPeriapsis(  3.115405933133050e+02 *deg ),
         Conic.trueAnomaly(  1.209168361915043e+02 *deg ),
         ],
      Controls = [
         [ -4.103825466760000e+02 *sec, 'Cosmic/Cosmic/CP07/TIME',  6.789617453324000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP08',
      Active = True,
      Propagator = 'DIVA',
      Time = '23-JUN-2045 23:07:14.208560867 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.760863478427415e+02 *km ),
         Conic.eccentricity(  8.553914039779035e-03 ),
         Conic.inclination(  9.332465565959443e+01 *deg ),
         Conic.longitudeOfNode( -6.069963855294673e+01 *deg ),
         Conic.argumentOfPeriapsis(  1.656092430129327e+02 *deg ),
         Conic.trueAnomaly( -5.550000506326776e+01 *deg ),
         ],
      Controls = [
         [ -4.342085608670000e+02 *sec, 'Cosmic/Cosmic/CP08/TIME',  6.765791439133000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP09',
      Active = True,
      Propagator = 'DIVA',
      Time = '24-JUN-2045 23:24:31.707070372 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.999315280434796e+02 *km ),
         Conic.eccentricity(  7.642711657271620e-02 ),
         Conic.inclination(  1.018772547968503e+02 *deg ),
         Conic.longitudeOfNode(  3.833460588689411e+01 *deg ),
         Conic.argumentOfPeriapsis(  1.784499542604011e+02 *deg ),
         Conic.trueAnomaly(  1.043493005801905e+01 *deg ),
         ],
      Controls = [
         [ -1.471707070372000e+03 *sec, 'Cosmic/Cosmic/CP09/TIME',  5.728292929628000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP10',
      Active = True,
      Propagator = 'DIVA',
      Time = '25-JUN-2045 23:27:21.167739673 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.825374030800735e+02 *km ),
         Conic.eccentricity(  3.267881689215893e-02 ),
         Conic.inclination(  9.639753729766424e+01 *deg ),
         Conic.longitudeOfNode(  1.394761498788563e+02 *deg ),
         Conic.argumentOfPeriapsis(  1.925764085852427e+02 *deg ),
         Conic.trueAnomaly(  3.772289237990044e+01 *deg ),
         ],
      Controls = [
         [ -1.641167739673000e+03 *sec, 'Cosmic/Cosmic/CP10/TIME',  5.558832260327000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP11',
      Active = True,
      Propagator = 'DIVA',
      Time = '26-JUN-2045 23:12:52.988813868 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.739740526326042e+02 *km ),
         Conic.eccentricity(  9.993488253658916e-06 ),
         Conic.inclination(  9.214487626722062e+01 *deg ),
         Conic.longitudeOfNode( -1.221323112370238e+02 *deg ),
         Conic.argumentOfPeriapsis(  2.452130769547535e+02 *deg ),
         Conic.trueAnomaly(  2.071636552942160e+01 *deg ),
         ],
      Controls = [
         [ -7.729888138680000e+02 *sec, 'Cosmic/Cosmic/CP11/TIME',  6.427011186132000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP12',
      Active = True,
      Propagator = 'DIVA',
      Time = '27-JUN-2045 23:18:45.155681545 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.882418909696012e+02 *km ),
         Conic.eccentricity(  4.507845573977087e-02 ),
         Conic.inclination(  9.943150751922387e+01 *deg ),
         Conic.longitudeOfNode( -2.963866948235722e+01 *deg ),
         Conic.argumentOfPeriapsis(  5.927383764220171e+00 *deg ),
         Conic.trueAnomaly( -5.610568018843722e+01 *deg ),
         ],
      Controls = [
         [ -1.125155681545000e+03 *sec, 'Cosmic/Cosmic/CP12/TIME',  6.074844318455000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP13',
      Active = True,
      Propagator = 'DIVA',
      Time = '28-JUN-2045 23:24:33.942153503 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.850368302567767e+02 *km ),
         Conic.eccentricity(  3.970762157707873e-02 ),
         Conic.inclination(  9.836026393814537e+01 *deg ),
         Conic.longitudeOfNode(  7.327038244714316e+01 *deg ),
         Conic.argumentOfPeriapsis(  3.636790031203886e+02 *deg ),
         Conic.trueAnomaly(  9.029076234668038e+00 *deg ),
         ],
      Controls = [
         [ -1.473942153503000e+03 *sec, 'Cosmic/Cosmic/CP13/TIME',  5.726057846497000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP14',
      Active = True,
      Propagator = 'DIVA',
      Time = '29-JUN-2045 23:21:16.119214904 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.868080251110255e+02 *km ),
         Conic.eccentricity(  4.346099406986192e-02 ),
         Conic.inclination(  9.956500292832223e+01 *deg ),
         Conic.longitudeOfNode(  1.749117826567594e+02 *deg ),
         Conic.argumentOfPeriapsis( -8.322034722130443e+00 *deg ),
         Conic.trueAnomaly(  6.098699603091532e+01 *deg ),
         ],
      Controls = [
         [ -1.276119214904000e+03 *sec, 'Cosmic/Cosmic/CP14/TIME',  5.923880785096000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP15',
      Active = True,
      Propagator = 'DIVA',
      Time = '30-JUN-2045 23:26:39.415615272 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.745826090166839e+02 *km ),
         Conic.eccentricity(  1.751371791440925e-02 ),
         Conic.inclination(  9.051252228536652e+01 *deg ),
         Conic.longitudeOfNode( -9.336046690918991e+01 *deg ),
         Conic.argumentOfPeriapsis(  9.520454900537065e+01 *deg ),
         Conic.trueAnomaly(  1.051320907945064e+01 *deg ),
         ],
      Controls = [
         [ -1.599415615272000e+03 *sec, 'Cosmic/Cosmic/CP15/TIME',  5.600584384728000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP16',
      Active = True,
      Propagator = 'DIVA',
      Time = '01-JUL-2045 23:24:11.392784454 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  4.000702043797144e+02 *km ),
         Conic.eccentricity(  7.322670060095804e-02 ),
         Conic.inclination(  1.022248356633796e+02 *deg ),
         Conic.longitudeOfNode(  4.205978977370402e+00 *deg ),
         Conic.argumentOfPeriapsis(  1.762155578957025e+02 *deg ),
         Conic.trueAnomaly( -2.815230430060966e+01 *deg ),
         ],
      Controls = [
         [ -1.451392784454000e+03 *sec, 'Cosmic/Cosmic/CP16/TIME',  5.748607215546000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP17',
      Active = True,
      Propagator = 'DIVA',
      Time = '02-JUL-2045 23:35:09.072677964 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.824981409484766e+02 *km ),
         Conic.eccentricity(  3.591018798441098e-02 ),
         Conic.inclination(  9.661666886130539e+01 *deg ),
         Conic.longitudeOfNode(  1.076706156472518e+02 *deg ),
         Conic.argumentOfPeriapsis(  2.125392493468430e+02 *deg ),
         Conic.trueAnomaly( -5.559282444884325e+00 *deg ),
         ],
      Controls = [
         [ -2.109072677964000e+03 *sec, 'Cosmic/Cosmic/CP17/TIME',  5.090927322036000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP18',
      Active = True,
      Propagator = 'DIVA',
      Time = '03-JUL-2045 23:24:41.133124468 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.791854701269284e+02 *km ),
         Conic.eccentricity(  2.660914644220505e-02 ),
         Conic.inclination(  9.717704156490531e+01 *deg ),
         Conic.longitudeOfNode( -1.519734767728459e+02 *deg ),
         Conic.argumentOfPeriapsis(  1.449000918196643e+02 *deg ),
         Conic.trueAnomaly(  1.026830052122992e+02 *deg ),
         ],
      Controls = [
         [ -1.481133124468000e+03 *sec, 'Cosmic/Cosmic/CP18/TIME',  5.718866875532000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP19',
      Active = True,
      Propagator = 'DIVA',
      Time = '04-JUL-2045 23:30:38.872850975 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.764843327036810e+02 *km ),
         Conic.eccentricity(  9.722474694396433e-03 ),
         Conic.inclination(  9.355132009246303e+01 *deg ),
         Conic.longitudeOfNode( -6.078207757262380e+01 *deg ),
         Conic.argumentOfPeriapsis(  3.391446842827818e+02 *deg ),
         Conic.trueAnomaly( -4.704265092188786e+01 *deg ),
         ],
      Controls = [
         [ -1.838872850975000e+03 *sec, 'Cosmic/Cosmic/CP19/TIME',  5.361127149025000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP20',
      Active = True,
      Propagator = 'DIVA',
      Time = '05-JUL-2045 23:45:28.772902915 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.998298609994976e+02 *km ),
         Conic.eccentricity(  7.528414617138635e-02 ),
         Conic.inclination(  1.018583144165152e+02 *deg ),
         Conic.longitudeOfNode(  3.870267407817597e+01 *deg ),
         Conic.argumentOfPeriapsis(  3.585267322599039e+02 *deg ),
         Conic.trueAnomaly(  9.078687768898853e+00 *deg ),
         ],
      Controls = [
         [ -2.728772902915000e+03 *sec, 'Cosmic/Cosmic/CP20/TIME',  4.471227097085000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP21',
      Active = True,
      Propagator = 'DIVA',
      Time = '06-JUL-2045 23:45:51.019845284 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.844161389826264e+02 *km ),
         Conic.eccentricity(  3.644487916704237e-02 ),
         Conic.inclination(  9.710231584955183e+01 *deg ),
         Conic.longitudeOfNode(  1.405858994192565e+02 *deg ),
         Conic.argumentOfPeriapsis(  1.545168205464789e+01 *deg ),
         Conic.trueAnomaly(  3.056550743073311e+01 *deg ),
         ],
      Controls = [
         [ -2.751019845284000e+03 *sec, 'Cosmic/Cosmic/CP21/TIME',  4.448980154716000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP22',
      Active = True,
      Propagator = 'DIVA',
      Time = '07-JUL-2045 23:50:59.078744952 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.751004814899285e+02 *km ),
         Conic.eccentricity(  9.836578276642950e-03 ),
         Conic.inclination(  9.276182995537015e+01 *deg ),
         Conic.longitudeOfNode( -1.276199482601473e+02 *deg ),
         Conic.argumentOfPeriapsis(  1.614851138413746e+02 *deg ),
         Conic.trueAnomaly( -5.465344121364425e+01 *deg ),
         ],
      Controls = [
         [ -3.059078744952000e+03 *sec, 'Cosmic/Cosmic/CP22/TIME',  4.140921255048000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP23',
      Active = True,
      Propagator = 'DIVA',
      Time = '08-JUL-2045 23:52:10.665812444 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.965264143099653e+02 *km ),
         Conic.eccentricity(  6.402916889297543e-02 ),
         Conic.inclination(  1.012511766462285e+02 *deg ),
         Conic.longitudeOfNode( -3.082362303743856e+01 *deg ),
         Conic.argumentOfPeriapsis(  1.760227194605311e+02 *deg ),
         Conic.trueAnomaly( -3.062434104147649e+01 *deg ),
         ],
      Controls = [
         [ -3.130665812444000e+03 *sec, 'Cosmic/Cosmic/CP23/TIME',  4.069334187556000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP24',
      Active = True,
      Propagator = 'DIVA',
      Time = '10-JUL-2045 00:08:19.855159316 ET',
      Body = 'encnf5',
      Center = 'Enceladus',
      Frame = 'IAU Enceladus Fixed',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.semiMajorAxis(  3.810338516805570e+02 *km ),
         Conic.eccentricity(  2.256497606706511e-02 ),
         Conic.inclination(  9.575833371064664e+01 *deg ),
         Conic.longitudeOfNode(  7.222147989542890e+01 *deg ),
         Conic.argumentOfPeriapsis(  2.090800999548818e+02 *deg ),
         Conic.trueAnomaly(  1.281430626271321e+01 *deg ),
         ],
      Controls = [
         [ -4.099855159316000e+03 *sec, 'Cosmic/Cosmic/CP24/TIME',  3.100144840684000e+03 *sec ],
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
      PosTol =  1.000000000000000e-06 *km,
      VelTol =  1.000000000000000e-09 *km/sec,
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
         Cartesian.x(  1.302627735635875e+02 *km ),
         Cartesian.y( -6.768699089951889e+01 *km ),
         Cartesian.z( -3.603170526412953e+02 *km ),
         Cartesian.dx( -1.239640123160565e-01 *km/sec ),
         Cartesian.dy(  1.147274737082190e-02 *km/sec ),
         Cartesian.dz( -4.985107650756116e-02 *km/sec ),
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
