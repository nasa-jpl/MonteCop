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
scName = 'Europa Clipper'
startTime = datetime.datetime.now()
stastOut='stat.out'

#--------------------------------------------------------------------------
# Setup problem tolerances:
#--------------------------------------------------------------------------
PosTolIn     = 1.0E-1
VelTolIn     = 1.0E-4
MajorFeasTol = 1e-5
MajorOptTol  = 1e-3

#--------------------------------------------------------------------------
# Setup SNOPT params for cosmicBatch() run:
#--------------------------------------------------------------------------
MaxIterLimit   = 300                 #Major iterations limit

#set SNOPT params (these values are overwritten on cosmicBatch.py calls):
MajorStepLimit = 2.0
LinesearchTol  = 0.001

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
      "$euclip_ephem/de430.boa",
      "/nav/common/import/ephem/mar097.boa",
      "$euclip_ephem/jup310.boa",
      # EOP File
      "$euclip_eop/simulated.eop",
      # Lockfile. Always load last
      #"$euclip_mdlock/euclip_MDLockfile_jup310_v1.boa",
      "/nav/euclip/common/lockfile/euclip_lockfile_jup310_v6.boa",
   ]
)

scName = 'Europa Clipper'

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
      },
   )

NewDivaPropagator(
    Name     = "DIVA", # Has to be the same name as the one used in defaults above
    StateTol = 1e-10,
    MassTol  = 1e-6,
    MinStep  = 1e-6 *sec,
    MaxStep  = 100 *day,
    Forces   = [
       "Gravity",
       "Impulse Burn",
       "Finite Burn"
       ]
    )
#==============================================================================
#
# Optimizer inputs
#
#==============================================================================

Plugins = []

#===========================================================================
#
# Constraints inputs
#
#===========================================================================
# Create dictionary of user defined tolerances for translatScaler
userTol = {}

#===========================================================================
#
# Optimizer inputs
#
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

# Build mean equator and equinox frames for the Jupiter bodies
#
NewMeanEquatorEquinoxFrame(
   Frame = "Jupiter Mean Equator And Equinox of Epoch",
   Body = "Jupiter",
   )

fb = [ "Callisto", "Ganymede", "Europa", "Io" ]
for body in fb:
   NewMeanEquatorEquinoxFrame(
      Frame = "%s Mean Equator And Equinox of Epoch" % body,
      PoleBody = body,
      OrbitBody = "Jupiter",
      )




#===========================================================================
#
# Cosmic Timeline
#
#===========================================================================

Timeline = [
   ControlPoint(
      Name = 'E04',
      Active = True,
      Propagator = 'DIVA',
      Time = '07-AUG-2032 01:09:54.886926996 ET',
      Body = 'Europa Clipper',
      Center = 'Europa',
      Frame = 'Europa Mean Equator And Equinox of Epoch',
      Mass =  0.000000000000000e+00 *kg,
      State = [
         Conic.periapsisAltitude(  2.752425993202381e+01 *km ),
         Conic.bPlaneTheta( -7.051017444799449e+01 *deg ),
         Conic.timeFromPeriapsis(  3.185508969791012e-13 *sec ),
         Conic.vInfinity(  4.492743026608124e+00 *km/sec ),
         Conic.inboundDec( -5.265169653570763e+00 *deg ),
         Conic.inboundRA(  8.685712010684547e+01 *deg ),
         ],
      Controls = [
         [ -1.728222604432380e+05 *sec, 'TIME',  1.727777395567620e+05 *sec ],
         [  2.500000000000000e+01 *km, 'Conic.periapsisAltitude',  1.000000000000124e+02 *km ],
         [ -3.600000000000000e+02 *deg, 'Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.042587551739654e+00 *km/sec, 'Conic.vInfinity',  4.940940341015134e+00 *km/sec ],
         [ None, 'Conic.inboundDec', None ],
         [ None, 'Conic.inboundRA', None ],
         ],
      ),
   OptImpulseBurn(
      Body = 'Europa Clipper',
      Name = 'E04-CU',
      Start = NewCosmicEvent(
         Name = 'E04',
         Delta = '3/00:00:00',
         ),
      Frame = 'EME2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.616054368416382e-08 *km/sec,
      # DeltaVelRA  =  3.329645866155159e+01 *deg,
      # DeltaVelDec = -2.927710699744100e+01 *deg,
      DeltaVel = [
          1.178224658440818e-08 *km/sec,
          7.738446989630684e-09 *km/sec,
         -7.903054838519385e-09 *km/sec,
         ],
      Controls = [
         CI( -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         ],
      ),
   BreakPoint(
      Name = 'E04->APO-E04',
      Mode = 'MATCH_ELEM',
      Time = '12-AUG-2032 09:00:20.516558323 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  1.000000000000000e-03 *kg,
      ),
   ControlPoint(
      Name = 'APO-E04',
      Active = True,
      Propagator = 'DIVA',
      Time = '17-AUG-2032 16:54:55.425043654 ET',
      Body = 'Europa Clipper',
      Center = 'Jupiter',
      Frame = 'EMO2000',
      Mass =  0.000000000000000e+00 *kg,
      State = [
         Conic.period(  1.840027476035891e+06 *sec ),
         Conic.periapsisRange(  6.643781966524181e+05 *km ),
         Conic.inclination(  5.552028157050645e+00 *deg ),
         Conic.longitudeOfNode( -2.124079646828037e+01 *deg ),
         Conic.argumentOfPeriapsis(  1.644042570581340e+02 *deg ),
         Conic.trueAnomaly(  1.804647801756878e+02 *deg ),
         ],
      Controls = [
         [ -1.730270184107670e+05 *sec, 'TIME',  1.725729815892330e+05 *sec ],
         [  1.498725157165174e+06 *sec, 'Conic.period',  2.248087735747761e+06 *sec ],
         [  3.323518240823094e+05 *km, 'Conic.periapsisRange',  9.970554722469284e+05 *km ],
         [  1.000000000000000e-02 *deg, 'Conic.inclination',  1.000000000000000e+01 *deg ],
         [ None, 'Conic.longitudeOfNode', None ],
         [ None, 'Conic.argumentOfPeriapsis', None ],
         [ None, 'Conic.trueAnomaly', None ],
         ],
      ),
   OptImpulseBurn(
      Body = 'Europa Clipper',
      Name = 'E05-TRG',
      Start = NewCosmicEvent(
         Name = 'APO-E04',
         Delta = '00:00:01',
         ),
      Frame = 'EME2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  3.886564582615740e-03 *km/sec,
      # DeltaVelRA  =  5.130867971903196e+01 *deg,
      # DeltaVelDec =  2.072613530508853e+01 *deg,
      DeltaVel = [
          2.272350156587584e-03 *km/sec,
          2.837237368336545e-03 *km/sec,
          1.375461063212682e-03 *km/sec,
         ],
      Controls = [
         CI( -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         ],
      ),
   BreakPoint(
      Name = 'APO-E04->E05',
      Mode = 'MATCH_ELEM',
      Time = '23-AUG-2032 00:41:56.296707453 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  1.000000000000000e-03 *kg,
      ),
   ControlPoint(
      Name = 'E05',
      Active = True,
      Propagator = 'DIVA',
      Time = '28-AUG-2032 08:32:18.575019646 ET',
      Body = 'Europa Clipper',
      Center = 'Europa',
      Frame = 'Europa Mean Equator And Equinox of Epoch',
      Mass =  0.000000000000000e+00 *kg,
      State = [
         Conic.periapsisAltitude(  2.650097941621978e+01 *km ),
         Conic.bPlaneTheta( -4.781717953427032e+01 *deg ),
         Conic.timeFromPeriapsis( -1.592754484895502e-13 *sec ),
         Conic.vInfinity(  4.502429463606267e+00 *km/sec ),
         Conic.inboundDec( -1.503311655055047e+01 *deg ),
         Conic.inboundRA(  9.027109235733543e+01 *deg ),
         ],
      Controls = [
         [ -1.727743882376280e+05 *sec, 'TIME',  1.728256117623720e+05 *sec ],
         [  2.500000000000000e+01 *km, 'Conic.periapsisAltitude',  1.000000000000004e+02 *km ],
         [ -3.600000000000000e+02 *deg, 'Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.042587551739654e+00 *km/sec, 'Conic.vInfinity',  4.940940341015132e+00 *km/sec ],
         [ None, 'Conic.inboundDec', None ],
         [ None, 'Conic.inboundRA', None ],
         ],
      ),
   OptImpulseBurn(
      Body = 'Europa Clipper',
      Name = 'E05-CU',
      Start = NewCosmicEvent(
         Name = 'E05',
         Delta = '3/00:00:00',
         ),
      Frame = 'EME2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  7.798916616466324e-07 *km/sec,
      # DeltaVelRA  =  1.658096490763530e+01 *deg,
      # DeltaVelDec =  1.219879471786661e+01 *deg,
      DeltaVel = [
          7.305843325889575e-07 *km/sec,
          2.175323858929155e-07 *km/sec,
          1.647944112196427e-07 *km/sec,
         ],
      Controls = [
         CI( -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         ],
      ),
   BreakPoint(
      Name = 'E05->APO-E05',
      Mode = 'MATCH_ELEM',
      Time = '02-SEP-2032 16:23:49.205616340 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  1.000000000000000e-03 *kg,
      ),
   ControlPoint(
      Name = 'APO-E05',
      Active = True,
      Propagator = 'DIVA',
      Time = '08-SEP-2032 00:33:58.005022792 ET',
      Body = 'Europa Clipper',
      Center = 'Jupiter',
      Frame = 'EMO2000',
      Mass =  0.000000000000000e+00 *kg,
      State = [
         Conic.period(  1.839974686498191e+06 *sec ),
         Conic.periapsisRange(  6.704639820118916e+05 *km ),
         Conic.inclination(  7.355070133142943e+00 *deg ),
         Conic.longitudeOfNode( -2.198183504428530e+01 *deg ),
         Conic.argumentOfPeriapsis(  1.696440972421559e+02 *deg ),
         Conic.trueAnomaly(  1.803739674995933e+02 *deg ),
         ],
      Controls = [
         [ -1.739437805721300e+05 *sec, 'TIME',  1.716562194278700e+05 *sec ],
         [  1.506893877980984e+06 *sec, 'Conic.period',  2.260340816971476e+06 *sec ],
         [  3.354057820633669e+05 *km, 'Conic.periapsisRange',  1.006217346190101e+06 *km ],
         [  1.000000000000000e-02 *deg, 'Conic.inclination',  1.000000000000000e+01 *deg ],
         [ None, 'Conic.longitudeOfNode', None ],
         [ None, 'Conic.argumentOfPeriapsis', None ],
         [ None, 'Conic.trueAnomaly', None ],
         ],
      ),
   OptImpulseBurn(
      Body = 'Europa Clipper',
      Name = 'E06-TRG',
      Start = NewCosmicEvent(
         Name = 'APO-E05',
         Delta = '00:00:01',
         ),
      Frame = 'EME2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  4.085355148311090e-03 *km/sec,
      # DeltaVelRA  =  5.514856512754878e+01 *deg,
      # DeltaVelDec =  2.354642601727973e+01 *deg,
      DeltaVel = [
          2.140193450357918e-03 *km/sec,
          3.073443786843779e-03 *km/sec,
          1.632066779267430e-03 *km/sec,
         ],
      Controls = [
         CI( -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         ],
      ),
   BreakPoint(
      Name = 'APO-E05->E06',
      Mode = 'MATCH_ELEM',
      Time = '13-SEP-2032 08:05:59.243284985 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  1.000000000000000e-03 *kg,
      ),
   ControlPoint(
      Name = 'E06',
      Active = True,
      Propagator = 'DIVA',
      Time = '18-SEP-2032 15:43:14.027125138 ET',
      Body = 'Europa Clipper',
      Center = 'Europa',
      Frame = 'Europa Mean Equator And Equinox of Epoch',
      Mass =  0.000000000000000e+00 *kg,
      State = [
         Conic.periapsisAltitude(  2.500000000000000e+01 *km ),
         Conic.bPlaneTheta( -2.158470004234120e+01 *deg ),
         Conic.timeFromPeriapsis( -3.981886212238773e-14 *sec ),
         Conic.vInfinity(  4.488740906428987e+00 *km/sec ),
         Conic.inboundDec( -2.226401546782742e+01 *deg ),
         Conic.inboundRA(  9.789112232951099e+01 *deg ),
         ],
      Controls = [
         [ -1.719697650058300e+05 *sec, 'TIME',  1.736302349941700e+05 *sec ],
         [  2.500000000000000e+01 *km, 'Conic.periapsisAltitude',  1.000000000000178e+02 *km ],
         [ -3.600000000000000e+02 *deg, 'Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.042587551739655e+00 *km/sec, 'Conic.vInfinity',  4.940940341015135e+00 *km/sec ],
         [ None, 'Conic.inboundDec', None ],
         [ None, 'Conic.inboundRA', None ],
         ],
      ),
   OptImpulseBurn(
      Body = 'Europa Clipper',
      Name = 'E06-CU',
      Start = NewCosmicEvent(
         Name = 'E06',
         Delta = '3/00:00:00',
         ),
      Frame = 'EME2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  3.416369483011598e-08 *km/sec,
      # DeltaVelRA  =  3.585500543564655e+01 *deg,
      # DeltaVelDec =  1.764872372908820e+01 *deg,
      DeltaVel = [
          2.638647093510562e-08 *km/sec,
          1.906905801675198e-08 *km/sec,
          1.035776145650533e-08 *km/sec,
         ],
      Controls = [
         CI( -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         ],
      ),
   BreakPoint(
      Name = 'E06->APO-E06',
      Mode = 'MATCH_ELEM',
      Time = '23-SEP-2032 23:48:23.986388926 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  1.000000000000000e-03 *kg,
      ),
   ControlPoint(
      Name = 'APO-E06',
      Active = True,
      Propagator = 'DIVA',
      Time = '29-SEP-2032 07:14:12.515125629 ET',
      Body = 'Europa Clipper',
      Center = 'Jupiter',
      Frame = 'EMO2000',
      Mass =  0.000000000000000e+00 *kg,
      State = [
         Conic.period(  1.840700065603949e+06 *sec ),
         Conic.periapsisRange(  6.760641565974976e+05 *km ),
         Conic.inclination(  8.055485085239731e+00 *deg ),
         Conic.longitudeOfNode( -2.318926210583562e+01 *deg ),
         Conic.argumentOfPeriapsis(  1.767484327387147e+02 *deg ),
         Conic.trueAnomaly( -1.799394499020320e+02 *deg ),
         ],
      Controls = [
         [ -1.712688044670860e+05 *sec, 'TIME',  1.743311955329140e+05 *sec ],
         [  1.503627685507271e+06 *sec, 'Conic.period',  2.255441528260907e+06 *sec ],
         [  3.389763232426511e+05 *km, 'Conic.periapsisRange',  1.016928969727953e+06 *km ],
         [  1.000000000000000e-02 *deg, 'Conic.inclination',  1.000000000000000e+01 *deg ],
         [ None, 'Conic.longitudeOfNode', None ],
         [ None, 'Conic.argumentOfPeriapsis', None ],
         [ None, 'Conic.trueAnomaly', None ],
         ],
      ),
   OptImpulseBurn(
      Body = 'Europa Clipper',
      Name = 'E07-TRG',
      Start = NewCosmicEvent(
         Name = 'APO-E06',
         Delta = '00:00:01',
         ),
      Frame = 'EME2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  3.522103556376202e-03 *km/sec,
      # DeltaVelRA  =  6.028382932637025e+01 *deg,
      # DeltaVelDec =  2.584179536064048e+01 *deg,
      DeltaVel = [
          1.571329942504119e-03 *km/sec,
          2.753028960466335e-03 *km/sec,
          1.535241745283636e-03 *km/sec,
         ],
      Controls = [
         CI( -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         ],
      ),
   BreakPoint(
      Name = 'APO-E06->E07',
      Mode = 'MATCH_ELEM',
      Time = '04-OCT-2032 15:31:03.434928161 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  1.000000000000000e-03 *kg,
      ),
   ControlPoint(
      Name = 'E07',
      Active = True,
      Propagator = 'DIVA',
      Time = '09-OCT-2032 23:09:40.942565713 ET',
      Body = 'Europa Clipper',
      Center = 'Europa',
      Frame = 'Europa Mean Equator And Equinox of Epoch',
      Mass =  0.000000000000000e+00 *kg,
      State = [
         Conic.periapsisAltitude(  2.500000000000000e+01 *km ),
         Conic.bPlaneTheta(  4.552301146025020e+00 *deg ),
         Conic.timeFromPeriapsis(  1.592754484895500e-13 *sec ),
         Conic.vInfinity(  4.494934007350437e+00 *km/sec ),
         Conic.inboundDec( -2.568151744887630e+01 *deg ),
         Conic.inboundRA(  1.080238673603783e+02 *deg ),
         ],
      Controls = [
         [ -1.720377833679350e+05 *sec, 'TIME',  1.735622166320650e+05 *sec ],
         [  2.500000000000000e+01 *km, 'Conic.periapsisAltitude',  1.000000000000000e+02 *km ],
         [ -3.600000000000000e+02 *deg, 'Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.042587551739656e+00 *km/sec, 'Conic.vInfinity',  4.940940341015136e+00 *km/sec ],
         [ None, 'Conic.inboundDec', None ],
         [ None, 'Conic.inboundRA', None ],
         ],
      ),
   OptImpulseBurn(
      Body = 'Europa Clipper',
      Name = 'E07-CU',
      Start = NewCosmicEvent(
         Name = 'E07',
         Delta = '3/00:00:00',
         ),
      Frame = 'EME2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  2.258593146335174e-07 *km/sec,
      # DeltaVelRA  =  1.914692156435872e+01 *deg,
      # DeltaVelDec =  1.968598939967166e+01 *deg,
      DeltaVel = [
          2.008943677694251e-07 *km/sec,
          6.975021160383033e-08 *km/sec,
          7.608410468991603e-08 *km/sec,
         ],
      Controls = [
         CI( -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         ],
      ),
   BreakPoint(
      Name = 'E07->APO-E07',
      Mode = 'MATCH_ELEM',
      Time = '15-OCT-2032 07:13:54.142928988 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  1.000000000000000e-03 *kg,
      ),
   ControlPoint(
      Name = 'APO-E07',
      Active = True,
      Propagator = 'DIVA',
      Time = '20-OCT-2032 14:26:06.538853347 ET',
      Body = 'Europa Clipper',
      Center = 'Jupiter',
      Frame = 'EMO2000',
      Mass =  0.000000000000000e+00 *kg,
      State = [
         Conic.period(  1.840619266826253e+06 *sec ),
         Conic.periapsisRange(  6.764556721621931e+05 *km ),
         Conic.inclination(  7.766593147535714e+00 *deg ),
         Conic.longitudeOfNode( -2.349970135027947e+01 *deg ),
         Conic.argumentOfPeriapsis(  1.833026432044318e+02 *deg ),
         Conic.trueAnomaly( -1.801785403689352e+02 *deg ),
         ],
      Controls = [
         [ -1.704414121931490e+05 *sec, 'TIME',  1.751585878068510e+05 *sec ],
         [  1.505908290987992e+06 *sec, 'Conic.period',  2.258862436481989e+06 *sec ],
         [  3.391985927046147e+05 *km, 'Conic.periapsisRange',  1.017595778113844e+06 *km ],
         [  1.000000000000000e-02 *deg, 'Conic.inclination',  1.000000000000000e+01 *deg ],
         [ None, 'Conic.longitudeOfNode', None ],
         [ None, 'Conic.argumentOfPeriapsis', None ],
         [ None, 'Conic.trueAnomaly', None ],
         ],
      ),
   OptImpulseBurn(
      Body = 'Europa Clipper',
      Name = 'E08-TRG',
      Start = NewCosmicEvent(
         Name = 'APO-E07',
         Delta = '00:00:01',
         ),
      Frame = 'EME2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  4.025773707329879e-03 *km/sec,
      # DeltaVelRA  =  6.680434012071160e+01 *deg,
      # DeltaVelDec =  2.576308694943016e+01 *deg,
      DeltaVel = [
          1.428026391789061e-03 *km/sec,
          3.332532354937485e-03 *km/sec,
          1.749806466520595e-03 *km/sec,
         ],
      Controls = [
         CI( -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         ],
      ),
   BreakPoint(
      Name = 'APO-E07->E08',
      Mode = 'MATCH_ELEM',
      Time = '25-OCT-2032 22:56:56.110391408 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  1.000000000000000e-03 *kg,
      ),
   ControlPoint(
      Name = 'E08',
      Active = True,
      Propagator = 'DIVA',
      Time = '31-OCT-2032 06:36:24.077309154 ET',
      Body = 'Europa Clipper',
      Center = 'Europa',
      Frame = 'Europa Mean Equator And Equinox of Epoch',
      Mass =  0.000000000000000e+00 *kg,
      State = [
         Conic.periapsisAltitude(  2.500000000000000e+01 *km ),
         Conic.bPlaneTheta(  3.042612909814056e+01 *deg ),
         Conic.timeFromPeriapsis(  7.963772424477523e-14 *sec ),
         Conic.vInfinity(  4.490450084283217e+00 *km/sec ),
         Conic.inboundDec( -2.439915450332769e+01 *deg ),
         Conic.inboundRA(  1.186749207291335e+02 *deg ),
         ],
      Controls = [
         [ -1.720769831865360e+05 *sec, 'TIME',  1.735230168134640e+05 *sec ],
         [  2.500000000000000e+01 *km, 'Conic.periapsisAltitude',  1.000000000000068e+02 *km ],
         [ -3.600000000000000e+02 *deg, 'Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.042587551739655e+00 *km/sec, 'Conic.vInfinity',  4.940940341015135e+00 *km/sec ],
         [ None, 'Conic.inboundDec', None ],
         [ None, 'Conic.inboundRA', None ],
         ],
      ),
   OptImpulseBurn(
      Body = 'Europa Clipper',
      Name = 'E08-CU',
      Start = NewCosmicEvent(
         Name = 'E08',
         Delta = '3/00:00:00',
         ),
      Frame = 'EME2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  2.948247039460569e-08 *km/sec,
      # DeltaVelRA  =  3.518079095780865e+02 *deg,
      # DeltaVelDec =  2.200435176222134e+01 *deg,
      DeltaVel = [
          2.705590456781950e-08 *km/sec,
         -3.895002905520696e-09 *km/sec,
          1.104640398282762e-08 *km/sec,
         ],
      Controls = [
         CI( -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         ],
      ),
   BreakPoint(
      Name = 'E08->APO-E08',
      Mode = 'MATCH_ELEM',
      Time = '05-NOV-2032 14:40:05.772220353 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  1.000000000000000e-03 *kg,
      ),
   ControlPoint(
      Name = 'APO-E08',
      Active = True,
      Propagator = 'DIVA',
      Time = '10-NOV-2032 22:20:44.303609224 ET',
      Body = 'Europa Clipper',
      Center = 'Jupiter',
      Frame = 'EMO2000',
      Mass =  0.000000000000000e+00 *kg,
      State = [
         Conic.period(  1.840643226081334e+06 *sec ),
         Conic.periapsisRange(  6.740302722362081e+05 *km ),
         Conic.inclination(  6.459451015647628e+00 *deg ),
         Conic.longitudeOfNode( -2.378496793261530e+01 *deg ),
         Conic.argumentOfPeriapsis(  1.890776992402474e+02 *deg ),
         Conic.trueAnomaly( -1.802835233902248e+02 *deg ),
         ],
      Controls = [
         [ -1.721398532911360e+05 *sec, 'TIME',  1.734601467088640e+05 *sec ],
         [  1.501216384686046e+06 *sec, 'Conic.period',  2.251824577029069e+06 *sec ],
         [  3.372526301365082e+05 *km, 'Conic.periapsisRange',  1.011757890409525e+06 *km ],
         [  1.000000000000000e-02 *deg, 'Conic.inclination',  1.000000000000000e+01 *deg ],
         [ None, 'Conic.longitudeOfNode', None ],
         [ None, 'Conic.argumentOfPeriapsis', None ],
         [ None, 'Conic.trueAnomaly', None ],
         ],
      ),
   OptImpulseBurn(
      Body = 'Europa Clipper',
      Name = 'E09-TRG',
      Start = NewCosmicEvent(
         Name = 'APO-E08',
         Delta = '00:00:01',
         ),
      Frame = 'EME2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.854257083810411e-03 *km/sec,
      # DeltaVelRA  =  7.290441675817766e+01 *deg,
      # DeltaVelDec =  2.288483684892728e+01 *deg,
      DeltaVel = [
          5.021848238787347e-04 *km/sec,
          1.632825138923527e-03 *km/sec,
          7.210837685223124e-04 *km/sec,
         ],
      Controls = [
         CI( -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         ],
      ),
   BreakPoint(
      Name = 'APO-E08->E09',
      Mode = 'MATCH_ELEM',
      Time = '16-NOV-2032 06:23:23.128415823 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  1.000000000000000e-03 *kg,
      ),
   ControlPoint(
      Name = 'E09',
      Active = True,
      Propagator = 'DIVA',
      Time = '21-NOV-2032 14:03:30.821350115 ET',
      Body = 'Europa Clipper',
      Center = 'Europa',
      Frame = 'Europa Mean Equator And Equinox of Epoch',
      Mass =  0.000000000000000e+00 *kg,
      State = [
         Conic.periapsisAltitude(  2.500000000000000e+01 *km ),
         Conic.bPlaneTheta(  5.556782862234355e+01 *deg ),
         Conic.timeFromPeriapsis( -1.990943106119368e-13 *sec ),
         Conic.vInfinity(  4.478093816822052e+00 *km/sec ),
         Conic.inboundDec( -1.885069933637145e+01 *deg ),
         Conic.inboundRA(  1.274084180247702e+02 *deg ),
         ],
      Controls = [
         [ -1.721090148365570e+05 *sec, 'TIME',  1.734909851634430e+05 *sec ],
         [  2.500000000000000e+01 *km, 'Conic.periapsisAltitude',  1.000000000000000e+02 *km ],
         [ -3.600000000000000e+02 *deg, 'Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.042587551739654e+00 *km/sec, 'Conic.vInfinity',  4.940940341015134e+00 *km/sec ],
         [ None, 'Conic.inboundDec', None ],
         [ None, 'Conic.inboundRA', None ],
         ],
      ),
   OptImpulseBurn(
      Body = 'Europa Clipper',
      Name = 'E09-CU',
      Start = NewCosmicEvent(
         Name = 'E09',
         Delta = '3/00:00:00',
         ),
      Frame = 'EME2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  3.039962275103695e-08 *km/sec,
      # DeltaVelRA  =  1.843395517646975e+02 *deg,
      # DeltaVelDec = -4.993510586553806e+01 *deg,
      DeltaVel = [
         -1.951076855465541e-08 *km/sec,
         -1.480567356927601e-09 *km/sec,
         -2.326531525927886e-08 *km/sec,
         ],
      Controls = [
         CI( -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         ],
      ),
   BreakPoint(
      Name = 'E09->APO-E09',
      Mode = 'MATCH_ELEM',
      Time = '26-NOV-2032 22:06:41.661392856 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  1.000000000000000e-03 *kg,
      ),
   ControlPoint(
      Name = 'APO-E09',
      Active = True,
      Propagator = 'DIVA',
      Time = '02-DEC-2032 06:03:54.874682028 ET',
      Body = 'Europa Clipper',
      Center = 'Jupiter',
      Frame = 'EMO2000',
      Mass =  0.000000000000000e+00 *kg,
      State = [
         Conic.period(  1.841194616142313e+06 *sec ),
         Conic.periapsisRange(  6.707554501768809e+05 *km ),
         Conic.inclination(  4.320257446557174e+00 *deg ),
         Conic.longitudeOfNode( -2.411530143957883e+01 *deg ),
         Conic.argumentOfPeriapsis(  1.931121800779058e+02 *deg ),
         Conic.trueAnomaly( -1.803639153465596e+02 *deg ),
         ],
      Controls = [
         [ -1.731333584098750e+05 *sec, 'TIME',  1.724666415901250e+05 *sec ],
         [  1.491379503435723e+06 *sec, 'Conic.period',  2.237069255153585e+06 *sec ],
         [  3.348321406646205e+05 *km, 'Conic.periapsisRange',  1.004496421993861e+06 *km ],
         [  1.000000000000000e-02 *deg, 'Conic.inclination',  1.000000000000000e+01 *deg ],
         [ None, 'Conic.longitudeOfNode', None ],
         [ None, 'Conic.argumentOfPeriapsis', None ],
         [ None, 'Conic.trueAnomaly', None ],
         ],
      ),
   OptImpulseBurn(
      Body = 'Europa Clipper',
      Name = 'E10-TRG',
      Start = NewCosmicEvent(
         Name = 'APO-E09',
         Delta = '00:00:01',
         ),
      Frame = 'EME2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.090567122589725e-06 *km/sec,
      # DeltaVelRA  =  2.623499209907778e+02 *deg,
      # DeltaVelDec = -1.932702498066663e+00 *deg,
      DeltaVel = [
         -1.450966573670136e-07 *km/sec,
         -1.080245729548502e-06 *km/sec,
         -3.678005827586400e-08 *km/sec,
         ],
      Controls = [
         CI( -2.000000000000000e-02 *km/sec, 'DX',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DY',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         CI( -2.000000000000000e-02 *km/sec, 'DZ',  2.000000000000000e-02 *km/sec, Scale=0.001 ),
         ],
      ),
   BreakPoint(
      Name = 'APO-E09->E10',
      Mode = 'MATCH_ELEM',
      Time = '07-DEC-2032 13:50:01.371151451 ET',
      Frame = 'EMO2000',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  1.000000000000000e-03 *kg,
      ),
   ControlPoint(
      Name = 'E10',
      Active = True,
      Propagator = 'DIVA',
      Time = '12-DEC-2032 21:31:58.245565871 ET',
      Body = 'Europa Clipper',
      Center = 'Europa',
      Frame = 'Europa Mean Equator And Equinox of Epoch',
      Mass =  0.000000000000000e+00 *kg,
      State = [
         Conic.periapsisAltitude(  8.946542966261600e+01 *km ),
         Conic.bPlaneTheta(  8.020945057898476e+01 *deg ),
         Conic.timeFromPeriapsis( -1.329386852414805e-12 *sec ),
         Conic.vInfinity(  4.457009032947158e+00 *km/sec ),
         Conic.inboundDec( -1.014543070582893e+01 *deg ),
         Conic.inboundRA(  1.325340986018890e+02 *deg ),
         ],
      Controls = [
         [ -1.722170195351230e+05 *sec, 'TIME',  1.733829804648770e+05 *sec ],
         [  2.500000000000000e+01 *km, 'Conic.periapsisAltitude',  5.872838801347247e+02 *km ],
         [ -3.600000000000000e+02 *deg, 'Conic.bPlaneTheta',  3.600000000000000e+02 *deg ],
         [  4.042587551739655e+00 *km/sec, 'Conic.vInfinity',  4.940940341015135e+00 *km/sec ],
         [ None, 'Conic.inboundDec', None ],
         [ None, 'Conic.inboundRA', None ],
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
