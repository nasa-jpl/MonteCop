#===========================================================================
#
# Section 392 Navigation and Mission Design
#
# Copyright (C) 2018, California Institute of Technology.
# U.S. Government Sponsorship under NASA Contract NAS7-03001 is acknowledged.
#
#===========================================================================

"""
    Description: COSMIC  Tempate for Earth-Moon System Problems.

    by: Ricardo L. Restrepo

    last edition: June 23 2022
"""

#--------------------------------------------------------------------------
# Setup Params:

scName = "mySC"
#scName = "lro"


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


# SNOPT Params:
MajorStepLimit = 0.1
LinesearchTol  = 0.99
MaxIterLimit   = 300                 #Major iterations limit

#--------------------------------------------------------------------------
# Configure the BOA database.
DefaultData(
   Inputs = [
      "time",
      "body",
      "frame/IAU 2000",
      "frame/inertial",
      "ephem/planet/de430",
      ]
   )

# Set up the default input values.
Defaults(
   Body = scName,
   Mass = 1000 *kg,
   Propagator = "DIVA",
   PosTol     = PosTolIn *km,
   VelTol     = VelTolIn *km/sec,
   MassTol = 0 *kg,
   )

# Set up the optimizer: SNOPT
NewSnOpt(
   Name = "SNOPT",
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
   } )

# Set up the optimizer: DBLSE
NewDblseOpt(
   Name             = 'DBLSE',
   MaxIterations    = 100,
   DblseTol         = 1.0E-9,
   CostRelTol       = 1.0E-2,
   CostAbsTol       = 1.0E-2,
   ConstraintRelTol = 1.0E-4,
)

# Create the Cosmic application.
NewCosmic(
   Optimizer = "SNOPT",
   #Optimizer = "DBLSE",     # USE iter(100,0.3) to converge (stepSize  0.1<x< 0.4 works)
   Cost = "MIN DV",
   MassModel = False,
   )

# Set up the propagator to use.
NewDivaPropagator(
   Name = "DIVA",
   StateTol = 1e-8,
   MassTol = 1e-6,
   MinStep = 1e-3 *sec,
   MaxStep = 100 *day,
   Forces = [
      "Gravity",
      'Impulse Burn',
      ],
   )

# Set up GUI:
#guiON = True
guiON = False
if guiON:
    Plugins = [
       'Bounds GUI',
       #'Cost Plot',
       NewTrajViewer(
          Name = 'COSMIC 3D Viewer',
          Bodies = [ 'Earth','Moon'],
          #Bodies = [ 'Earth'],
          ),
       ]

#Define Frames:
M.BodyVelDirFrame(boa, 'VUW_Earth', 'EMO2000', M.TimeInterval(), scName, 'Earth')
M.BodyVelDirFrame(boa, 'VUW_Moon', 'EMO2000', M.TimeInterval(), scName, 'Moon')

# NewBodyVelDirFrame(
#    Frame = 'VUW_Earth',
#    BaseFrame = CoordName.EME2000,
#    Body = scName,
#    Center = 'Earth',
#    )

# Set up the force model(s) to use for the propagation.
NewGravityBasic(
   PointMass = [ "Sun", "Earth", "Moon" ],
   )

# Add Spherical Harmonics here:

#--------------------------------------------------------------------------
# Set up the timeline:
#--------------------------------------------------------------------------
Timeline = [
   ControlPoint(
      Name = 'CP00',
      Active = True,
      Propagator = 'DIVA',
      Time = '27-NOV-2024 10:12:04.0465 ET',
      Body = 'mySC',
      Center = 'Moon',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Cartesian.x(  2.701978227288027e+01 *km ),
         Cartesian.y(  3.274051792782518e+02 *km ),
         Cartesian.z( -1.722349366402676e+03 *km ),
         Cartesian.dx( -4.422574927160464e-01 *km/sec ),
         Cartesian.dy(  1.604987497420431e+00 *km/sec ),
         Cartesian.dz(  2.981575795177492e-01 *km/sec ),
         ],
      Controls = [
         ],
      ),
   BreakPoint(
      Name = 'CP00->CP-DV01',
      Mode = 'MATCH_ELEM',
      Time = '27-NOV-2024 10:40:54.0465 ET',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'mySC',
      Name = 'DV01',
      Start = NewCosmicEvent(
         Name = 'CP-DV01',
         Delta = '-00:00:10',
         ),
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  1.883262244484874e-02 *km/sec,
      # DeltaVelRA  =  2.908721218389679e+02 *deg,
      # DeltaVelDec = -1.019214283981652e+01 *deg,
      DeltaVel = [
          6.603870200571526e-03 *km/sec,
         -1.731910761027349e-02 *km/sec,
         -3.332428259991626e-03 *km/sec,
         ],
      Controls = [
         [ -4.781665207040106e-02 *km/sec, 'DX',  4.781665207040106e-02 *km/sec ],
         [ -4.781665207040106e-02 *km/sec, 'DY',  4.781665207040106e-02 *km/sec ],
         [ -4.781665207040106e-02 *km/sec, 'DZ',  4.781665207040106e-02 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'CP-DV01',
      Active = True,
      Propagator = 'DIVA',
      Time = '27-NOV-2024 11:09:10.5926 ET',
      Body = 'mySC',
      Center = 'Moon',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  9.849769717733449e+01 *km ),
         Conic.apoapsisAltitude(  9.909894311859073e+01 *km ),
         Conic.inclination(  8.635117170926463e+01 *deg ),
         Conic.longitudeOfNode(  1.048162191996463e+02 *deg ),
         Conic.argumentOfPeriapsis( -4.088509080070889e+01 *deg ),
         Conic.trueAnomaly(  1.418263203412904e+02 *deg ),
         ],
      Controls = [
         [ -1.796654606341500e+04 *sec, 'Cosmic/Cosmic/CP-DV01/TIME',  1.803345393658500e+04 *sec ],
         [ None, 'Cosmic/Cosmic/CP-DV01/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP-DV01/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP-DV01/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP-DV01->CP-DV02',
      Mode = 'MATCH_ELEM',
      Time = '27-NOV-2024 14:03:39.0465 ET',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'mySC',
      Name = 'DV02',
      Start = NewCosmicEvent(
         Name = 'CP-DV02',
         Delta = '-00:00:10',
         ),
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  6.504721250764025e-01 *km/sec,
      # DeltaVelRA  =  2.888888955923936e+02 *deg,
      # DeltaVelDec =  5.375069704730667e+00 *deg,
      DeltaVel = [
          2.096540188475240e-01 *km/sec,
         -6.127367630453759e-01 *km/sec,
          6.093305420843258e-02 *km/sec,
         ],
      Controls = [
         [ -1.301583208878472e+00 *km/sec, 'DX',  1.301583208878472e+00 *km/sec ],
         [ -1.301583208878472e+00 *km/sec, 'DY',  1.301583208878472e+00 *km/sec ],
         [ -1.301583208878472e+00 *km/sec, 'DZ',  1.301583208878472e+00 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'CP-DV02',
      Active = True,
      Propagator = 'DIVA',
      Time = '27-NOV-2024 16:56:55.7285 ET',
      Body = 'mySC',
      Center = 'Moon',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  9.849906344147692e+01 *km ),
         Conic.apoapsisAltitude(  7.327695178676526e+04 *km ),
         Conic.inclination(  8.647915281279938e+01 *deg ),
         Conic.longitudeOfNode(  1.060650572366789e+02 *deg ),
         Conic.argumentOfPeriapsis(  8.459538061439049e+01 *deg ),
         Conic.trueAnomaly(  4.536042358002669e-01 *deg ),
         ],
      Controls = [
         [ -1.796168199555800e+04 *sec, 'Cosmic/Cosmic/CP-DV02/TIME',  1.803831800444200e+04 *sec ],
         [ None, 'Cosmic/Cosmic/CP-DV02/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP-DV02/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP-DV02/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP-DV02->CP-DV03',
      Mode = 'MATCH_ELEM',
      Time = '27-NOV-2024 22:58:29.0465 ET',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   OptImpulseBurn(
      Body = 'mySC',
      Name = 'DV03',
      Start = NewCosmicEvent(
         Name = 'CP-DV03',
         Delta = '-00:00:10',
         ),
      Frame = 'EMO2000',
      DeltaMass =  0.000000000000000e+00 *kg,
      # DeltaVelMag =  6.390323207318507e-02 *km/sec,
      # DeltaVelRA  =  1.552399251918780e+02 *deg,
      # DeltaVelDec = -2.793600660287099e+01 *deg,
      DeltaVel = [
         -5.126658929441335e-02 *km/sec,
          2.364515276765519e-02 *km/sec,
         -2.993771270674398e-02 *km/sec,
         ],
      Controls = [
         [ -1.287954551886715e-01 *km/sec, 'DX',  1.287954551886715e-01 *km/sec ],
         [ -1.287954551886715e-01 *km/sec, 'DY',  1.287954551886715e-01 *km/sec ],
         [ -1.287954551886715e-01 *km/sec, 'DZ',  1.287954551886715e-01 *km/sec ],
         ],
      ),
   ControlPoint(
      Name = 'CP-DV03',
      Active = True,
      Propagator = 'DIVA',
      Time = '28-NOV-2024 05:09:17.0186 ET',
      Body = 'mySC',
      Center = 'Moon',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Conic.periapsisAltitude(  1.711030493192100e+03 *km ),
         Conic.apoapsisAltitude(  8.411163709159468e+04 *km ),
         Conic.inclination(  9.304609604842031e+01 *deg ),
         Conic.longitudeOfNode(  1.182974967877064e+02 *deg ),
         Conic.argumentOfPeriapsis(  9.472486593659633e+01 *deg ),
         Conic.trueAnomaly(  1.473884525667292e+02 *deg ),
         ],
      Controls = [
         [ -1.859297211624500e+04 *sec, 'Cosmic/Cosmic/CP-DV03/TIME',  1.740702788375500e+04 *sec ],
         [ None, 'Cosmic/Cosmic/CP-DV03/Conic.longitudeOfNode', None ],
         [ None, 'Cosmic/Cosmic/CP-DV03/Conic.argumentOfPeriapsis', None ],
         [ None, 'Cosmic/Cosmic/CP-DV03/Conic.trueAnomaly', None ],
         ],
      ),
   BreakPoint(
      Name = 'CP-DV03->CP-END',
      Mode = 'MATCH_ELEM',
      Time = '28-NOV-2024 07:35:44.0465 ET',
      PosTol =  1.000000000000000e-01 *km,
      VelTol =  1.000000000000000e-04 *km/sec,
      MassTol =  0.000000000000000e+00 *kg,
      ),
   ControlPoint(
      Name = 'CP-END',
      Active = True,
      Propagator = 'DIVA',
      Time = '28-NOV-2024 10:12:04.0465 ET',
      Body = 'mySC',
      Center = 'Moon',
      Frame = 'EMO2000',
      Mass =  1.000000000000000e+03 *kg,
      State = [
         Cartesian.x(  5.088807786296963e+03 *km ),
         Cartesian.y( -1.330420333434784e+04 *km ),
         Cartesian.z( -3.402716706443689e+04 *km ),
         Cartesian.dx( -2.223871518028797e-02 *km/sec ),
         Cartesian.dy( -6.366729417369155e-03 *km/sec ),
         Cartesian.dz( -3.908462527759457e-01 *km/sec ),
         ],
      Controls = [
         ],
      ),
   ]



# =============================================================================
# runCosmic(): function added to be able to run cosmic with step=0.3 by default
#              !!** => USE ONLY WITH DBLSE AS YOUR OPTIMIZER **
# ============================================================================
# # These commands will run automatically in non-interactive mode
# def runCosmic():
#     Cosmic.run(step=0.3) #NOTE: step=0.1-0.4 works for this problem. Optimal was 0.3
#     # Cosmic.iter(500,step=0.3) #NOTE: step=0.1-0.4 works for this problem. Optimal was 0.3
#     # Cosmic.saveChkPt( "chkpt-1.py", allowOverwrite=True )
#     # # ... can add more functions here.
