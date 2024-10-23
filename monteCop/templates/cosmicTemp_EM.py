# ===========================================================================
# Section 392 Navigation and Mission Design
#
# Copyright (C) 2021, California Institute of Technology.
# U.S. Government Sponsorship under NASA Contract NAS7-03001 is acknowledged.
#
# ===========================================================================

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
Timeline = []

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
