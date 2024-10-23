# ===========================================================================
# Section 392 Navigation and Mission Design
#
# Copyright (C) 2018, California Institute of Technology.
# U.S. Government Sponsorship under NASA Contract NAS7-03001 is acknowledged.
#
# ===========================================================================

""" Cosmic utilities """

from __future__ import print_function

__version__ = "0.1"
__author__ = "Ricardo L. Restrepo (392M)"

# ===========================================================================
# imports here:

# import mpylab
import Monte as M
from mpy.units import s, sec, hour, m, km, day, deg, rad, kg

from cosmicUtils import *

# import bdaLib.util.cosmicUtils as bdalib
# import atnLib.utilities.trajInspector as insp

import mmath
# ===========================================================================


#==============================================================================
# Defaults:
#==============================================================================

scName  = 'mySC'
sc_ID = -99999   # TODO: Move to default (-99999)

centerDefault = 'Earth'    # TODO: Move to default
frameDefault = 'EME2000'

# ADD:
# NewDivaPropagator
# - ...
# - ..
# ===========================================================================
# Functions:
# ===========================================================================

# ----------------------------------------------------------------------------
# quiet the manager output
def silenceMgr(mgr):
    """quiet the manager output

    def silenceMgr(mgr):
    """

    # Quiet managers
    mgr.outputMgr.quiet = True
    mgr._quiet = True

# ----------------------------------------------------------------------------
# undo silenceMgr
def unsilenceMgr(mgr):
    """undo silenceMgr

    def unsilenceMgr(mgr):
    """

    # undo quiet managers
    mgr.outputMgr.quiet = False
    mgr._quiet = False

# ----------------------------------------------------------------------------

# check if manager is silenced
def mgrIsSilent(mgr):
    """check if manager is silenced

    def mgrIsSilent(mgr):
    """

    # get manager quiet condition
    cond1 = mgr.outputMgr.quiet
    cond2 = mgr._quiet
    isQuiet = cond1 & cond2

    return isQuiet

# ===========================================================================
# SCAN BSP TRAJECTORY::

# ---------------------------------------------------------------------------
#  Find Peri/Apos

# ---------------------------------------------------------------------------
# Find DV's discontinuities

# ---------------------------------------------------------------------------
# Find DV's discontinuities


# ===========================================================================
# ADD CONTROL POINTS (BPs) and MANEUVERS (DVs)

# ---------------------------------------------------------------------------
# Add DVs




# ===========================================================================
# append control point given in Classical Orbital Elements
# OE: sma + ecc
def appendCpCoe(
    mgr,
    cpName,
    cpTime,
    sma,
    ecc,
    inc,
    aop,
    raan,
    tra,
    mass=1000*kg,
    center=centerDefault,
    frame=frameDefault,
    body=scName,
    propagator="DIVA",
):
    """ append control point given in Classical Orbital Elements

    def appendCpCoe(
        mgr,
        cpName,
        cpTime,
        sma,
        ecc,
        inc,
        argp,
        raan,
        tra,
        mass=1000*kg,
        center=centerDefault,
        frame=frameDefault,
        body=scName,
        propagator="DIVA",
    ):

    ex: appendCpCoe(mgr,'RP1',t0, ecc,rpAlt,inc,aop,raan,tra)
    """

    # silence manager temporarily for quiet execution
    silenceMgr(mgr)

    newState = [
        M.Conic.semiMajorAxis(sma),
        M.Conic.eccentricity(ecc),
        M.Conic.inclination(inc),
        M.Conic.longitudeOfNode(raan),
        M.Conic.argumentOfPeriapsis(aop),
        M.Conic.trueAnomaly(tra),
    ]

    # create Control Point and add it to the problem
    cp = M.ControlPoint(
        mgr.boa,
        cpName,
        cpTime,
        body,
        center,
        frame,
        newState,
        mass,
        propagator,
    )
    mgr.cosmic.timeline().append(cp)

    # Create new control list
    newControls = M.OptControlList(mgr.boa)

    lbound = -mmath.inf
    ubound = mmath.inf

    # add controls in
    baseStr = "Cosmic/Cosmic/{0}/".format(cpName) #may not need this
    # baseStr = "" #may not need this
    newControls.add(baseStr + "TIME", -3600 * s, 3600 * s)
    newControls.add(baseStr +  "Conic.semiMajorAxis", 0.9 * sma, 1.1 * sma)
    newControls.add(baseStr + "Conic.eccentricity", 1e-5, ubound)
    newControls.add(baseStr + "Conic.inclination", 0.95 * inc, 1.05 * inc)
    #newControls.add(baseStr + "Conic.inclination", 88*deg, 108*deg)
    newControls.add(baseStr + "Conic.longitudeOfNode", lbound, ubound)
    newControls.add(baseStr + "Conic.argumentOfPeriapsis", lbound, ubound)
    newControls.add(baseStr + "Conic.trueAnomaly", lbound, ubound)

    # Add the new controls list to the new control point
    mgr.cp[cpName].controls().append(newControls)

    # un silence manager
    unsilenceMgr(mgr)

    # print new CP
    #print(mgr.cp[cpName])

# ===========================================================================

# ===========================================================================
# append control point given in Classical Orbital Elements
# OE_V2: ecc + rpAlt
def appendCpCoe_V2(
    mgr,
    cpName,
    cpTime,
    ecc,
    rpAlt,
    inc,
    aop,
    raan,
    tra,
    mass=1000*kg,
    center=centerDefault,
    frame=frameDefault,
    body=scName,
    propagator="DIVA",
):
    """ append control point given in Classical Orbital Elements

    def appendCpCoe(
        mgr,
        cpName,
        cpTime,
        ecc,
        rpAlt,
        inc,
        argp,
        raan,
        tra,
        mass=1000*kg,
        center=centerDefault,
        frame=frameDefault,
        body=scName,
        propagator="DIVA",
    ):

    ex: appendCpCoe(mgr,'RP1',t0, ecc,rpAlt,inc,aop,raan,tra)
    """

    # silence manager temporarily for quiet execution
    silenceMgr(mgr)

    newState = [
        #M.Conic.semiMajorAxis(sma),
        M.Conic.eccentricity(ecc),
        M.Conic.periapsisAltitude(rpAlt),
        M.Conic.inclination(inc),
        M.Conic.longitudeOfNode(raan),
        M.Conic.argumentOfPeriapsis(aop),
        M.Conic.trueAnomaly(tra),
    ]

    # create Control Point and add it to the problem
    cp = M.ControlPoint(
        mgr.boa,
        cpName,
        cpTime,
        body,
        center,
        frame,
        newState,
        mass,
        propagator,
    )
    mgr.cosmic.timeline().append(cp)

    # Create new control list
    newControls = M.OptControlList(mgr.boa)

    lbound = -mmath.inf
    ubound = mmath.inf

    # add controls in
    baseStr = "Cosmic/Cosmic/{0}/".format(cpName) #may not need this
    #newControls.add( "TIME", -3600 * s, 3600 * s)
    newControls.add(baseStr + "TIME", -86400 * s, 86400 * s)
    newControls.add(baseStr + "Conic.eccentricity", 1e-5, 0.07)
    #newControls.add(baseStr + "Conic.periapsisAltitude", 10*km, 35*km)
    newControls.add(baseStr + "Conic.periapsisAltitude", 10*km, ubound)
    newControls.add(baseStr + "Conic.inclination", 85*deg, 105*deg),
    newControls.add(baseStr + "Conic.longitudeOfNode", lbound, ubound)
    newControls.add(baseStr + "Conic.argumentOfPeriapsis", lbound, ubound)
    newControls.add(baseStr + "Conic.trueAnomaly", lbound, ubound)


    # Add the new controls list to the new control point
    mgr.cp[cpName].controls().append(newControls)

    # un silence manager
    unsilenceMgr(mgr)

    # print new CP
    #print(mgr.cp[cpName])

# ===========================================================================
# append control point given in Classical Orbital Elements
# OE3: rpAlt, raAlt, inc, ...
def appendCpCoe_V3(
    mgr,
    cpName,
    cpTime,
    stQuery,
    mass=1000*kg,
    center=centerDefault,
    frame=frameDefault,
    body=scName,
    propagator="DIVA",
    #controls =[]
):
    """ append control point given in Classical Orbital Elements

    def appendCpCoe(
        mgr,
        cpName,
        cpTime,
        stQuery,
        mass=1000*kg,
        center=centerDefault,
        frame=frameDefault,
        body=scName,
        propagator="DIVA",
        #controls =[]
    )

    stateParams{
        rpAlt,
        raAlt,
        inc,
        argp,
        raan,
        tra,
    }

    ex: appendCpCoe(mgr,- ,-,-,...)
    """

    # silence manager temporarily for quiet execution
    silenceMgr(mgr)

    #Read state from Query
    cpState= stQuery.state(cpTime).rotate(frame)

    newState = [
        # Seems redundant, but needed in order to make newState a 'Coordinate' list
        M.Conic.periapsisAltitude(M.Conic.periapsisAltitude(cpState)),
        M.Conic.apoapsisAltitude(M.Conic.apoapsisAltitude(cpState)),
        M.Conic.inclination(M.Conic.inclination(cpState)),
        M.Conic.longitudeOfNode(M.Conic.longitudeOfNode(cpState)),
        M.Conic.argumentOfPeriapsis(M.Conic.argumentOfPeriapsis(cpState)),
        M.Conic.trueAnomaly(M.Conic.trueAnomaly(cpState)),
    ]

    # create Control Point and add it to the problem
    cp = M.ControlPoint(
        mgr.boa,
        cpName,
        cpTime,
        body,
        center,
        frame,
        newState,
        mass,
        propagator,
    )
    mgr.cosmic.timeline().append(cp)

    # Create new control list
    newControls = M.OptControlList(mgr.boa)
    lbound = -mmath.inf
    ubound = mmath.inf

    # add controls in
    baseStr = "Cosmic/Cosmic/{0}/".format(cpName) #may not need this
    # baseStr = "" #may not need this
    newControls.add(baseStr + "TIME", -5*3600 * s, 5*3600 * s)
    #newControls.add(baseStr +  "Conic.periapsisAltitude", 0.9 * sma, 1.1 * sma)
    #newControls.add(baseStr + "Conic.apoapsisAltitude", 1e-5, ubound)
    #newControls.add(baseStr + "Conic.inclination", 0.95 * inc, 1.05 * inc)
    #newControls.add(baseStr + "Conic.inclination", 88*deg, 108*deg)
    newControls.add(baseStr + "Conic.longitudeOfNode", lbound, ubound)
    newControls.add(baseStr + "Conic.argumentOfPeriapsis", lbound, ubound)
    newControls.add(baseStr + "Conic.trueAnomaly", lbound, ubound)

    # Add the new controls list to the new control point
    mgr.cp[cpName].controls().append(newControls)

    # un silence manager
    unsilenceMgr(mgr)

# ===========================================================================

def appendCpCart(
    mgr,
    cpName,
    cpTime,
    xx,
    mass=1000*kg,
    center=centerDefault,
    frame=frameDefault,
    body=scName,
    propagator="DIVA",
    fixCP = False
):
    """ append control point given in Classical Orbital Elements

    def appendCpCoe(
        mgr,
        cpName,
        cpTime,
        xx,
        mass=1000*kg,
        center=centerDefault,
        frame=frameDefault,
        body=scName,
        propagator="DIVA",
    ):

    ex: appendCpCoe(mgr,'RP1',t0, xx)
    """

    # silence manager temporarily for quiet execution
    silenceMgr(mgr)

    newState = [
        M.Cartesian.x(xx[0]),
        M.Cartesian.y(xx[1]),
        M.Cartesian.z(xx[2]),
        M.Cartesian.dx(xx[3]),
        M.Cartesian.dy(xx[4]),
        M.Cartesian.dz(xx[5]),
    ]

    # create Control Point and add it to the problem
    cp = M.ControlPoint(
        mgr.boa,
        cpName,
        cpTime,
        body,
        center,
        frame,
        newState,
        mass,
        propagator,
    )
    mgr.cosmic.timeline().append(cp)

    if not fixCP:
        # Create new control list
        newControls = M.OptControlList(mgr.boa)

        lbound = -mmath.inf
        ubound = mmath.inf

        posBOX = 20*km
        velBOX = 0.01*km/s

        # # add controls in
        baseStr = "Cosmic/Cosmic/{0}/".format(cpName) #may not need this
        newControls.add(baseStr + "TIME", lbound, ubound)
        newControls.add(baseStr + "Cartesian.x", xx[0]-posBOX, xx[0]+posBOX)
        newControls.add(baseStr + "Cartesian.y", xx[1]-posBOX, xx[1]+posBOX)
        newControls.add(baseStr + "Cartesian.z", xx[2]-posBOX, xx[2]+posBOX)
        newControls.add(baseStr + "Cartesian.dx", xx[3]-velBOX, xx[3]+velBOX)
        newControls.add(baseStr + "Cartesian.dy", xx[4]-velBOX, xx[4]+velBOX)
        newControls.add(baseStr + "Cartesian.dz", xx[5]-velBOX, xx[5]+velBOX)

        # # add controls in
        # baseStr = "Cosmic/Cosmic/{0}/".format(cpName) #may not need this
        # newControls.add(baseStr + "TIME", lbound, ubound)
        # newControls.add(baseStr + "Cartesian.x", 0.9*x[0]*km, 1.1*x[0]*km)
        # newControls.add(baseStr + "Cartesian.y", 0.9*x[1]*km, 1.1*x[1]*km)
        # newControls.add(baseStr + "Cartesian.z", 0.9*x[2]*km, 1.1*x[2]*km)
        # newControls.add(baseStr + "Cartesian.dx", 0.9*x[3]*km/sec, 1.1*x[3]*km/sec)
        # newControls.add(baseStr + "Cartesian.dy", 0.9*x[4]*km/sec, 1.1*x[4]*km/sec)
        # newControls.add(baseStr + "Cartesian.dz", 0.9*x[5]*km/sec, 1.1*x[5]*km/sec)


        # Add the new controls list to the new control point
        mgr.cp[cpName].controls().append(newControls)

    # un silence manager
    unsilenceMgr(mgr)

    # print new CP
    #print(mgr.cp[cpName])
# ===========================================================================

# ===========================================================================
# append control point using Flyby params (as in Europa Clipper)
def appendCpFBs(
    mgr,
    cpName,
    cpTime,
    stQuery,
    mass=1000*kg,
    center=centerDefault,
    frame=frameDefault,
    body=scName,
    propagator="DIVA",
):
    """ append control point given in Classical Orbital Elements

    def appendCpCoe(
        mgr,
        cpName,
        cpTime,
        stQuery,
        mass=1000*kg,
        center=centerDefault,
        frame=frameDefault,
        body=scName,
        propagator="DIVA",
    ):

    stateParams{
        periAlt,
        bTheta,
        vInf,
        vInfDec,
        vInfRA,
        tra,
        }
    """

    # silence manager temporarily for quiet execution
    silenceMgr(mgr)

    #Read state from Query
    cpState= stQuery.state(cpTime).rotate(frame)

    newState = [
        # Seems redundant, but needed in order to make newState a 'Coordinate' list
        M.Conic.periapsisAltitude(M.Conic.periapsisAltitude(cpState)),
        M.Conic.bPlaneTheta(M.Conic.bPlaneTheta(cpState)),
        M.Conic.vInfinity(M.Conic.vInfinity(cpState)),
        M.Conic.inboundDec(M.Conic.inboundDec(cpState)),
        M.Conic.inboundRA(M.Conic.inboundRA(cpState)),
        M.Conic.trueAnomaly(M.Conic.trueAnomaly(cpState)),
    ]

    # create Control Point and add it to the problem
    cp = M.ControlPoint(
        mgr.boa,
        cpName,
        cpTime,
        body,
        center,
        frame,
        newState,
        mass,
        propagator,
    )
    mgr.cosmic.timeline().append(cp)

    # Create new control list
    newControls = M.OptControlList(mgr.boa)

    lbound = -mmath.inf
    ubound = mmath.inf

    # add controls in
    baseStr = "Cosmic/Cosmic/{0}/".format(cpName) #may not need this
    #newControls.add( "TIME", -3600 * s, 3600 * s)
    newControls.add(baseStr + "TIME", -86400 * s, 86400 * s)
    newControls.add(baseStr + "Conic.bPlaneTheta", -360*deg, 360*deg)
    newControls.add(baseStr + "Conic.periapsisAltitude", 40*km, 200*km)
    newControls.add(baseStr + "Conic.vInfinity", 0.2*km/sec, 1.0*km/sec),
    newControls.add(baseStr + "Conic.inboundDec", -360*deg, 360*deg)
    newControls.add(baseStr + "Conic.inboundRA", -360*deg, 360*deg)
    #newControls.add(baseStr + "Conic.trueAnomaly", lbound, ubound)

    # Add the new controls list to the new control point
    mgr.cp[cpName].controls().append(newControls)
    # un silence manager
    unsilenceMgr(mgr)

    #TODO:
    #- input stateParams as strings
    #- input control values -> {min,max,type: REL, ABS, FRAC}
    #   Check Units always


# ===========================================================================

# ===========================================================================
# Taken from Brian Anderson:
# ===========================================================================

# ===========================================================================
# add a time-based cosmic burn
def addTimeBurn(
    mgr,
    burnName,
    burnTime,
    frame="EMO2000",
    dmass=0 * kg,
    dvel=M.Dbl3Vec([1e-5] * 3),
    dtBound=3600*s,
    dvBound=0.02*km/s,
    silent=True,
):
    """add a CP-event based cosmic burn

    def addCpBurn(
        mgr,
        burnName,
        burnTime,
        frame="EMO2000",
        dmass=0 * kg,
        dvel=M.Dbl3Vec([1e-5] * 3),
        dtBound=3600*s,
        dvBound=0.02*km/s,
        silent=True,
    ):
    """

    # silence manager temporarily for quiet execution
    if silent:
        if not mgrIsSilent(mgr):
            silenceMgr(mgr)
            alreadySilent = False
        else:
            alreadySilent = True

    # create burn manager
    burnMgr = M.ImpulseBurnMgrBoa.read(mgr.boa, mgr.tl.body())
    burnBody = mgr.tl.body()

    # create new burn
    newBurn = M.ImpulseBurn(burnName, burnTime, frame, dmass, dvel)
    newOptBurn = M.OptBurnControl(mgr.boa, "IMPULSE", burnBody, burnName)
    newOptBurn.controls().add("TIME", -dtBound, dtBound)
    newOptBurn.controls().add("DX", -dvBound, dvBound)
    newOptBurn.controls().add("DY", -dvBound, dvBound)
    newOptBurn.controls().add("DZ", -dvBound, dvBound)

    # insert burn
    burnMgr.insert(newBurn)
    mgr.cosmic.burns().add(newOptBurn)

    # unsilence manager if it was temporarily silenced
    if silent and not alreadySilent:
        unsilenceMgr(mgr)


# ===========================================================================

# ===========================================================================
# add a CP-event based cosmic burn
def addCpBurn(
    mgr,
    burnName,
    cpName,
    tDelta=1 * s,
    frame="EMO2000",
    dmass=0 * kg,
    dvel=M.Dbl3Vec([1e-5] * 3),
    dvBound=0.02*km/s,
    silent=True,
):
    """add a CP-event based cosmic burn

    def addCpBurn(
    mgr,
    burnName,
    cpName,
    tDelta=1 * s,
    frame="EMO2000",
    dmass=0 * kg,
    dvel=M.Dbl3Vec([1e-5] * 3),
    dvBound=0.02*km/s,
    silent=True,
    ):
    """

    # silence manager temporarily for quiet execution
    if silent:
        if not mgrIsSilent(mgr):
            silenceMgr(mgr)
            alreadySilent = False
        else:
            alreadySilent = True

    # create burn manager
    burnMgr = M.ImpulseBurnMgrBoa.read(mgr.boa, mgr.tl.body())
    burnBody = mgr.tl.body()

    # create CP event for burn
    newEvent = M.CosmicEvent(mgr.boa, cpName, tDelta)

    # create new burn
    newBurn = M.ImpulseBurn(burnName, newEvent, frame, dmass, dvel)
    newOptBurn = M.OptBurnControl(mgr.boa, "IMPULSE", burnBody, burnName)
    newOptBurn.controls().add("DX", -dvBound, dvBound)
    newOptBurn.controls().add("DY", -dvBound, dvBound)
    newOptBurn.controls().add("DZ", -dvBound, dvBound)

    # insert burn
    burnMgr.insert(newBurn)
    mgr.cosmic.burns().add(newOptBurn)

    # unsilence manager if it was temporarily silenced
    if silent and not alreadySilent:
        unsilenceMgr(mgr)

# ===========================================================================
# ===========================================================================





# ===========================================================================
# TO ADD:
# -findCenterBodies()   #"Natural Center Bodies"
# - findDvsDiscon()
# --> Todo:  for dt in  M.TrajSetBoa.read(mgr.boa).intervals('mySC')
#                   check at each interval.
# - tList = Epoch.range(dt.begin(),dt.end(), timeStep)
# search for Natural Center Bodies
# Search for dv discontinuities

# Find dvDisc:
#   dv =(sq.state(t0+1*sec).vel() - sq.state(t0).vel()).mag()
#   relDV = (sq.state(t0+timeStep).vel() - sq.state(t0).vel()).mag()/sq.state(t0).vel().mag()

#  Try: just only VelMag differences:
# dvMag = (sq.state(t0+timeStep).vel().mag() - sq.state(t0).vel().mag())

# SCRIPT ---> :
# import matplotlib.pyplot as plt
# mgr.tl.createTraj(boa,mgr.problem)
# dt = M.TrajSetBoa.read(mgr.boa).intervals('mySC')[0]
# sq = TrajQuery( boa,'mySC','Earth','EME2000')
# timeStep = 10*sec
# tList =Epoch.range(dt.begin(),dt.end(), timeStep)
# dvMag_list = []
# for tt in tList[:-1]:
#     dvMag = (sq.state(tt+timeStep).vel() - sq.state(tt).vel()).mag()/sq.state(tt).vel().mag()
#     dvMag_list.append(dvMag)
# xx = range(0,len(dvMag_list))
# plt.plot(xx,dvMag_list)
# plt.show()
## <--|

# dvMag = (sq.state(tt+timeStep).vel().mag() - sq.state(tt).vel().mag())
# dvMag = (sq.state(tt+timeStep).vel() - sq.state(tt).vel()).mag()/sq.state(tt).vel().mag()

# ===========================================================================
