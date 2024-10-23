#!/usr/bin/env mpython_q

#===========================================================================
#
# Section 392 Navigation and Mission Design
#
# Copyright (C) 2020, California Institute of Technology.
# U.S. Government Sponsorship under NASA Contract NAS7-03001 is acknowledged.
#
#===========================================================================

"""
 Plot r (distance from central body) vs time for cosmic timeline
"""

#__author__  = 'Ricardo Restrepo (392M)'

#===========================================================================
# Place all imports after here.
#
import argparse
from mpy.units import s, sec, minute, hour, m, km, day, deg, rad, kg

import mpylab
import matplotlib.pyplot as plt

import Monte as M
from mpy.opt.cosmic import Manager

import math

# Place all imports before here.
#===========================================================================

# Params Set up:
#center = 'Enceladus'
#BodyID = "encnf5"
#frame = "IAU Enceladus Fixed"

# Create and populate the parser
parser = argparse.ArgumentParser(description = ' Plot r (distance from central body) vs time for cosmic timeline ')
parser.add_argument("inputFile", metavar="input.py", help="Cosmic input file to be used")
parser.add_argument('-c', '--Center', default= "Enceladus", help = "Body Center. Default: 'Enceladus' ")
parser.add_argument('-f', '--Frame', default='IAU Enceladus Fixed', help = "Ploting Frame. Ex = 'EMO2000' ")
parser.add_argument('-sc', '--scName', default="encnf5", help = "Spacecraft Name. Default ='encnf5' ")
parser.add_argument("-to","--tiniOffset", default="0", help="t0 offset- days. Default='' (use t0 from bsp)")
parser.add_argument("-tl","--tlInterval", default="0", help="Trajectory time line interval from t0_offset. Default='' (use bsp time span)")
parser.add_argument('-dt', '--timeStep', default=10, help = "Step time for plotting (in minutes). Default: 10 (min)" )
#parser.add_argument('-s','--savePlot', action='store_true', help='save plots')

# TODO: Add Epoch Calendar format: ...

args = parser.parse_args()

scName = args.scName
bodyCenter = args.Center
plotFrame = args.Frame

# TODO:
# Add optional Geographic, etc.
# Move to encLibrlr.py

args = parser.parse_args()

# raise Exception('raise')

M.Epoch.setFormat('full')

# Load Cosmic input
mgr = Manager(boa)
mgr.loadInput(args.inputFile)
# Propagate the traj to get the burns in the "right places"
mgr.tl.createTraj(mgr.boa, mgr.problem)

# Define t0:
if args.tiniOffset: t0_offset = args.tiniOffset
t0_offset =  float(args.tiniOffset)*day  # Default = 0
t0 = mgr.tl.interval().begin() + t0_offset
# Define tf:
tl_interval = float(args.tlInterval)*day  # Default = 0: If == 0 -> use tl_end()
if not tl_interval:
    tf = mgr.tl.interval().end()
else:
    tf = t0 + tl_interval
# Define dt:
dt = 1*minute
tt_arr = Epoch.range(t0,tf,dt)

# Create queary:
#query  = TrajQuery( boa, scName, bodyCenter, 'IAU ' + bodyCenter +' Fixed' )
query  = TrajQuery( boa, scName, bodyCenter, plotFrame )

# ----------------------------------------------------------------------------
# --> Find min range points:
# Search Paramaeters:
srchInterval = TimeInterval(t0,tf)
rpSrchStep = 10*minute
raSrchStep = 10*minute

ev = ApsisEvent( TrajQuery( boa, scName, "Enceladus" ), ApsisEvent.PERIAPSIS )
rp_srch = ev.search(srchInterval, rpSrchStep )
rp_arr = [r.value() for r in rp_srch]
rp_tarr = [r.time() for r in rp_srch]


ev = ApsisEvent( TrajQuery( boa, scName, "Enceladus" ), ApsisEvent.APOAPSIS )
ra_srch = ev.search(srchInterval, raSrchStep )
ra_arr = [r.value() for r in ra_srch]
ra_tarr = [r.time() for r in ra_srch]

#Find APOs at apo:
aop_arr=[]
ran_arr=[]
for tt in rp_tarr:
    aop_arr.append(Conic.argumentOfPeriapsis(query.state(tt)))
    ran_arr.append(Conic.longitudeOfNode(query.state(tt)))

# Plot Ra/Rp: Separate figures

# fig1, ax1 = plt.subplots()
# ax1.plot(rp_tarr, rp_arr)
# ax1.set(ylabel='Rp (km)')
# ax1.grid()
# fig2, ax2 = plt.subplots()
# ax2.plot(ra_tarr, ra_arr)
# ax2.set(ylabel='Ra (km)')
# ax2.grid()
# mpylab.show()

fig0, ax0 = plt.subplots()
# ax1.plot(rp_tarr[23:], rp_arr[23:],'-ro')
# ax1.plot(ra_tarr[23:], ra_arr[23:],'-bo')
ax0.plot(rp_tarr, rp_arr,'-ro', label = 'Rp')
ax0.plot(ra_tarr, ra_arr,'-bo', label = 'Ra')

# ax0.scatter(rp_tarr, rp_arr, c='r', label='Rp')
# ax0.scatter(ra_tarr, ra_arr, c='b', label='Ra')
ax0.legend()
ax0.set(ylabel='Ra/Rp (km)')
ax0.grid()


#New Plots: filter rp.
rpIntv = 6
RpId = 7
iniCount =  0

rp_arrNew = [rp for ii,rp in enumerate(rp_arr) if ii >= iniCount and (ii+RpId)%rpIntv == 0 ]
rp_tarrNew = [tt for ii,tt in enumerate(rp_tarr) if ii >= iniCount and (ii+RpId)%rpIntv == 0 ]
aop_arrNew = [aop for ii,aop in enumerate(aop_arr) if ii >= iniCount and (ii+RpId)%rpIntv == 0 ]
ran_arrNew = [ran for ii,ran in enumerate(ran_arr) if ii >= iniCount and (ii+RpId)%rpIntv == 0 ]

print(' ... No or rp selected: '+ str(len(rp_arrNew)))

fig2, ax2 = plt.subplots()
#plot Peris All and sel;ected:
ax2.plot(rp_tarr, rp_arr,'-r.')
#ax2.scatter(rp_tarr, rp_arr,'-ro', marker = ".")
ax2.scatter(rp_tarrNew, rp_arrNew, c='g', label='Rp')
ax2.legend()
ax2.grid()
fig3, ax3 = plt.subplots()
ax3.scatter(rp_tarrNew, aop_arrNew, c='r', label='AOP')
ax3.legend()
ax3.grid()
fig4, ax4 = plt.subplots()
ax4.scatter(rp_tarrNew, ran_arrNew, c='b', label='RAN')
ax4.legend()
ax4.grid()

#mpylab.show(block=False)
mpylab.show(block=True)


print(' ... TODO: UPDATE with bsp2cosmic_RpCP.py -> Add plots')


#mpylab.close('all')

# ----------------------------------------------------------------------------
