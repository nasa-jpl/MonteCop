#===========================================================================
#
# Section 392 Navigation and Mission Design
#
# Copyright (C) 2018, California Institute of Technology.
# U.S. Government Sponsorship under NASA Contract NAS7-03001 is acknowledged.
#
#===========================================================================

 monteCop: a Monte-Copernicus interface. 

 Description: monteCop is python module composed by a set of scripts that allows the transfer of
 		Copernicus solutions into COSMIC timelines (COSMIC is the Monte Optimizer module),	
		as well as to convert Monte trajectories ( in form of *.bsp, *.boa, *.py files)	
		into Copernicus idecks.  Additional scripts to convert generic trajectories in form of SPK kernels
		(*.bsp files) into Cosmic timelines and Copernicus idecks are also available.  
		
Scripts:	
		bsp2cosmic.py:  	converts a *.bsp file into a Cosmic timeline
		bsp2ideck.py:	  	converts a *.bsp file into a Copernicus ideck
		bsp2visualCol.py: 	converts a *.bsp file into a Copernicus ideck (Visualization only)
		cosmic2ideck.py: 	converts a Cosmic timeline into a Copernicus ideck
		csv2ideck.py		converts a *.csv file into a Copernicus ideck
                 

**Note: The user need to have robocoppy, the Copernicus Python interface, in the python path.  


 by: Ricardo L. Restrepo. NASA/JPL 392M
       ricardo.l.restrepo@jpl.nasa.gov		

		
 version:    1.0.1
 Last edit: Nov 15, 2023		
		
#==============================================================
#==============================================================
