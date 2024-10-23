# MonteCop
The module MonteCop consists of a set of python scripts that serves as an interface to transfer spacecraft trajectory solutions between different NASA flight mechanics tools, in particular between Monte (JPL ) and Copernicus (JSC). The module is accompanied by a set of use cases, including a detailed description per case (See monteCop_UseCases folder). 

These products derived from an NESC (NASA Engineering and Safety Center) effort called “Flight Mechanics Analysis Tools Interoperability and Component Sharing”, described in the technical report, NASA/TM-20230006507 NESC-RP-18-01313.  

About Flight Mechanics Analysis Tools Interoperability and Component Sharing:
Several NASA centers have developed independent flight mechanics tools to meet the science needs of missions. This assessment sought to explore the ways to increase the interoperability of three specific tools: Copernicus from Johnson Space Center (JSC), the General Mission Analysis Tool (GMAT) from Goddard Spaceflight Flight Center (GSFC), and the Mission-Analysis Operations Navigation Toolkit Environment (MONTE) from the Jet Propulsion Laboratory (JPL). All of these tools are utilized in various spaceflight regimes and mission lifecycles (e.g., design, analysis, operations) to generate a variety of products (e.g., maneuver planning, orbit determination, error analysis, trade study, flight products). Each tool was built over the years with specific goals unique to each center, which were based on the science missions they supported, and naturally lead to variations in their capabilities. Before this assessment, these tools were not integrated and could not easily share data, models, or components. The goal of the assessment was to improve interoperability and component sharing of these flight mechanic tools to increase Agency efficiency and reduce cost.

# Monte-Copernicus interface. 
 Description: monteCop is python module composed by a set of scripts that allows the transfer of Copernicus solutions into COSMIC timelines (COSMIC is the Monte Optimizer module),  as well as to convert Monte trajectories ( in form of *.bsp, *.boa, *.py files)	into Copernicus idecks.  Additional scripts to convert generic trajectories in form of SPK kernels (*.bsp files) into Cosmic timelines and Copernicus idecks are also available.  
		
Scripts:  
          *bsp2cosmic.py*:  	converts a *.bsp file into a Cosmic timeline    
          *bsp2ideck.py*:	converts a *.bsp file into a Copernicus ideck    
          *bsp2visualCop.py*: 	converts a *.bsp file into a Copernicus ideck (Visualization only)    
          *cosmic2ideck.py*: 	converts a Cosmic timeline into a Copernicus ideck    
          *csv2ideck.py*:		converts a *.csv file into a Copernicus ideck    
                 

*Note: The user needs to have robocoppy, the Copernicus Python interface, in its python path.  

# Uses Cases
This folder contains the use cases and supporting data files described in the technical report, NASA/TM-20230006507 NESC-RP-18-01313 “Flight Mechanics Analysis Tools Interoperability and Component Sharing”:

**HLS**   
Problem Description: This use case represents the Copernicus to MONTE data transfer scenario of a Low Lunar Orbit (LLO) to a Near Rectilinear Halo Orbit (NRHO) trajectory. During the Artemis III mission, one of the phases will involve returning the astronauts from the south pole of the Moon to the Orion capsule at the NRHO Gateway orbit.

**Monte 2 Copernicus Visualization**  
Problem Description: This use case represents the SPK to Copernicus data transfer scenario where complex trajectories developed in Monte or any other tool, are visualized in Copernicus.  The python script used to generate trajectory visualizations in Copernicus from SPK kernels is called bsp2visualCop.py. The script follows the trajectory recovery process outlined in the previous section. Here, the trajectory contained in the SPK kernel is loaded into Copernicus as an ephemeris file and a segment attached to it as a “static point trajectory” is used for visualization purposes (no optimization is applied here). The python function call requires as input the SPK kernel (.bsp) file. Optional argument inputs, as the spacecraft ID, the central body, the frame to visualize the trajectory are possible. For more complex trajectory visualizations, an optional JSON input file can be passed to the function call to generate a more personalized output with specific frames, central body, number of segments, colors, etc. 

**Icy Moons Tour** 
Problem Description: This use case represents the SPK to Monte/Cosmic data transfer scenario where trajectories developed in Copernicus or any other astrodynamics capable to store trajectories in form of SPK kernels, are recreated in Monte.  The use case considered is the design of Moon tours to explore the Icy worlds.  A Saturn Moon tour that enables an Enceladus capture orbit is used.   A traditional tour consists of a sequence of flybys of the moons Titan, Rhea, Dione, Tethys, and Enceladus to pump down the energy of the initial Saturn post-capture orbit, enabling an affordable ∆V Enceladus capture insertion. For simplicity, in the consider example, only the Enceladus phase, the last part of the tour, es considered, although the process of generating other phases of the tour is similar.  For details refer to Trajectory Reverse Engineering paper†.

**Enceladus Orbiter** 
Problem Description: The following use case demonstrates the use of the Monte-Copernicus interface (Trajectory Reverse Engineering strategy) for two orbiter trajectories around Enceladus: an NRHO and an Enceladus Low Orbit (ELO). Enceladus is a celestial body known for its complex gravitational environment and designing trajectories under this environment is particularly challenging.  The initial trajectories are generated using the Circular Restricted Three Body Problem (CR3BP) model in Copernicus. The Copernicus solutions serves as initial condition for reconstruction in Monte under a higher fidelity dynamical model.  For details refer to Trajectory Reverse Engineering paper†. 

**Europa Flyby**  
Problem Description: This use case represents a direct trajectory transfer strategy from MONTE to Copernicus, using as a test case a Europa Clipper flyby sequence. The script used is cosmic2cop.py. The steps taken for this use case are to 1) Design the trajectory in MONTE in form of a Cosmic Timeline 2) identify all the relevant control points (e.g., maneuvers) and events (e.g., apoapsis/periapsis), 3) convert the trajectory into a Copernicus file.  4) utilize Copernicus’ GUI to adjust the segments and/or constraints/functions, optimizer, etc. To generate an ideck from a Cosmic timeline, a similar process to the one implemented for bsp2cosmic.py is performed, but the multiple-shooting strategy is manually set up on the python script by creating segments that are propagated forward and backward on time while imposing continuity true state constraints.  For details refer to Trajectory Reverse Engineering paper†.


 by: Ricardo L. Restrepo. NASA/JPL 392M
     ricardo.l.restrepo@jpl.nasa.gov		
