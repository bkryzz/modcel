# modcel <img width="200" src="https://github.com/user-attachments/assets/54fb86e4-f129-4bd0-8092-902baa056f78" />

Agent-based modeling of biological tissues 

<b>MODCEL v2.5 is an AGENT-BASED biophysical simulation model</b> aimed at mathematically reconstructing and modeling tissue physiology starting from a 2D image (for backward compatibility with versions ≤2.0, it also includes older modules for 3D simulation of e.g. spheroids or organoids, with central symmetry.)

  The main target application is to input 2D representations of
  biological tissues, as typically resulting from a histological
  section or a clinical biopsy. The user can either input a tissue map
  from a real sample (e.g., reconstructed via <a href="https://github.com/CellProfiler">CellProfiler</a>), or
  build their own tissue model, to imitate a realistic pattern or create an entirely
  arbitrary condition.

  The main module MODLOG treats a population of cells/agents that
  evolve in time and space according to a set of Markovian rules and 
  internal chemical evolution of a prescribed set of metabolites. 
  Chemical species can be assigned fixed, or transported and diffused by the fluids.
  Cell-cell and cell-fluid xchanges are regulated by a basic model of advection-diffusion
  equations.
  
  The following features decribe the cell/agent population:

  - <b>Discrete:</b> models an assembly of individual cells/agents living on a predefined lattice

  - <b>Continuous:</b> models Navier-Stokes fluid flow in blood and other 
  fluid capillaries (e.g., bile, lymph...) 

  - <b>Multistate:</b> each cell/agent evolves by probabilistic Markov chain through discrete states,
  e.g., from “healthy” to "diseased" to “dead”.

MODCEL is written with a mix of Python for the user-interaction parts
(creating and managing input files and cell configs, creating and
managing 2D maps and xy-plots of interesting quantities, etc...), and
Fortran90 for the computationally-intensive parts.

----- VERSION 2.0 was 2D/3D WITH SQ/CUB or VORONOI space. It included an attempt 
at OpenMP parallelization.

----- VERSION 2.4 included the <a href="https://github.com/julesghub/FEM2D">FEM2D Navier-Stokes solver</a> by
      J. Burkard to describe fluid flow in capillaries and ducts.

----- VERSION 2.4.2 opened to compatibility with <a href="https://github.com/opencor/opencor">OpenCOR</a> chemical models according to CellML standards.

----- Starting with version >2.5 the model has turned yet more general:

  1) You can select an arbitrary number of cell types:
     - the first nr_cel (ITYPE=1,2,3...) are solid cells (epithelial,
       hepato, adipo, etc,)
     - the next fl_cel (ITYPE=nr_cel+1, +2 ...) are fluid cells: can be
       used to define the capillary flows (blood, bile, lymph etc)
       Note that fl_cel can be =0, but nr_cel should be ≥1 (see below).
     - NUMF=SUM_nr_cel+SUM_fl_cel
     - NCEL=NUMF(total solid/fluid agents) + NDIS (boundary agents)
     - After building the first Voronoi map, all cells are set to
       parenchyma ITYPE=1 (whatever you chose it to represent...). If both nr_cel and fl_cel are =0, an empty Voronoi lattice is generated, to be filled by hand by putting down agents one by one.
     - Note that fluid cells are selected after the Voronoi map
       of the whole tissue has been constructed. This can be done
       by hand-picking with the mouse, or by asking the code to
       provide a random guess of vases distribution.
     - Usually, fl_cel=1 is venous blood, fl_cel=2 is arterial,
       fl_cel=3 can be bile or lymph, and so on...
       
<b>HOWEVER:</b> such cell definitions are just a mnemonic, not necessarily attached to a fixed real cell type. How a cell-agent actually behaves is fixed by the chemistry, see 2) below. 

<b>IN PRACTICE:</b> even if you can use a potentially illimited number of variants, try to use not more than 2-3 different cell types, and 1 or 2 fluids, otherwise the model could become very difficult to analyse and interpret. Make the chemistry as rich as possible, instead.

  2) You can select an arbitrary number of metabolites, with the FIXED convention that: 
     1=Insulin, 2=glucose, 3=oxygen, 4=FFA (free fatty acids).
     Then, the indices 5,6,7... can be used for any other metabolite.
     - Chemical equations for metabolite cycling (the 4 basic
       + any other) can be coded in CellML language. The easiest
       way to provide them is to use OpenCOR to script the eqs
       and then export them to Fortran; they will be linked by
       MODLOG and used normally. (Alternatively, you can code them yourself in the corresponding subroutine.)
     - Be sure to create a 'compartment' for each of the nr_cel
       and fl_cel cell types. For solid cells, you should define
       how metabolites are consumed, created and exchanged. For
       blood cells you must decide how nutrients, oxygen etc.
       are delivered, and/or how fluids (bile, lymph,...) are
       recovered from the solid cells and transported.

  3) Boolean logic rules determine how cells evolve, according to the summary of their chemical activity during time (a typical time step used in MODCEL simulations is 1 minute, adjustable with the parameter ITCONV). That is, how cells can go around their phases (G1, S, G2, M), enter the G0, move, duplicate, trans-differentiate, or differentiate (if they are stems), go into apoptose or necrose, and more. 


------------------------------------------------------------------
<b>To run the main program MODLOG requires a F90 compiler (GNU is ok) 
and a standard Python installation (v.3.7 or higher recommended).</b>

<b>The library GUIzero must be installed. On Windows/Mac just type:</b>
<i>pip3 install guizero</i>

<b>On Linux you may need to install tkinter first, e.g.:</b>
<i>sudo apt install python3-tk</i>

***

0) Always run your simulations in the main directory (i.e., you
must see there ./source/ ./models/ ./bin/ etc as subdirectories). 
At the very first use, run 'make install'. (NOTE: if you do a
'make clean', e.g. to get rid of old stuff, you will have to run
again 'make install'.)

1) MODCEL main program can be started with the script
                     ./modcel 

2) The first step is to provide the code with a tissue metabolic
model. This can be written in CellML format, or more easily can
be generated by OpenCOR. Run it through OpenCOR to generate the
corresponding Fortran code. 
MANDATORY: the models must be stored in the local ./models/ 
directory, where the ./modcel script will look for finding them.

3) The script recompiles the main program MODLOG by including your
chosen metabolic model. Then it opens a dialog window, where you
can specify the running parameters for a new simulation, or ask
the program to read back a previously stored configuration (this
is the contents of the binary file BACKUP; MODLOG keeps track of 
the last configuration in a file BACKUP_OLD whenever a new conf
is generated).

4) When the "start simulation" is pressed, the F90 code
   starts. It creates a new python module "voro.py" to:
   3.1) Build a new Voronoi 2D cell map
   3.2) Read an existing 2D cell map and modify it
   3.3) Read an histological section and obtain the
        corresponding cell map (experimental feature!)

   If the map includes fluids (e.g. blood sinusoids,
   bile ducts) the next step is to build triangular mesh
   for the Navier-Stokes solver. MODLOG understands it
   and creates the python module "triang.py", which will
   propose different 2D triangulations.

   <i>You will notice that in the console window, messages and comments coming from the Fortran code appear in</i> UPPERCASE, <i>while messages and comments from the Python modules appear in</i> lowercase.

6) Ensuite, the complete simulation can start

Different types of cells are included in 2.5 version. 
E.g. for the liver model with nr_cel=2 and fl_cel=3 they could be:
   ITYP=1 normal hepatocyte
   ITYP=2 diseased hepatocyte (fibrotic etc)
   ITYP=3 venous capillary
   ITYP=4 arterial capillary
   ITYP=5 bile duct,
plus the special "boundary agents" defining the interface between
the solid and fluid regions (e.g., including macrophages, stellate cells...)

-------------------------------------------------------
  Damage: each cell has a number of tags describing the 
  accumulation of individual defects (e.g., SSB, DSB, 
  base-excision, base-polymerization etc.)
  
----- In V1.0 only SSB and DSB damages to DNA were implemented

  Repair: each cell has an individual probability for 
  repairing the different defects.

  Inheritance: each cell has an individual probability 
  for duplication, the two sister cells inherit the DNA 
  of the mother cell.

  Geometry & topology: each cell in the model occupies a
  site of a regular lattice and “knows” its connectivity 
  to neighbor cells, allowing distributed effects
  
----- In V1.0 only 2D was implemented with diamond symmetry

----- In V2.0 the 2D/3D model introduced Voronoi tessellation

----- In V2.4 the 2D/3D model restored also square/cube symmetry

  Chemistry: reaction-diffusion equations describing the 
  flux of chemicals can be coupled to the spatial matrix 
  
----- In V1.0 only one tracer 2D field wass implemented

----- In V2.0 glucose/oxygen/activ/necro implemented

----- In V2.4 the interface to CellML was initiated

***

<b>BUILD AND COMPILE</b>

The main directory contains this README file, a ./source/ folder,
./bin/ folder where the executables are stored, ./models/ which
contains chemical examples, ./images/ for logos, a generic Makefile 
and a sample input file. 

Little modification is needed for the Makefile, usually only the
choice of the Fortran compiler needs to be updated for your 
local machine.

Once the 'make' command is issued in the main directory
compilation starts. 

Two executable files are created in the main directory:

- MODLOG, the executable code

- modmin, an utility to convert xyz files (for old versions, now
  replaced by Python utilities)

<b>RUNNING THE CODE</b>

MODLOG looks for an input file named 'input' in the same
directory from which it is launched. the sample input
provided should be enough self-explanatory. Starting from v2.4,
this input file is interactively generated by the GUI windows.

<b>ANALYZING THE SIMULATION</b>

<b>OUTPUT FILES</b>


ACCESSORY DIRECTORIES (for old versions before v2.0)

tcl_scripts contains some tcl scripts for treating graphical
            outputs

gnuplot_scripts contains scripts to accelerate the production
            of gnuplot elaborate graphics

LIST        contains a (not updated, sorry) list of the code
            variables and array names

---------------------------------------------------------------------------

The version 1.0 of MODLOG was at the origin of this work:
M. Tomezak, C. Abbadie, E. Lartigau, F. Cleri, <a href="https://www.sciencedirect.com/science/article/abs/pii/S0022519315005160">A biophysical model of cell evolution after cytotoxic treatments: Damage, repair and cell response</a>, J. Theor. Biol. <bf>389</bf>, 149-158 (2016)

The version 2.0 of MODLOG was used in this work:
F. Cleri, <a href="https://link.springer.com/article/10.1140/epje/i2019-11878-7">Agent-based model of multicellular tumor spheroid evolution including cell metabolism</a>, Eur. J. Phys. E: Soft Matter Biol. Phys. <bf>42</bf>, 112 (2019)
 
