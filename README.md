<img width="130" src=https://github.com/user-attachments/assets/1e76c787-8494-41d4-93ee-a799ae8675c1 /> 
<!img width="200" src="https://github.com/user-attachments/assets/54fb86e4-f129-4bd0-8092-902baa056f78" />


Agent-based modeling of biological tissues (from the Latin <i><b>MOD</b>ulamen <b>CEL</b>lularum</i>, or "the rhythm of cells").

<b>MODCEL v2.5 is an AGENT-BASED biophysical simulation model</b> aimed at mathematically reconstructing and modeling tissue physiology, possibly starting from 2D images (for backward compatibility with versions ≤2.0, it also includes older modules for 3D simulation of e.g. spheroids or organoids with central symmetry.)

  The key intended target application is to input 2D representations of
  biological tissues, as typically resulting from a histological
  section or a clinical biopsy. The user can either input a tissue map
  from a real sample (e.g., slides reconstructed via segmentation tools like <a href="https://github.com/MouseLand/cellpose">CellPose</a> or <a href="https://github.com/CellProfiler">CellProfiler</a>), or
  build their own tissue model, to imitate a realistic pattern or create an entirely
  arbitrary condition.

  The main module MODLOG (its number-crunching F90 engine) treats a population of cells/agents that
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
      J. Burkard to describe fluid flow in capillaries and ducts. In the most recent version this part has been now superseded by the <i>multigrid</i> structure and the related OpenMP-based solvers.

----- VERSION 2.4.2 opened to compatibility with <a href="https://github.com/opencor/opencor">OpenCOR</a> chemical models according to CellML scripting. We plan to include a more general parser to read more various formats (e.g. from COPASI).

----- Starting with version >2.5 the model has turned yet more general: arbitrary number of cell types, metabolites, border agents, Boolean logic, and more.

---------------------------------------------------------------------------

The version 1.0 of MODLOG was at the origin of this work:
M. Tomezak, C. Abbadie, E. Lartigau, F. Cleri, <a href="https://www.sciencedirect.com/science/article/abs/pii/S0022519315005160">A biophysical model of cell evolution after cytotoxic treatments: Damage, repair and cell response</a>, J. Theor. Biol. <bf>389</bf>, 149-158 (2016)

The version 2.0 of MODLOG was used in this work:
F. Cleri, <a href="https://link.springer.com/article/10.1140/epje/i2019-11878-7">Agent-based model of multicellular tumor spheroid evolution including cell metabolism</a>, Eur. J. Phys. E: Soft Matter Biol. Phys. <bf>42</bf>, 112 (2019)

A variant thereof was used in this work:
L. Terrassoux et al., <a href="https://onlinelibrary.wiley.com/doi/full/10.1002/smll.202505343">Novel Diffuse Midline Glioma-on-Chip Recapitulating Tumor Biophysical Microenvironment to Assess the Heterogeneity of Response to Therapies</a>, Small e05343, (2025)
 
