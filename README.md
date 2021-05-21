Source Data:
Raw data for figures in Nature Communications publication NCOMMS-21-03415A.

Julia Files:

These files simulates the performance of DASH, F-SHARP, IMPACT, a genetic algorithm, a partitioning algorithm, and a continuous sequential algorithm usign Julia version 1.5.3, instructions for inatallation can be found at www.julialang.org. 
Installation time is generally less than five minutes. All simulations are run from the DASH_main file, with companion functions in the DASH_functions file.

In the simulation, a random phase scatterer is assumed, which is optically in the same plane ("conjugate") as the correction device (SLM).
To run the simulation, just execute the entire Julia file. Execution should take less than 1 minute and the behavior shown in Fig. 1a of the paper will automatically be plotted.
If desired, some parameters can be changed in the section "user variables".

simulation of the following methods: 
DASH    -   our method introduced in this paper
F-SHARP -   introduced by Papadopoulos et al. Nat. Phot. 2017
IMPACT  -   introduced by Tang et al. PNAS 2012, the version here is somewhat modified: we use a smaller number of total measurements because we perform phase-stepping (see suppl. document)
CSA     -   continuous sequential algorithm, see Vellekoop & Mosk Opt. Commun 2008 
PA      -   partitioning algorithm, see Vellekoop & Mosk Opt. Commun 2008 
GA      -   genetic algorithm, see Conkey et al. Opt. Express 2012

Assumptions: 
-) random 2D phase scatterer with (N*infl)^2 pixels in a squared pupil.
-) SLM with N^2 pixels, conjugated to the objective pupil and phase scatterer
-) 2D fluorescent sample; structure can be chosen ("layer", "bead" (a single pixel is fluorescent) or multiple beads at random positions in the focal plane)


Python File: 

This file simulates the basic working principle of DASH usign Python version 3.9.0, instructions can be found at www.python.org. 
Installation time is generally less than five minutes.

In the simulation, a random phase scatterer is assumed, which is optically in the same plane ("conjugate") as the correction device (SLM).
To run the simulation, just execute the entire Python file. Execution should take less than 1 minute and the behavior shown in Fig. 1a of the paper will automatically be plotted.
If desired, some parameters can be changed in the section "user parameters".
