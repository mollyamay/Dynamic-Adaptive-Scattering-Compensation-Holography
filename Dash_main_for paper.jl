"""
Scattering correction simulation 
 

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

"""

include("DASH_functions_for paper.jl")

## --------- user variables --------------

#choose the simulated method by activating the corresponding line below

method = "DASH"  
# method = "F-SHARP"  #see Papadopoulos et al. Nat. Phot. 2017
# method = "IMPACT"  #see Cui et. al. 
# method = "CSA"  #single pixel Basis
# method = "PA"   #partitioning algorithm (see Vellekoop et al. 2008; works with sets of random phase pattern)
# method = "GA"   #Genetic Algorithm (see Conkey et. al. Opex 2012)
    gen = 300 # for GA: no. of generations 
    N_gensize = 50 #for GA: number of patterns in each generation 
    R0 = [0.02, 0.01] #for GA: initial and end fraction of mutatable pixels 
    λ = 10  #for GA: decay rate of mutation rate
    pow = 5 #for GA: exponent of probability-law to select parents

iter = 10 #number of algorithm iterations (does not apply to GA -> here the variable "gen" takes the role of the iteration no.)
cyc = 1 #no of simulated experiments cycles to average over
t_meas = 200e-6 #time per measurement in seconds (affects the number of photons accumulated in each measurement)
bg = 0 #background counts per second

N = 16 #side length of SLM in pixels
infl = 2 #inflation factor -> together with N determines the number of scattering pixels. The number of scattering pixels is (infl * N)^2

f = 0.3  #energy fraction going into the scan beam
type_specimen = "layer"  #"layer" or "bead" or a number between 0 and 100 (in this case, the number is interpreted as a percentage of pixels which are "on")


#------Execution-------------

B = basics(N, infl, type_specimen = type_specimen)
RUN(B,method, gen, N_gensize, R0, λ, pow, iter, cyc, t_meas, bg, N, infl, f, type_specimen)















