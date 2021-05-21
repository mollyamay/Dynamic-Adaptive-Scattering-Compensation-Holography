# F-sharp
"""
Created on Fri Jul 31 16:33:42 2020

This file simulates the basic working principle of DASH
A random phase scatterer is assumed, which is optically in the same plane ("conjugate") as the correction device (SLM)
Just execute this Python file. 
If desired, some parameters can be changed in the section below "user parameters"


@author: q004aj
"""

#importing required packages
from numpy import * 
from matplotlib.pyplot import *
import random

#%% user parameters

# method = "F-SHARP"  #choose either "DASH" or "F-SHARP"
method = "F-SHARP"

NN = 16 #side length of the SLM in pixels; the number of correctable modes is thus N x N
N = 32 #number of scattered modes = N^2

iter = 10 #no. of iterations, i.e. full mode cycles

eta = 100 #efficiency of 2-photon signal generation; vary e.g. between 1 and 100

sampletype = 'layer'  #sample can be a point source: "bead" or a 2D fluorescent plane: "layer"

trials = 5 #run multiple trials to find standard deviation
all_data = zeros([NN**2*iter, trials])
    
#%% calculating plane-wave modes / definition of required functions

N_modes = NN*NN
M = np.zeros((N, N, N_modes)) #initializing array of modes
k = fft.fftfreq(NN, 1/NN)
Kx, Ky = np.meshgrid(k,k) #pixel coordinates

mm = 0
for m in range(NN):      
    for n in range(NN):            
        gx = -pi + (n+1)*2*pi/NN #grating vector in x
        gy = -pi + (m+1)*2*pi/NN #grating vector in x
        M[:,:,mm] = np.resize((gx*Kx + gy*Ky), (N,N))
        mm += 1

def TPEF(E_in, scat, sample): #two photon fluorescence signal
    "calculation of two photon signal based on pupil field"
    I2ph = int(eta*1e3*sum(np.abs(sample * (np.fft.fft2(scat * E_in)/N/N))**4))
    return np.random.poisson(I2ph)



#%% EXECUTE OPTIMIZATION
for idx in range(trials):
    scat = np.exp(1j*2*pi*np.random.rand(N,N)) #create a random phase scatterer
    theta = array([0, 2*pi/3, 4*pi/3]) #reference phase values
    f = 0.3 #energy fraction going into the testing M
    
    #initializations
    E_in = ones((N,N), dtype="complex")
    I2ph = zeros((3))
    a = zeros(N_modes, dtype = "complex") #initializing mode correction amplitudes
    C = zeros((N,N), dtype = "complex") #initializing correction mask
    signals = zeros((iter, N_modes))
    
    if sampletype == 'layer':
        sample = np.ones((N,N)) #fluorescent plane
    elif sampletype == 'bead':
        sample = zeros((N,N)); sample[0,0] = 1;    #single beacon
    
    for i in range(iter):
        
        print(iter-i)
        
        for m in range(N_modes): #mode-stepping loop
             
            for p in range(size(theta)): #phase-stepping loop
                           
                if method == "F-SHARP":
                    E_SLM = sqrt(f) * exp(1j*M[:,:,m] + 1j*theta[p]) + sqrt(1-f) * exp(1j*angle(C))   #F-SHARP overlaps two "real" beams 
                elif method == "DASH":
                    E_SLM = exp(1j * angle(sqrt(f) * exp(1j*M[:,:,m] + 1j*theta[p]) + sqrt(1-f) * exp(1j*angle(C)))) #in DASH, the two beams are created by a single phase-only SLM
                
                I2ph[p] = TPEF(E_SLM, scat, sample) #get 2-photon signal
            
            signals[i, m] = np.mean(I2ph) #mean signal over all phase steps
            a[m] = np.sum(sqrt(I2ph) * exp(+1j*theta)) / size(theta) #retrieving a = |a|*exp(1j*phi_m), we multiply with exp(+1j*theta) instead of exp(-1j*theta), because we want to directly calculate the correction phase )
                
            if method == "DASH":
                C += a[m] * exp(1j*M[:,:,m]) #for DASH, immediately update the correction mask
            
        if method == "F-SHARP":
            C = np.sum(a * exp(1j * M), 2) #for F-SHARP, a correction pattern is calculated after every full mode iteration
    all_data[:,idx] = ravel(signals)

#-------------- display results------------------

#display singal trend
figure(1)        
plot(mean(all_data, axis = 1))
xlabel('measurement no.')
ylabel('2 photon signal / photons')
grid(1)        
title(method + ", " + str(N_modes) + " modes")

#display focus
# =============================================================================
# I_scat = np.abs(fft.ifftshift(fft.fft2(scat)))**2
# I_corr = np.abs(fft.ifftshift(fft.fft2(scat * exp(1j*angle(C)))))**2
# 
# figure(2)
# imshow(I_scat, cmap = "jet")
# title("scattered focus irradiance")
# colorbar()
# 
# figure(3)
# imshow(I_corr, cmap = "jet")
# title("corrected focus irradiance")
# colorbar()
# =============================================================================

print("minimum / maximum signals per measurement: " + str(np.int(np.min(signals))) + " / " + str(np.int(np.max(signals))) + " photons")



 