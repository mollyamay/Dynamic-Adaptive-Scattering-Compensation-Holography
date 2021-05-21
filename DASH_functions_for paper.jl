"""DASH functions"""

#loading packages
begin 

    using FFTW
    using PyPlot
    pygui(true)
    using Distributed
    #using ImageView
    #using Plots
    #gr()
    #pyplot() #backend
    #default(show = true) #for Plots
    using Noise
    using SparseArrays
    using Random
    using Statistics
    using StatsBase
    #using Images
    using LinearAlgebra
    using LsqFit
    using Images

end

#---define strcture containing basic constants/quantities
struct basics
    N::Int
    infl::Int
    N2::Int
    type_specimen::String
    NL::Int
    pupil::Array{Bool,2}
    specimen::Array{Any,2}
    x
    η
    FT

    function basics(N, infl; type_specimen = "layer")
        N2 = infl * N
        x = fftfreq(N,N)
        η = 5e4 #efficiency of two-photon process
        FT = plan_fft!(0im .+ ones(N2,N2))
        NL = 2 #nonlinearity
        pupil = ones(N,N)

        if type_specimen == "layer"
            specimen = ones(N2,N2)
        elseif type_specimen == "bead"
            specimen = zeros(N2,N2)
            specimen[1] = 1
        else #if the input is a number, it is interpreted as a percentage of pixels which act as fluorescent points
            idx = sample(1:N2*N2, Int(trunc.(parse(Int, type_specimen)/100 * N2*N2)), replace = false)
            specimen = zeros(N2,N2)
            specimen[idx] .= 1
        end

        new(N, infl, N2, type_specimen, NL, pupil, specimen, x, η, FT)
    end

end


##---defining some functions and fit models ---------------


"""
create a scatterer
"""
#random white noise
function scatterer(N2::Int)
    2π * rand(N2,N2)
end

#random with gaussian envelope
function scatterer(N2::Int, σ::Real)
    x = fftfreq(N2,N2)
    w = [exp.(-(x[m]^2 + x[n]^2)/2/σ^2) for m = 1:N2, n = 1:N2]
    angle.(fft(w .* exp.(1im*2π*rand(N2,N2))))
end

#with Zernikes: a = vector of zernikes: 1st column: mode no. 2nd column: magnitudes
function scatterer(N2::Int, a::Array{Float64,2})
    x = fftfreq(N2,2)
    Z = evaluateZernike(x, Int.(a[:,1]), a[:,2], index=:Noll)
end

"""
calculate two-photon signal
"""
function TPEF(B, scat::Array{Float64,2}, holo::AbstractArray{Complex{Float64},2}, t_meas::Float64; bg = 0, noise = true)
    

    If2 = abs.(B.specimen .* (B.FT*(exp.(1im*scat) .* inflate(B.pupil .* holo, B.infl))/B.N2^2)).^(B.NL^2)
    
    #-----activate next 2 lines for denser spatial sampling in the sample plane------
    # infl_res = 4
    # If2 = abs.(inflate(B.specimen, infl_res) .* fft(fftshift(embed(ifftshift(inflate(B.pupil .* holo, B.infl) .* exp.(1im*scat)), infl_res*[B.N2, B.N2])))/(B.N2)^2).^4 
    
    if noise
        PMT = poisson(B.η^2*[sum(If2)*t_meas] .+ bg)[1] 
    else
        PMT = B.η^2 * sum(If2) * t_meas
    end
    return PMT, If2
end

"""
evaluation of a_n and phi_n
"""
@. model_cos(x, p) = p[1] + p[2]*cos(x + p[3])    #@. --> to execute on each ROW of the data -> this means that x and p must be COLUMS 

function get_phase(PMT, ϕ_R; method = "simple")
    
    P = length(PMT)
    
    if method == "simple"
        a = sum(PMT.^(1/B.NL) .* exp.(-1im*ϕ_R))/P

    elseif method == "fit"
        p0 = [1, 0.1, 0]  
        fit = curve_fit(model_cos, ϕ_R, I2ph.^(1/B.NL), p0)
        (stderror(fit)[3] > pi/5) ? E_ϕ = π : E_ϕ = stderror(fit)[3]
        c = coef(fit)[2] .* exp(1im*coef(fit)[3])

    end

end

"""
2D Gaussian fit
"""
function model_gauss(lincoord, p)
    N = Int(sqrt(length(lincoord)))
    x = (lincoord .-1).%N .+1 #row coordinate
    y = ceil.(lincoord./N) #column coordinate
    p[1] .+ p[2].*exp.(-(x .- N/2 .- p[3]).^2 ./ 2p[5]^2 .- (y .- N/2 .- p[4]).^2 ./ 2p[6]^2)
end

    ##inflating a matrix by an integer factor
"""inflating a matrix by an integer factor
inflate!(E_small, E_large, infl::Int)
"""
function inflate!(B::AbstractArray, A::AbstractArray, n::Int)     
    niA, njA = size(A)   
    niB, njB = size(B)
    @assert niA == niB*n
    @assert njA == njB*n
    @inbounds for j in 1:njB
        for i in 1:niB
            for k in 1:n
                for l in 1:n                    
                    A[n*(i-1)+k,n*(j-1)+l] = B[i,j]
                end
            end
        end
    end    
    A
end

"""inflating a matrix by an integer factor
inflate(E_small, E_large, infl::Int)
"""    
function inflate(B::AbstractArray, n::Int)
    niB, njB = size(B)    
    A = similar(B, (niB*n, njB*n))
    inflate!(B, A, n)
end

"""
embedding a 2D array in to a larger canvas of zeros
"""
function embed(E_in, N_target)

    #---test inputs
    #   E_in = ones(30,30)
    #   N_target = (29, 29)
    
    N = size(E_in)

    if ndims(E_in) == 2   #if output should be a 2D array

        ΔN = N_target .- N

        #process dimension 1
        if ΔN[1] > 0  #padding
            pad_top = Int64(ceil((ΔN[1]-1)/2))
            pad_bottom =  Int64(floor((ΔN[1]+1)/2))
            tmp = padarray(E_in,Fill(0,(pad_top,0),(pad_bottom,0)))
        else #cropping
            idx = Int64(floor(-ΔN[1]/2)) + 1  #start index
            tmp = E_in[idx .+ (1:N_target[1]),:]
        end

        #process dimension 2
        if ΔN[2] > 0 #padding
            pad_left =  Int64(ceil((ΔN[2]-1)/2))
            pad_right =  Int64(floor((ΔN[2]+1)/2))
            E_out = collect(padarray(tmp,Fill(0,(0,pad_left),(0,pad_right))))
        else #cropping
            idx = Int64(floor(-(ΔN[2])/2)) + 1
            E_out = tmp[:, idx .+ (1:N_target[2])]
        end

    end

return E_out
end


##-------------- main DASH function -----------------------------------------------

"""general scattering correction function; uses phase stepping""" 
function SCATCORR(B, scat::Array{Float64,2}, t_meas::Float64; P::Int = 3, method::String="IMPACT", f::Float64=0.25, iter::Int=5, update::Int=1, bg::Int=0)

    #-----choosing basis-------
    if method == "DASH" || method == "F-SHARP" 
        N_modes = B.N^2
        type = "complex"
        
        #-----plane wave modes-----
        modes = angle.(fft(reshape((I(B.N^2)), B.N^2, B.N, B.N), (2,3))) #creating plane wave modes
        freq = [B.x[m]^2 + B.x[n]^2 for m = 1:B.N, n = 1:B.N]# spatial frequencies of modes
                 
    elseif method == "IMPACT"
        N_modes = Int(sum(B.pupil))
        type = "complex"
        freq = (0:N_modes-1)/N_modes #frequencies of all pixels
        ri, ci, _ = findnz(sparse(B.pupil))
        modes = zeros(N_modes, B.N, B.N)
        ind = shuffle(1:N_modes)
        for t = 1:N_modes, m = 1:N_modes 
            modes[t, ri[ind[m]], ci[ind[m]]] = 2π*freq[m]*t
        end
   
    elseif method == "CSA" #continuous sequential algorithm
        N_modes = Int(sum(B.pupil))
        type = "phase"
        modes = reshape((I(B.N^2)), B.N^2, B.N, B.N)
   
    elseif method == "PA" 
        N_modes = Int(sum(B.pupil)) * iter #for GA, this is the pupulation size 
        iter = 1 #only one iteration is measured
        type = "phase"
        modes = rand(Bool,N_modes,B.N,B.N)
    
    end

    #initial definitions
    begin

        #idx = 1:N_modes
        idx = randperm(N_modes)
        #idx = sortperm(freq[:])   #sort according to frequency

        if type == "complex"
            C = zeros(Complex{Float64}, (B.N,B.N)) #init. correction field
            holo = zeros(Complex{Float64},P,B.N,B.N)
        elseif type == "phase"
            C = ones(Complex{Float64}, (B.N,B.N)) #init. correction field
            holo = ones(Complex{Float64},P,B.N,B.N)
        end
        Corr = copy(C)
        C_iter = zeros(Complex{Float64}, (iter, B.N, B.N))
        ϕ_R = range(0, step = 2π/P, length=P)
        a = zeros(Complex, N_modes) #complex mode amplitudes
        a_old = copy(a)
       
        #signals = zeros(iter,N_modes)
        sig = zeros(N_modes)
        scat_fit = zeros(B.N, B.N)
        signals = Any[]
        gain = Any[] #enhancement factor 
                     
        PMT = fill(0.0,P)
        scat_model = zeros(B.N,B.N)

    end
    
    #running correction iterations
    for i = 1:iter

        println(i)
       
        m = 0
        while m < N_modes
            m += 1
                                     
            a_old[idx[m]] = a[idx[m]] #store old value of a

            if type == "complex"
                refbeam = √(1-f) * exp.(1im*angle.(C))
            end

            for p = 1:P #phase stepping loop

                if type == "complex"
                    holo[p,:,:] .= √f *  exp.(1im *(modes[idx[m],:,:] .+ ϕ_R[p])) .+ refbeam

                    if (method == "F-SHARP") || (method == "IMPACT")  #"complex" field updating
                        PMT[p] = TPEF(B, scat, holo[p,:,:], t_meas, bg = bg)[1]
                    elseif method == "DASH" #phase-only field shaping
                        PMT[p] = TPEF(B, scat, exp.(1im.*angle.(holo[p,:,:])), t_meas, bg = bg)[1]
                    end

                elseif type == "phase"
                    PMT[p] = TPEF(B, scat, C .* exp.(1im * modes[idx[m],:,:] .* ϕ_R[p]), t_meas, bg = bg)[1]
                end

            end
                                
            a[idx[m]] = get_phase(PMT, ϕ_R, method = "simple") #retrieving amplitude and phase of individual mode
            sig[m], I2ph = TPEF(B, scat, exp.(1im*angle.(C)), t_meas, noise = false) #current signal level when correction is applied

            if type == "complex"

                if (method == "F-SHARP") || (method == "IMPACT") 
                    wa = 1      
                    w_old = 1            
                elseif method == "DASH"
                    wa =  1 
                    w_old = 0
                end 
                
                Corr .+= conj(a[idx[m]] - w_old*a_old[idx[m]])/wa .* exp.(1im * modes[idx[m],:,:])                
            
            elseif type == "phase"            
                Corr .*= exp.(-1im*angle(a[idx[m]]) .* modes[idx[m],:,:])
            end

            #update of correction pattern
            if (method == "DASH" || type == "phase") && ((m%update == 0) || (m == N_modes))
                C .= Corr #update correcting reference wave
            end

            if (method == "F-SHARP" || method == "IMPACT") && (m == N_modes)
                C .= Corr  
            end   
        end
        
        C_iter[i,:,:] .= C #storing the correction pattern after each iteration
        
        append!(signals, sig)
        append!(gain, TPEF(B, scat, exp.(1im*angle.(C)), 1.0)[1] / TPEF(B, scat, 0im .+ ones(B.N,B.N), 1.0)[1]  )

    end

    C_iter, N_modes, signals, modes, gain
end

"""genetic algorithm"""
function GA(B, scat::Array{Float64,2}, t_meas::Float64; gen::Int=1000, N_modes::Int = 50, R = [0.1 0.01], λ = 200, pow::Int = 3, bg::Int=0)
    
    modes = 2π * rand(N_modes,B.N,B.N)
    PMT = [TPEF(B, scat, exp.(1im * modes[m,:,:]), t_meas, bg = bg)[1] for m = 1:N_modes] #evaluating modes
    idx = sortperm(PMT, rev = true) #sorting in descending order according to PMT signal
    
    #init.
    kids = zeros(Int(N_modes/2), B.N, B.N)
    PMT_kids = fill(0.0, Int(N_modes/2))
    sig = fill(0.0, gen)

    for g = 1:gen #generations
                  
        #weight = Weights(PMT[idx].^pow)
        weight = Weights((PMT[idx] .- minimum(PMT)).^pow)

        for m = 1 : Int(N_modes/2)
            
            #breeding
            mum, dad = sample(idx, weight, 2, replace = false) #selecting parents, likelihood according to PMT signal
            T = rand(Bool,B.N,B.N) #selection mask 
            kid = modes[mum,:,:] .* T .+ modes[dad,:,:] .* .!T 
            
            #mutation
            R0 = (R[1] - R[2]) * exp(-g/λ) + R[2] #fraction of mutating pixels
            mut_idx = rand(1:B.N^2,Int(round(R0*B.N^2))) #selection mask: which pixels are mutated?        
            kid[mut_idx] .= 2π*rand(length(mut_idx))    #mutation    

            kids[m,:,:] = kid
            PMT_kids[m] = TPEF(B, scat, exp.(1im * kid), t_meas, bg = bg)[1]
        
        end
        
        modes[idx[Int(N_modes/2)+1 : end],:,:] = kids #replace worse half of modes with kids, regardless of their PMT signal
        PMT[idx[Int(N_modes/2)+1 : end]] = PMT_kids #replace also the PMT signals
        idx = sortperm(PMT, rev = true) #sorting in descending order according to PMT signal
    
        sig[g] = TPEF(B, scat, exp.(1im * modes[idx[1],:,:]), t_meas)[1] #evaluating modes
        
    end

    return sig, modes[idx[1],:,:]
end

"""execution of the selected iterative algorithm"""
function RUN(B,method, gen, N_gensize, R0, λ, pow, iter, cyc, t_meas, bg, N, infl, f, type_specimen)
    
    # NA = 0.8 #NA of objective lens 
    # RI = 1.33 #refractive index of sample and immersion medium
    # λ0 = 800e-9 #excitation wavelength
    # use_noise = true #use poissonian noise?
    # NL = 2 #order of nonlinearity
    # ux = λ0/2/NA  #size of a simulated pixel 
    # N2 = infl * N
    # B = basics(N, infl, NA, RI, λ0, ux = ux, NL = NL, type_specimen = type_specimen, type_pupil = "square")
    S = [scatterer(B.N2) for m = 1:cyc] #create scatterers
    P = 3 #no. of phase steps
    update = 1 #mode no. after which SLM is updated; not used for F-SHARP and IMPACT


    #----------------- EXECUTE OPTIMIZATION ROUTINE --------------------------

    if method == "GA"

        R = [begin sig, Φ = GA(B, S[i], t_meas, gen = gen, N_modes = N_gensize, R = R0, λ = λ, pow = pow, bg = bg); (sig, Φ); end for i = 1:cyc] #execute GA
        sig, Φ = zip(R...)
        
        #plotting
        begin
            no_meas = N_gensize + gen * Int(N_gensize/2)
            println(" ")
            println("$(1e6*t_meas) µs pixel dwell time")
        
            figure(2)
            errorbar(N_gensize .+ (1:gen).*Int(N_gensize/2), mean(sig)', yerr = std(collect(sig))'/sqrt(cyc))
            title(method*", "*type_specimen*" sample, N = $N_gensize / $(B.N2^2)")
            xlabel("meas. no.")
            ylabel("photons per measurement")

        end


    else
        R = [begin 
                C, N_modes, signals, modes, gain  = SCATCORR(B, S[i], t_meas, P = P, method = method, f = f, iter = iter, update = update, bg = bg)
                (C, N_modes, signals, modes, gain) 
            end
            for i = 1:cyc]
            C, _, sig, modes, gain = zip(R...)

        #plotting
        begin
            N_modes = R[1][2]
            modes = R[1][5]
            no_meas = length(sig[1])
            maxgain = 0.5 * N_modes^2/B.N2^2 #see considerations on OneNote  (the factor 0.5 is my own "guess" supported by simulations 
            max_signal = B.η^2*maximum(B.specimen)*t_meas #expectancy value of signal for full correction (all scattered modes corrected)
            println(" ")
            println("$(1e6*t_meas) µs pixel dwell time")
            println("$N_modes out of $(B.N2^2) modes corrected")
            println("$no_meas measurements in total")
            println("$(round(sum(sig[1]), digits = 3)) photons collected ")
            println("$bg photons / sec. background")
            println("$(round(minimum(sig[1]))) photons = min signal (mean over all phase steps)")
            println("$(round(maximum(sig[1]))) photons = max signal")
            println("")

            figure(2)
            errorbar(1:P:P*no_meas, mean(sig)', yerr = std(collect(sig))'/sqrt(cyc))
            title(method*", sample = "*type_specimen*", N = $N_modes / $(B.N2^2), bg = $bg /s")
            xlabel("meas. no.")
            ylabel("photons per measurement")
            
            Φ = [angle.(C[1])[end,:,:]]
        end

    end

    I_corr = sqrt.(TPEF(B, S[1], exp.(1im*Φ[1]), t_meas, noise = false)[2])
    I_scat = sqrt.(TPEF(B, S[1], 0im .+ ones(B.N,B.N), t_meas, noise = false)[2])

    figure(1)
    clf()
    imshow(ifftshift((I_corr)), cmap = "gray")
    title("corrected irradiance")
    colorbar(); 

    figure(5)
    clf()
    imshow(ifftshift(I_scat), cmap="gray")
    title("scattered irradiance")
    colorbar() 


end
