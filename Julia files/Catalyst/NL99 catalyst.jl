#Pkg.add("Graphviz")
#using Graphviz
using Catalyst
using OrdinaryDiffEq
using Plots
using Latexify

# This is a catalyst program to help recreate the DESPOTIC model of Neslon 1999 reaction network 
# NOTE: I only used this to help create the ODE's for me, see "Nelson ODEs" for the actual model, 
# I don't think this runs correctly, since I only needed it for the latexify function

# Observations:
# - DESPOTIC leaves H and H+ out lof the network
# - DESPOTIC treats O and OI the same as well as C and CI
# - DESPOTIC skips the photoreaction M -> M+  +  e


Nelson = @reaction_network begin
    # Cosmic-ray Ionization
    1.2e-17, H2 --> H3⁺ + e + H
    6.8e-18, He --> He⁺ + e

    # Ion-Molecule Reactions
    2e-9, H3⁺ + CI --> CHₓ + H2
    8e-10, H3⁺ + OI --> OHₓ + H2 
    1.7e-9, H3⁺ + CO --> HCO⁺ + H2 
    7e-15, He⁺ + H2 --> He + H + H⁺
    1.6e-9, He⁺ + CO --> C⁺ + O + He
    4e-16, C⁺ + H2 --> CHₓ + H
    1e-9, C⁺ + OHₓ --> HCO⁺

    # Neutral-Neutral Reactions
    # T = 60K? (See section 4.1)
    2e-10, OI + CHₓ --> CO + H
    5.8e-12 * (T)^(0.5), CI + OHₓ --> CO + H 

    # Electron recombination
    # T = 60K? (See section 4.1)
    9e-11 * (T)^(-0.64), He⁺ + e --> He
    1.9e-6 * (T)^(-.54), H3⁺ + e --> H + H2 
    1.4e-10 * (T)^(-0.61), C⁺ + e --> CI
    3.3e-5 * (T)^(-1.0), HCO⁺ + e --> CO + H 
    3.8e-10 * (T)^(-0.65), M⁺ + e --> M

    # Charge-Transfer Reactions
    2e-9, H3⁺ + M --> M⁺ + e + H2 

    # Photoreactions FIX G₀ and FIX Aᵥ
    # G₀ = 1 (See section 4.1)
    # Aᵥ is the visual extinction, see glover pg3 for more?
    3e-10 * G₀ * exp(-3 * Aᵥ), CI --> C⁺ + e 
    1e-9 * G₀ * exp(-1.5 * Aᵥ), CHₓ --> CI + H 
    10e-10 * (CO12) * G₀ * exp(-3 *  Aᵥ), CO --> CI + O # FIX 
    5e-10 * G₀ * exp(-1.7 * Aᵥ), OHₓ --> OI + H 
    # here lies the skipped photoreaction:
    2e-10 * G₀ * exp(-1.9 * Aᵥ), M --> M⁺ + e 
    1.5e-10 * G₀ * exp(-2.5 * Aᵥ), HCO⁺ --> CO + H 

    # I have deleted all of the cr's and hv's because the ODE's will turn out the same without them
    # ADD CONSERVATION CONDITIONS
end


# n_H is 611 in Despotic
n_h2 = 611

u0 = [
    :He⁺ => 1, 
    :H3⁺ => 1, 
    :OHₓ => 1, 
    :CHₓ => 0, 
    :CO => 0, 
    :CI => 0, 
    :C⁺ => 10e-4 * (n_h2), 
    :HCO⁺ => 0, 
    :OI => 2e-4 * (n_h2), 
    :M => 10e-7 * (n_h2), 
    :M⁺ => 10e-7 * (n_h2), 
    :e => 10^-4 * (n_h2) + 1 + 10e-7,

    :H2 => 50, 
    :He => 0.28 * (n_h2), 
    :O => 2e-4 * (n_h2), 
    :H => 1, 
    :H⁺ => 1
    ]
# T_dust = 10-30k starts looking different above 22k Bergin et al 1995
# T_d = 10k see nelson 1990 section 2.2.1
# for number density, n_h >~ 10e7 cm^-3 see nelson pg 2
# T = 60       # in Kelvin, (See section 4.1)
# G₀ = 1       # (See section 4.1)
# Aᵥ = 10      # Incorrect, Aᵥ is the visual extinction, see glover pg3 for more?
# CO12 = 880-1600 Aᵒ # See Discheock and Black pg 4

tspan = (0.001, 10)
params = [:T => 10.0, :G₀ => 1.0, :Aᵥ => 10.0, :CO12 => 1300]
oprob = ODEProblem(Nelson, u0, tspan, params; combinatoric_ratelaws = false)
sol = solve(oprob, Tsit5())
Plots.plot(sol)

#Plots.plot(sol, vars = (0,8), xscale=:log10)
#Plots.plot(sol, xscale =:log10; idxs = [:C⁺, :CO])
#Plots.plot(sol, vars = (0,11), xscale=:log10, yscale=:log10)


#latexify(Nelson; form = :ode)


#= SPECIES LABELS
1: H2
2: H3+
3: e
4: H
5: He
6: He+
7: CI
8: CHₓ
9: OI
10: OHₓ
11: CO
12: HCO+
13: H+
14: C+
15: O
16: M+
17: M
=#