#Pkg.add("Graphviz")
#using Graphviz
using Catalyst
using OrdinaryDiffEq
using Plots
using Latexify
# Despotic uses updated rates, I am just planning to use what despotic uses
Go = 1.7
Av = 2
shield_CO = 1
shield_H2 = 1


Nelson = @reaction_network begin

    # Cosmic-ray Ionization
    6.0e-18, H --> H⁺ + el
    1.2e-17, H2 --> H3⁺ + H + el 
    6.8e-18, He --> He⁺ + el

    # Photoreactions
    3.3e-11 * shield_H2 * Go * exp(-3.74 * Av), H2 --> H + H  # added reaction
    1.2e-10 * shield_CO * Go * exp(-3.53 * Av), CO --> C + O #updated rate: 1e-10 -> 1.2e-10
    1.8e-10 * Go * exp(-3 * Av), C --> C⁺ + el # updated rate: 3e-10 -> 1.8e-10
    1e-9 * Go * exp(-1.5 * Av), CHx --> C + H 
    5e-10 * Go * exp(-1.7 * Av), OHx --> O + H 
    2e-10 * Go * exp(-1.9 * Av), M --> M⁺ + el 
    1.5e-10 * Go * exp(-2.5 * Av), HCO⁺ --> CO + H⁺  # is nelson wrong? Nelson days its HCO+ --> CO + H


    # Two Body Reactions
    2e-9, H3⁺ + C -> CHx + H2 
    8e-10, H3⁺  + O   -> OHx + H2
    1.7e-9, H3⁺  + CO  -> HCO⁺ + H2
    7e-15, He⁺  + H2  -> He + H + H⁺
    1.4e-9*((T/300)^(-0.5)), He⁺  + CO  -> C⁺ + O + He  #Update rate: 1.6e-9 -> 1.4e-9*(T/300.)**-0.5
    4e-16, C⁺   + H2  -> CHx + H
    1e-9, C⁺ + OHx -> HCO⁺
    2e-10, O + CHx -> CO + H
    5.8e-12 * (T)^(0.5), C + OHx -> CO + H
    (1.0e-11/sqrt(T)) * (11.19 - 1.676*log(T) - 0.2852*log(T)^2 + 0.04433*log(T)^3), He⁺ + el  -> He #Updated rate: 9e-11 * (T)^(-0.64) -> 1.0e-11/np.sqrt(T) * \ (11.19 - 1.676*logT - 0.2852*logT**2 + 0.04433*logT**3)
    2.34e-8 * (T/300)^(-0.52), H3⁺ + el  -> H2 + H # Updated rate: 1.9e-6 * (T)^(-.54) ->2.34e-8 * (T/300.)**-0.52
    4.36e-8 * (T/300)^(-0.52), H3⁺ + el  -> 3H # new reaction
    4.67e-12 * (T/300)^(-0.6), C⁺ + el  -> C # Updated rate: 1.4e-10 * (T)^(-0.61) -> 4.67e-12 * (T/300.)**-0.6
    2.76e-7 * (T/300)^(-0.64), HCO⁺ + el  -> CO + H # updated rate: 3.3e-5 * (T)^(-1.0) -> 2.76e-7 * (T/300.)**-0.64
    2.753e-14*(315614/T)^(1.5) * (1.0+(115188/T)^(0.407))^(-2.242), H⁺ + el  -> H # New reaction
    #15 below
    1, H2 + H   -> 3H # New reaction
    1, H2 + H2  -> H2 + 2H # New reaction
    1, H + el  -> H⁺ + 2el # New reaction
    1, He⁺ + H2  -> H⁺ + He + H # New reaction
    1, H + H   -> H2 # Be careful with grain = true
    1, H⁺ + el  -> H # Be careful with grain = true    
    # finish here
    
    3.8e-10 * (T)^(-0.65), M⁺ + el  -> M
    2e-9, H3⁺ + M   -> M⁺ + H + H2
end


n_h2 = 100

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

tspan = (0.001, 10)
params = [:T => 10.0, :G₀ => 1.0, :Aᵥ => 10.0, :CO12 => 1300]
oprob = ODEProblem(Nelson, u0, tspan, params; combinatoric_ratelaws = false)
sol = solve(oprob, Tsit5())
Plots.plot(sol)

#latexify(Nelson; form = :ode)
