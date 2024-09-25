using Catalyst
two_state_model = @reaction_network begin
    (k1,k2), X1 <--> X2
end
diffusion_rx = @transport_reaction D X1
lattice = CartesianGrid((5,5))
lrs = LatticeReactionSystem(two_state_model, [diffusion_rx], lattice)

u0 = [:X1 => rand(5, 5), :X2 => 2.0]
tspan = (0.0, 10.0)
ps = [:k1 => 1.0, :k2 => 2.0, :D => 0.2]
oprob = ODEProblem(lrs, u0, tspan, ps)
using OrdinaryDiffEq
sol = solve(oprob)
import CairoMakie
lattice_animation(sol, :X1, lrs, "lattice_simulation_2d.mp4")