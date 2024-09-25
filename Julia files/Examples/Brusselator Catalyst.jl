# %%
using Catalyst
using OrdinaryDiffEq
using Plots
using Latexify

Brusselator = @reaction_network begin
    #@species X(t)=1.0 Y(t)=1.0
    # @parameters A [isconstantspecies=true] B [isconstantspecies=true] E [isconstantspecies=true] D [isconstantspecies=true]
    1.0, A --> X
    1.0, B + X --> Y + D 
    1.0, 2*X + Y --> 3*X 
    1.0, X --> E
end

function B2_jac(J,u,p,t)
    J[1,1] = -1.0 * B(t) + 2 * 1.0 * X(t) * Y(t) - 1.0
    J[1,2] = 1.0 * X(t)^2
    J[2,1] = 1.0 * B(t) - 2.0 * 1.0 * X(t) * Y(t)
    J[2,2] = -1.0 * X(t)^2
    nothing
end

# %%
u0 = [:X => 1.0, :Y => 1.0]
tspan = (0.0, 20.0)
params = [:A => 1.0, :B => 3.0, :E => 1.0, :D => 1.0]
oprob = ODEProblem(Brusselator, u0, tspan, params; jac = true, combinatoric_ratelaws = false)
sol = solve(oprob, Tsit5())
Plots.plot(sol)



#osys = convert(ODESystem, Brusselator; combinatoric_ratelaws = false)
#Brusselator2 = ODEFunction(osys, jac=B2_jac)
# Not sure how to define custom jacobian, I want to
# do it the same way as before, but my brusselator function
# is not technically ad ODE function, its a list of Reactions




