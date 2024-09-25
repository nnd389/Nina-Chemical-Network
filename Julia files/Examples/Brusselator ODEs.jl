using Catalyst
using OrdinaryDiffEq
using Plots
using Latexify


function Bruss(du,u,p,t)
    du[1] = 1.0*1.0 - 
            1.0*3.0*u[1] + 
            1.0*(u[1]^2)*u[2] - 
            1.0*u[1]

    du[2] = 1.0*3.0*u[1] - 
            1.0*(u[1]^2)*u[2]
end

function B_jac(J,u,p,t)
    J[1,1] = -1.0 * B(t) + 2 * 1.0 * X(t) * Y(t) - 1.0
    J[1,2] = 1.0 * X(t)^2
    J[2,1] = 1.0 * B(t) - 2.0 * 1.0 * X(t) * Y(t)
    J[2,2] = -1.0 * X(t)^2
    nothing
end

u0 = [1.0 ; 1.0]
tspan = (0.0, 20.0)
#params = [:A => 1, :B => 3]

ff = ODEFunction(Bruss; jac=B_jac)
prob = ODEProblem(ff, u0, tspan)
sol = solve(prob, Tsit5())
plot(sol)



