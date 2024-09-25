using OrdinaryDiffEq, CUDA, LinearAlgebra
u0 = cu(rand(1000))
A = cu(randn(1000, 1000))
f(du, u, p, t) = mul!(du, A, u)
prob = ODEProblem(f, u0, (0.0f0, 1.0f0)) # Float32 is better on GPUs!
sol = solve(prob, Tsit5())