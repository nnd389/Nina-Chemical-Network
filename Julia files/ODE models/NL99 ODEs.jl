@time begin

using OrdinaryDiffEq
using Plots
using Sundials
using CUDA


using Catalyst
using Symbolics
using DiffEqDevTools
using Sundials, ODEInterface, ODEInterfaceDiffEq, LSODA

#=
This code implements the Nelson and Langer (1999) network as a system of ODEs
I was meant to be compared to DESPOTIC's model for their predefined NL99 network, 
so the parameters should be closely defined to that of DESPOTIC's
=#

    T = 10
    Av = 2 # This is what depotic uses for visual extinction
    Go = 1.7 # I think despotic uses 1.7
    n_H = 611
    shield = 1

    function Nelson_ODE(du,u,p,t)
        # 1: H2
        #= du[1] = -1.2e-17 * u[1] + 
                n_H * (1.9e-6 * u[2] * u[3]) / (T^0.54) - 
                n_H * 4e-16 * u[1] * u[12] - 
                n_H * 7e-15 * u[1] * u[5] + 
                n_H * 1.7e-9 * u[10] * u[2] + 
                n_H * 2e-9 * u[2] * u[6] + 
                n_H * 2e-9 * u[2] * u[14] + 
                n_H * 8e-10 * u[2] * u[8] =#
        du[1] = 0

        # 2: H3+
        du[2] = 1.2e-17 * u[1] + 
                n_H * (-1.9e-6 * u[3] * u[2]) / (T^0.54) - 
                n_H * 1.7e-9 * u[10] * u[2] - 
                n_H * 2e-9 * u[2] * u[6] - 
                n_H * 2e-9 * u[2] * u[14] - 
                n_H * 8e-10 * u[2] * u[8]

        # 3: e
        du[3] = n_H * (-1.4e-10 * u[3] * u[12]) / (T^0.61) - 
                n_H * (3.8e-10 * u[13] * u[3]) / (T^0.65) - 
                n_H * (3.3e-5 * u[11] * u[3]) / T + 
                1.2e-17 * u[1] - 
                n_H * (1.9e-6 * u[3] * u[2]) / (T^0.54) + 
                6.8e-18 * u[4] - 
                n_H * (9e-11 * u[3] * u[5]) / (T^0.64) + 
                3e-10 * Go * exp(-3 * Av) * u[6] +
                n_H * 2e-9 * u[2] * u[13] ## CHECK I added this extra term from a CR ionization reaction
                + 2.0e-10 * Go * exp(-1.9 * Av) * u[14] # this term was added as part of the skipped photoreaction
        
        
        # 4: He
        du[4] = n_H * (9e-11 * u[3] * u[5]) / (T^0.64) - 
                6.8e-18 * u[4] + 
                n_H * 7e-15 * u[1] * u[5] + 
                n_H * 1.6e-9 * u[10] * u[5]
        #du[4] = 0
        
        # 5: He+   6.8e-18 or 1.1
        du[5] = 6.8e-18 * u[4] - 
                n_H * (9e-11 * u[3] * u[5]) / (T^0.64) - 
                n_H * 7e-15 * u[1] * u[5] - 
                n_H * 1.6e-9 * u[10] * u[5]
        #u[5] = 0
        
        # 6: C
        du[6] = n_H * (1.4e-10 * u[3] * u[12]) / (T^0.61) - 
                n_H * 2e-9 * u[2] * u[6] - 
                n_H * 5.8e-12 * (T^0.5) * u[9] * u[6] + 
                1e-9 * Go * exp(-1.5 * Av) * u[7] - 
                3e-10 * Go * exp(-3 * Av) * u[6] + 
                1e-10 * Go * exp(-3 * Av) * u[10] * shield

        # 7: CHx
        du[7] = n_H * (-2e-10) * u[7] * u[8] + 
                n_H * 4e-16 * u[1] * u[12] + 
                n_H * 2e-9 * u[2] * u[6] - 
                1e-9 * Go * u[7] * exp(-1.5 * Av)
        
        # 8: O
        du[8] = n_H * (-2e-10) * u[7] * u[8] + 
                n_H * 1.6e-9 * u[10] * u[5] - 
                n_H * 8e-10 * u[2] * u[8] + 
                5e-10 * Go * exp(-1.7 * Av) * u[9] + 
                1e-10 * Go * exp(-3 * Av) * u[10] * shield
        
        # 9: OHx
        du[9] = n_H * (-1e-9) * u[9] * u[12] + 
                n_H * 8e-10 * u[2] * u[8] - 
                n_H * 5.8e-12 * (T^0.5) * u[9] * u[6] - 
                5e-10 * Go * exp(-1.7 * Av) * u[9]

        # 10: CO
        du[10] = n_H * (3.3e-5 * u[11] * u[3]) / T + 
                n_H * 2e-10 * u[7] * u[8] - 
                n_H * 1.7e-9 * u[10] * u[2] - 
                n_H * 1.6e-9 * u[10] * u[5] + 
                n_H * 5.8e-12 * (T^0.5) * u[9] * u[6] - 
                1e-10 * Go * exp(-3 * Av) * u[10] + 
                1.5e-10 * Go * exp(-2.5 * Av) * u[11] * shield
        
        # 11: HCO+
        du[11] = n_H * (-3.3e-5 * u[11] * u[3]) / T + 
                n_H * 1e-9 * u[9] * u[12] + 
                n_H * 1.7e-9 * u[10] * u[2] - 
                1.5e-10 * Go * exp(-2.5 * Av) * u[11]
        
        # 12: C+
        du[12] = n_H * (-1.4e-10 * u[3] * u[12]) / (T^0.61) - 
                n_H * 4e-16 * u[1] * u[12] - 
                n_H * 1e-9 * u[9] * u[12] + 
                n_H * 1.6e-9 * u[10] * u[5] + 
                3e-10 * Go * exp(-3 * Av) * u[6]
        
        # 13: M+
        du[13] = n_H * (-3.8e-10 * u[13] * u[3]) / (T^0.65) + 
                n_H * 2e-9 * u[2] * u[14] 
                + 2.0e-10 * Go * exp(-1.9 * Av) * u[14] # this term was added as part of the skipped photoreaction
        
        # 14: M
        du[14] = n_H * (3.8e-10 * u[13] * u[3]) / (T^0.65) - 
                n_H * 2e-9 * u[2] * u[14] 
                - 2.0e-10 * Go * exp(-1.9 * Av) * u[14] # this term was added as part of the skipped photoreaction

    end

    tspan = (0.0, 30 * 3.16e10) # ~30 thousand yrs
    u0 = [0.5 ;    # 1:  H2   yep?
        9.059e-9 ; # 2:  H3+  yep
        2.0e-4 ;   # 3:  e    yep
        0.1 ;              # 4:  He  SEE lines 535 NL99
        7.866e-7 ; # 5:  He+  yep? should be 2.622e-5
        0.0 ;      # 6:  C    yep
        0.0 ;      # 7:  CHx  yep
        0.0004 ;   # 8:  O    yep
        0.0 ;      # 9:  OHx  yep
        0.0 ;      # 10: CO   yep
        0.0 ;      # 11: HCO+ yep
        0.0002 ;   # 12: C+   yep
        2.0e-7 ;   # 13: M+   yep
        2.0e-7 ]   # 14: M    yep


    

    prob = ODEProblem(Nelson_ODE, u0, tspan)
    #Autotsit5 rosenbrock
    #sol = solve(prob, Euler(); dt = 1e6)
    #sol = solve(prob, CVODE_Adams())
    sol = solve(prob, Rodas4P()) 
    #sol = solve(prob, AutoTsit5(Rosenbrock23())) #incorrect HCO+

    #sol1 = solve(prob, Euler(); dt = 9e6)
    #dt = 1e6
    #sol2 = solve(prob, Rodas4P())   #does a great job!!! and is way faster than Euler()
    #sol3 = solve(prob, Rodas4())    #also does a good job
    #sol4 = solve(prob, Rodas5P())    #also does a good job

    #sol5 = solve(prob, FBDF())      #incorrect HCO+
    #sol6 = solve(prob, QNDF())      #incorrect HCO+
    #sol7 = solve(prob, CVODE_BDF()) #incorrect HCO+
    #sol8 = solve(prob, lsoda())     #incorrect HCO+
    



    plot(sol, vars = (0,11), label = "Reference Solution", linewidth = 3, lc=:orange, xlabel = "time (s)", ylabel = "abundance per H nucleus", title = "Abundance of HCO+")
    #plot(sol2, vars = (0,11), label = "Rodas4P", linewidth = 3, lc=:orange)
    #plot!(sol3, vars = (0,11), label = "Rodas4", linewidth = 3, lc=:yellow)
    #plot!(sol4, vars = (0,11), label = "Rodas5P", linewidth = 3, lc=:pink)
    #plot!(sol7, vars = (0,11), label = "CVODE_BDF", linewidth = 3, lc=:black)
    #plot!(sol6, vars = (0,11), label = "QNDF", linewidth = 3, lc=:blue)
    #plot!(sol5, vars = (0,11), label = "FBDF", linewidth = 3, lc=:green)
    #plot!(sol8, vars = (0,11), label = "lsoda", linewidth = 3, lc=:purple)

   

    #palette=:PRGn



# timetepping! make it smaller
# look into how the timestep is chosen
# try to see how fast julia runs
# figure out timestepping for both codes


# Lets start forming groups!
# C and C+
# groups were a bust! we're gonna have to do it one by one 

    
    #plot(sol, vars = (0,6), label = "C", linewidth = 3, lc=:mediumpurple3)
    #plot!(sol, vars = (0,12), label = "C+", linewidth = 3, lc=:forestgreen, xlabel = "time (s)", ylabel = "abundance per H nucleus", title = "Abundance of C and C+")
    

    

    #print(length(sol))
    #plot(sol, vars = (0,8), label = "O", linewidth = 3, lc=:green, xlabel = "time (s)", ylabel = "abundance per H nucleus", title = "Abundance of O")
    #plot(sol, vars = (0,5), label = "He+", linewidth = 3, lc=:pink, xlabel = "time (s)", ylabel = "abundance per H nucleus", title = "Abundance of He+")
    #plot(sol, vars = (0,13), label = "M+", linewidth = 3, lc=:mediumpurple3, xlabel = "time (s)", ylabel = "abundance per H nucleus", title = "Abundance of M+")
    #plot(sol, vars = (0,2), label = "H3+", linewidth = 3, lc=:orange, ylimits = (0, 2.75e-10), xlabel = "time (s)", ylabel = "abundance per H nucleus", title = "Abundance of H3+") 
    

    #xtick_labels = ["0", "6.01e3", "1.20e4", "1.80e4", "2.41e4"]
    #xticks_positions = collect(0:2e11:8e11)
    #gren = palette(:Greens)
#=
    OHX = plot(sol, vars = (0,9), label = "OHx",  dpi = 600, legend = false, titlefontsize = 12, xguidefontsize=9, yguidefontsize=9, tickfontsize=7, linewidth = 4, lc=gren[3],         xlabel = "Time (yrs)", title = "OHx")
    CO  = plot(sol, vars = (0,10), label = "CO",  dpi = 600, legend = false, titlefontsize = 12, xguidefontsize=9, yguidefontsize=9, tickfontsize=7, linewidth = 4, lc=gren[5],         xlabel = "", ylabel = "Abundance per H Nucleus", title = "CO")
    CHX = plot(sol, vars = (0,7),  label = "CHx", dpi = 600, legend = false, titlefontsize = 12, xguidefontsize=9, yguidefontsize=9, tickfontsize=7, linewidth = 4, lc=gren[7], xlabel = "Time (yrs)", ylabel = "Abundance per H Nucleus", title = "CHx") 
    M   = plot(sol, vars = (0,14), label = "M",   dpi = 600, legend = false, titlefontsize = 12, xguidefontsize=9, yguidefontsize=9, tickfontsize=7, linewidth = 4, lc=gren[9], xlabel = "", title = "M")
    
    abun = plot(CO, M, CHX, OHX, layout = (2,2), xticks=(xticks_positions, xtick_labels))
    display(abun)
    savefig(abun, "abun.svg")
=#
#=
    # These require the euler solver

    colors = palette(:acton, 5)

    p1 = plot(sol, vars = (0,11),  lc=colors[1], dpi = 600, yticks=0:(1.2e-13)/4:1.2e-13, legend = false,  titlefontsize = 12, xguidefontsize=8, yguidefontsize=8, tickfontsize=6, linewidth = 3, xlabel = "", title = "Reference solution")
    p2 = plot(sol5, vars = (0,11), lc=colors[2], dpi = 600, legend = false, titlefontsize = 12, xguidefontsize=8, yguidefontsize=8, tickfontsize=6, linewidth = 3, xlabel = "", title = "FBDF")
    p3 = plot(sol7, vars = (0,11), lc=colors[3], dpi = 600, legend = false, titlefontsize = 12, xguidefontsize=8, yguidefontsize=8, tickfontsize=6, linewidth = 3, xlabel = "", title = "CVODE_BDF")
    p4 = plot(sol8, vars = (0,11), lc=colors[4], dpi = 600, legend = false, titlefontsize = 12, xguidefontsize=8, yguidefontsize=8, tickfontsize=6, linewidth = 3, xlabel = "", title = "LSODA")
    
    #combined_plot = plot(p1, p2, p3, p4, layout=(4, 1), dpi = 600, pallete=:acton, xticks=(xticks_positions, xtick_labels))
    combined_plot = plot(p1, p2, p3, p4, layout=(4, 1), dpi = 600, pallete=:acton)

    # Display the combined plot
    display(combined_plot)
    savefig(combined_plot, "my_plot_600.svg")

=#
    # These are not in DESPOTIC
    #plot!(sol, vars = (0,3), label = "e", linewidth = 3, lc=:green, xlabel = "time (s)", ylabel = "abundance per H nucleus", title = "Abundance of He, He+, and e")
    #plot(sol, vars = (0,4), label = "He", linewidth = 3,  lc=:orange)
    #plot(sol, vars = (0,1), label = "H2", linewidth = 3, lc=:blue1) 
    #plot(sol, vars = (0,14), label = "M", legend = false, titlefontsize = 26, xguidefontsize=16, yguidefontsize=16, tickfontsize=11, linewidth = 3, lc=:mediumpurple3, xlabel = "time(s)", ylabel = "abundance per H nucleus", title = "Abundance of M")

#=
    # Run Benchmark

abstols = 1.0 ./ 10.0 .^ (7:13)
reltols = 1.0 ./ 10.0 .^ (4:10)

setups = [
          Dict(:alg=>Rodas4P()),
          Dict(:alg=>Rodas4()),
          Dict(:alg=>Rodas5P()),

          Dict(:alg=>FBDF()),
          Dict(:alg=>QNDF()),
          Dict(:alg=>CVODE_BDF()),
          Dict(:alg=>lsoda()),
          #Dict(:alg=>AutoTsit5(Rosenbrock23())),
          #Dict(:alg=>Euler(), :dt=>1e6)
          ]

bigprob = remake(prob, u0 = big.(prob.u0), tspan = big.(prob.tspan))
refsol = solve(bigprob, Euler(); dt = 1e6, abstol=1e-18, reltol=1e-18)

wp = WorkPrecisionSet(prob,abstols,reltols,setups;verbose=false,
                      save_everystep=false,appxsol=refsol,maxiters=Int(1e5),numruns=10
                      )
work_precision = plot(wp; palette=:PRGn_8, 
                xlabel = "Error", ylabel = "Time to solve (s)",
                xguidefontsize=16, yguidefontsize=16, tickfontsize=11)

display(work_precision)
savefig(work_precision, "work_pres.svg")

#palette=:acton10

=#



end
