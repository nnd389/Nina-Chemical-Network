using Catalyst
using OrdinaryDiffEq
using Plots

# I deleted all the gamma, cr, and gamma_cr, and H(s)

Glover = @reaction_network begin
    8.4668e-16, H + e --> H_m
    1.5e-09, H_m + H --> H2 + e
    4.8607e-18, H_p + H --> H2_p
    6.4e-10, H + H2_p --> H2 + H_p
    7.9689e-08, H_m + H_p --> H + H
    1e-08, H2_p + e --> H + H
    15852692190.9159, H2 + H_p --> H2_p + H
    5.3224e-52, H2 + e --> H + H + e
    1.8752e-38, H2 + H --> H + H + H
    3.873547154666856e+19, H2 + H2 --> H2 + H + H
    2.1277e-08, H + e --> H_p + e + e
    2.0211e-12, H_p + e --> H
    4.8183e-07, H_m + e --> H + e + e
    0.00056806, H_m + H --> H + H + e
    6.1496e-10, H_m + H_p --> H2_p + e
    2.1306e-08, He + e --> He_p + e + e
    2.0113e-11, He_p + e --> He
    1.689e-15, He_p + H --> He + H_p
    3.0049e-67, He + H_p --> He_p + H
    2.2677e-12, C_p + e --> C
    1.5629e-12, O_p + e --> O
    1.0802e-07, C + e --> C_p + e + e
    9.4858e-08, O + e --> O_p + e + e
    8.5052e-10, O_p + H --> O + H_p
    6.3993e-10, O + H_p --> O_p + H
    6.6617e-15, O + He_p --> O_p + He
    1.6985e-15, C + H_p --> C_p + H
    9.8592e-90, C_p + H --> C + H_p
    1.6014e-14, C + He_p --> C_p + He
    7.709e-46, H2 + He --> H + He
    4.705e-31, OH + H --> O + H + H
    3.8e-10, HOC_p + H2 --> HCO_p + H2
    4e-10, HOC_p + CO --> HCO_p + CO
    5.5071e-15, C + H2 --> CH + H
    1.2093e-10, CH + H --> C + H2
    7.8227e-11, CH + H2 --> CH2 + H
    6.59e-11, CH + C --> C2 + H
    6.6e-11, CH + O --> CO + H
    6.64e-11, CH2 + H --> CH + H2
    1.33e-10, CH2 + O --> CO + H + H
    8e-11, CH2 + O --> CO + H2
    9.1287e-11, C2 + O --> CO + C
    3.4728e-13, O + H2 --> OH + H
    2.8951e-13, OH + H --> O + H2
    2.2521e-12, OH + H2 --> H2O + H
    1e-10, OH + C --> CO + H
    3.5e-11, OH + O --> O2 + H
    6.1923e-12, OH + OH --> H2O + H
    4.5215e-15, H2O + H --> H2 + OH
    7.4909e-14, O2 + H --> OH + O
    9.8397e-20, O2 + H2 --> OH + OH
    3.1212e-11, O2 + C --> CO + O
    6.6007e-44, CO + H --> C + OH
    2.3062e-09, H2_p + H2 --> H3_p + H
    1.8209e-16, H3_p + H --> H2_p + H2
    2.4e-09, C + H2_p --> CH_p + H
    2e-09, C + H3_p --> CH_p + H2
    9.6577e-13, C_p + H2 --> CH_p + H
    7.5e-10, CH_p + H --> C_p + H2
    1.2e-09, CH_p + H2 --> CH2_p + H
    3.5e-10, CH_p + O --> CO_p + H
    1.4e-09, CH2 + H_p --> CH_p + H2
    8.4177e-13, CH2_p + H --> CH_p + H2
    1.6e-09, CH2_p + H --> CH_p + H2
    7.5e-10, CH2_p + O --> HCO_p + H
    1.8153e-14, CH3_p + H --> CH2_p + H2
    4e-10, CH3_p + O --> HCO_p + H2
    4.8e-10, C2 + O_p --> CO_p + C
    1.7e-09, O_p + H2 --> OH_p + H
    1.5e-09, O + H2_p --> OH_p + H
    8.4e-10, O + H3_p --> OH_p + H2
    1.3e-09, OH + H3_p --> H2O_p + H2
    7.7e-10, OH + C_p --> CO_p + H
    1.01e-09, OH_p + H2 --> H2O_p + H
    6.4e-10, H2O_p + H2 --> H3O_p + H
    5.9e-09, H2O + H3_p --> H3O_p + H2
    9e-10, H2O + C_p --> HCO_p + H
    1.8e-09, H2O + C_p --> HOC_p + H
    1e-11, H3O_p + C --> HCO_p + H2
    3.8e-10, O2 + C_p --> CO_p + O
    6.2e-10, O2 + C_p --> CO + O_p
    9.1e-10, O2 + CH2_p --> HCO_p + OH
    5.2e-11, O2_p + C --> CO_p + O
    2.7e-11, CO + H3_p --> HOC_p + H2
    1.7e-09, CO + H3_p --> HCO_p + H2
    1.1e-09, HCO_p + C --> CO + CH_p
    2.5e-09, HCO_p + H2O --> CO + H3O_p
    7.2e-15, H2 + He_p --> He + H2_p
    3.5727e-14, H2 + He_p --> He + H + H
    1.9e-09, CH + H_p --> CH_p + H
    1.4e-09, CH2 + H_p --> CH2_p + H
    7.5e-10, CH2 + He_p --> C_p + He + H2
    1.6e-09, C2 + He_p --> C_p + He + H2
    2.1e-09, OH + H_p --> OH_p + H
    1.1e-09, OH + He_p --> O_p + He + H
    6.9e-09, H2O + H_p --> H2O_p + H
    2.04e-10, H2O + He_p --> OH + He + H_p
    2.86e-10, H2O + He_p --> OH_p + He + H
    6.05e-11, H2O + He_p --> H2O_p + He
    2e-09, O2 + H_p --> O2_p + H
    3.3e-11, O2 + He_p --> O2_p + He
    1.1e-09, O2 + He_p --> O_p + O + He
    5.2e-11, O2_p + C --> O2 + C_p
    9995.0705, CO + He_p --> C_p + O + He
    7.6681e-17, CO + He_p --> C + O_p + He
    7.5e-10, CO_p + H --> CO + H_p
    1.2598e-07, C_m + H_p --> C + H
    1.2598e-07, O_m + H_p --> O + H
    1.2971e-07, He_p + H_m --> He + H
    1.2512e-08, H3_p + e --> H2 + H
    2.3313e-08, H3_p + e --> H + H + H
    3.8341e-08, CH_p + e --> C + H
    7.7695e-08, CH2_p + e --> CH + H
    1.9569e-07, CH2_p + e --> C + H + H
    3.7294e-08, CH2_p + e --> C + H2
    4.2448e-08, CH3_p + e --> CH2 + H
    1.0681e-07, CH3_p + e --> CH + H2
    1.2356e-07, CH3_p + e --> CH + H + H
    3.5348e-09, OH_p + e --> O + H
    1.6706e-07, H2O_p + e --> O + H + H
    2.1361e-08, H2O_p + e --> O + H2
    4.7104e-08, H2O_p + e --> OH + H
    5.9154e-08, H3O_p + e --> H + H2O
    3.2973e-08, H3O_p + e --> OH + H2
    1.4131e-07, H3O_p + e --> OH + H + H
    3.0672e-09, H3O_p + e --> O + H + H2
    8.395e-08, O2_p + e --> O + O
    1.4182e-07, CO_p + e --> C + O
    1.2772e-07, HCO_p + e --> CO + H
    1.1106e-08, HCO_p + e --> OH + C
    3.3e-08, HOC_p + e --> CO + H
    1e-09, H_m + C --> CH + e
    1e-09, H_m + O --> OH + e
    1e-10, H_m + OH --> H2O + e
    5e-10, C_m + H --> CH + e
    1e-13, C_m + H2 --> CH2 + e
    5e-10, C_m + O --> CO + e
    5e-10, O_m + H --> OH + e
    7e-10, O_m + H2 --> H2O + e
    5e-10, O_m + C --> CO + e
    1e-16, H2 + H_p --> H3_p
    2.25e-15, C + e --> C_m
    1e-17, C + H --> CH
    1e-17, C + H2 --> CH2
    5.65512683915318e+18, C + C --> C2
    2.1e+19, C + O --> CO
    1.3425e-17, C_p + H --> CH_p
    3.144e-16, C_p + H2 --> CH2_p
    2.5e-18, C_p + O --> CO_p
    1.5e-15, O + e --> O_m
    6.265291349000606e+19, O + H --> OH
    3.283556530235987e+21, O + O --> O2
    8963387880202277, OH + H --> H2O
    8.353721798667474e+31, H + H + H --> H2 + H
    4.437700938891119e+29, H + H + H2 --> H2 + H2
    4.353605676913333e+31, H + H + He --> H2 + He
    7.866450774189283e+34, C + C + M --> C2 + M
    1.510477208797841e+28, C + O + M --> CO + M
    96000, C_p + O + M --> CO_p + M
    96000, C + O_p + M --> CO_p + M
    1.299e+32, O + H + M --> OH + M
    2.303999999999999e+30, OH + H + M --> H2O + M
    2.76e+34, O + O + M --> O2 + M
    3.397e-11, O + CH --> HCO_p + e
    4.07118e-13, H --> H2
    7.1e-07, H_m --> H + e
    1.1e-09, H2_p --> H + H_p
    5.6e-11, H2 --> H + H
    4.9e-13, H3_p --> H2 + H_p
    4.9e-13, H3_p --> H2_p + H
    3.1e-10, C --> C_p + e
    2.4e-07, C_m --> C + e
    8.7e-10, CH --> C + H
    7.7e-10, CH --> CH_p + e
    2.6e-10, CH_p --> C + H_p
    7.1e-10, CH2 --> CH + H
    5.9e-10, CH2 --> CH2_p + e
    4.6e-10, CH2_p --> CH_p + H
    1e-09, CH3_p --> CH2_p + H
    1e-09, CH3_p --> CH_p + H2
    1.5e-10, C2 --> C + C
    2.4e-07, O_m --> O + e
    3.7e-10, OH --> O + H
    1.6e-12, OH --> OH_p + e
    1e-12, OH_p --> O + H_p
    6e-10, H2O --> OH + H
    3.2e-11, H2O --> H2O_p + e
    5e-11, H2O_p --> H2_p + O
    5e-11, H2O_p --> H_p + OH
    5e-11, H2O_p --> O_p + H2
    1.5e-10, H2O_p --> OH_p + H
    2.5e-11, H3O_p --> H_p + H2O
    2.5e-11, H3O_p --> H2_p + OH
    7.5e-12, H3O_p --> H2O_p + H
    2.5e-11, H3O_p --> OH_p + H2
    5.6e-11, O2 --> O2_p + e
    7e-10, O2 --> O + O
    2e-10, CO --> C + O
    1, H --> H_p + e
    1.1, He --> He_p + e
    0.037, H2 --> H_p + H + e
    0.22, H2 --> H + H
    0.00065, H2 --> H_p + H_m
    2, H2 --> H2_p + e
    3.8, C --> C_p + e
    5.7, O --> O_p + e
    6.5, CO --> CO_p + e
    2800, C --> C_p + e
    4000, CH --> C + H
    960, CH_p --> C_p + H
    2700, CH2 --> CH2_p + e
    2700, CH2 --> CH + H
    1300, C2 --> C + C
    2800, OH --> O + H
    5300, H2O --> OH + H
    4100, O2 --> O + O
    640, O2 --> O2_p + e
    525, CO --> C + O
end

u0 = [  
    :H_m => 1,
    :H2_p => 1,
    :H3_p => 1,
    :CH_p => 1,
    :CH2_p => 1, 
    :OH_p => 1,
    :H2O_p => 1,
    :H3O_p => 1,
    :CO_p => 1,
    :HOC_p => 1,
    :O_m => 1,
    :C_m => 1,
    :O2_p => 1, 
    :e => 1,
    :H_p => 1, 
    :H => 1,
    :H2 => 1, 
    :He => 1,
    :He_p => 1, 
    :C => 1,
    :C_p => 1, 
    :O => 1,
    :O_p => 1,
    :OH => 1,
    :H2O => 1, 
    :CO => 1,
    :C2 => 1,
    :O2 => 1,
    :HCO_p => 1, 
    :CH => 1,
    :CH2 => 1,
    :CH3_p => 1,


    
    :M => 1,
    ]


tspan = (0.0, 1.0)
#params = [:b => 1.0, :d => 0.2]
oprob = ODEProblem(Glover, u0, tspan; combinatoric_ratelaws = false)
sol = solve(oprob; adaptive = false, dt = 0.001, saveat = 0.1)
Plots.plot(sol)