## Packages
using OptimalControl
using Plots

### Initial conditions
#Initial
t0 = 0
# Final time
tf = 30
# Intial substrate concentration
s0 = 0.1
# Initial product concentration
p0 = 0.001
# Initial biomass concentration
r0 = 0.1
# Initial volume
V0 = 0.003


## Defintion of the optimal control problem
@def ocp begin

    t ∈ [ t0, tf ], time
    φ ∈ R⁴, state 
    u ∈ R, control

    s = φ₁
    p = φ₂ 
    r = φ₃ 
    V = φ₄ 

    φ(t0) == [ s0, p0, r0, V0 ]

    0 ≤ u(t) ≤ 1
    
    s(t) ≥ 0
    p(t) ≥ 0
    0 ≤ r(t) ≤ 1
    V(t) ≥ 0

   φ̇(t) == F0(φ(t)) + u(t) * F1(φ(t))

    V(tf) → max

end

# Defintions of the constants. For the moment, we use the same as in the paper, cf. [7, 27] line 202
const kᵣ = 1.1
const kₘ = 1.2
const Kᵣ = 1.3
const Kₘ = 1.4

# Michaelis-Menten kinetics functions
wᵣ(p) = kᵣ * p / (Kᵣ + p)
wₘ(s) = kₘ * s / (Kₘ + s)

# Defintion of the vector fields
F0 = VectorField( φ -> begin
    s, p, r, V = φ
    return [ -wₘ(s) * (1 - r) * V
              wₘ(s) * (1 - r) - wᵣ(p) * r * (p + 1)
             -wᵣ(p) * r^2
              wᵣ(p) * r * V ]
end )

F1 = VectorField( φ -> begin
    s, p, r, V = φ
    return [ 0, 0, wᵣ(p) * r, 0 ]
end )

# Solving the OCP
direct_sol = solve(ocp, grid_size=400)
# Plotting the solution
plt = plot(direct_sol, size=(600, 600))
savefig(plt, "figures/finiteTimeHorizonOptimalControlProblem.pdf")
