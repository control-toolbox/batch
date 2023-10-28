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
    x ∈ R⁴, state 
    u ∈ R, control

    s = x₁
    p = x₂
    r = x₃
    v = x₄

    x(t0) == [ s0, p0, r0, V0 ]
    
    s(t) ≥ 0
    p(t) ≥ 0
    0 ≤ r(t) ≤ 1
    v(t) ≥ 0

    ẋ(t) == F0(x(t)) + u(t) * F1(x(t))

    v(tf) → max

end;

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
direct_sol = solve(ocp, grid_size=100)
# Plotting the solution
plt = plot(direct_sol, size=(600, 600))
savefig(plt, "figures/finiteTimeHorizonOptimalControlProblem.pdf")