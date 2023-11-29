using OptimalControl

t0 = 0      # initial time
tf = 30    # final time
s0 = 0.1
p0 = 0.001
r0 = 0.1
V0 = 0.003

@def ocp begin # definition of the optimal control problem

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
    0 ≤ u(t) ≤ 1

    ẋ(t) == F0(x(t)) + u(t) * F1(x(t))

    v(tf) → max

end;

# Dynamics
const kᵣ = 1.1
const kₘ = 1.2
const Kᵣ = 1.3
const Kₘ = 1.4


wᵣ(p) = kᵣ * p / (Kᵣ + p)
wₘ(s) = kₘ * s / (Kₘ + s)

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

direct_sol1 = solve(ocp, grid_size=100)

direct_sol2 = solve(ocp, grid_size=1000)

plt1 = plot(direct_sol1, size=(600, 600))
plt2 = plot(direct_sol2, size=(600, 600))
