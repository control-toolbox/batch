# batch.jl
using OptimalControl

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

H0 = Lift(F0)
H1 = Lift(F1)

H01 = @Lie { H0, H1 }
H101 = @Lie { H1, H01 }

φ = [ 1, 1, 1, 1 ]
λ = [ 1, 1, 1, 1 ]

println("H101 = ", H101(φ, λ))

using LinearAlgebra

φ = rand(4) 
F01 = @Lie [ F0, F1 ]
F101 = @Lie [ F1, F01 ]
println("rank(F1(φ), F01(φ)) = ", rank([ F1(φ) F01(φ) ]))
println("rank(F1(φ), F01(φ), F101(φ) = ", rank([ F1(φ) F01(φ) F101(φ) ]))