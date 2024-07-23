using Symbolics
using Markdown

@variables s, p, r, v, k_r, k_m, K_r, K_m

wᵣ(p) = kᵣ * p / (Kᵣ + p)
wₘ(s) = kₘ * s / (Kₘ + s)

F0 = [-wₘ(s) * (1-r)*v, wₘ(s)*(1-r) - wᵣ(p) *r * (p+1), -wᵣ(p)*r^2, wᵣ(p)*r*v]
F1 = [0, 0, wᵣ(p)*r, 0]

F0_prime = Symbolics.jacobian(F0, [s, p, r, v])
F1_prime = Symbolics.jacobian(F1, [s, p, r, v])

F01 = F1_prime * F0 - F0_prime * F1

F01_prime = Symbolics.jacobian(F01, [s, p, r, v])

F101 = F01_prime * F1 - F1_prime * F01

print(F101)
print(F01)
