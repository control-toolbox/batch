using SymEngine

function jacobian(f, vars)
    return [diff(expr, var) for expr in f, var in vars]
end

@vars s, p, r, v, k_r, k_m, K_r, K_m

F0 = [-(k_m * s / (K_m + s)) * (1-r)*v, (k_m * s / (K_m + s))*(1-r) - (k_r * p / (K_r + p)) *r * (p+1), -(k_r * p / (K_r + p))*r^2, (k_r * p / (K_r + p))*r*v]
F1 = [0, 0, (k_r * p / (K_r + p))*r, 0]

F0_prime = jacobian(F0, [s, p, r, v])
F1_prime = jacobian(F1, [s, p, r, v])

F01 = F1_prime * F0 - F0_prime * F1

F01_prime = jacobian(F01, [s, p, r, v])

F101 = F01_prime * F1 - F1_prime * F01

# Test avec valeurs :
s = 1
p = 1
r = 1
v = 1
k_r = 1
k_m = 1
K_r = 1
K_m = 1

F0test = [-(k_m * s / (K_m + s)) * (1-r)*v, (k_m * s / (K_m + s))*(1-r) - (k_r * p / (K_r + p)) *r * (p+1), -(k_r * p / (K_r + p))*r^2, (k_r * p / (K_r + p))*r*v]
F1test = [0, 0, (k_r * p / (K_r + p))*r, 0]

sF01 = string(F01)
sF01 = string(F101)

F01test = [-p*r*s*v*k_m*k_r/((K_m + s)*(K_r + p)), -p*r*k_r*(-s*k_m/(K_m + s) - p*(1 + p)*k_r/(K_r + p))/(K_r + p), (r*k_r/(K_r + p) - p*r*k_r/(K_r + p)^2)*(s*(1 - r)*k_m/(K_m + s) - p*r*(1 + p)*k_r/(K_r + p)) + p^2*r^2*k_r^2/(K_r + p)^2, -p^2*r*v*k_r^2/(K_r + p)^2]

F101test = [-p^2*r*s*v*k_m*k_r^2/((K_m + s)*(K_r + p)^2), -p^2*r*k_r^2*(-s*k_m/(K_m + s) - p*(1 + p)*k_r/(K_r + p))/(K_r + p)^2, p*r*k_r*((k_r/(K_r + p) - p*k_r/(K_r + p)^2)*(s*(1 - r)*k_m/(K_m + s) - p*r*(1 + p)*k_r/(K_r + p)) + (r*k_r/(K_r + p) - p*r*k_r/(K_r + p)^2)*(-s*k_m/(K_m + s) - p*(1 + p)*k_r/(K_r + p)) + 2*p^2*r*k_r^2/(K_r + p)^2)/(K_r + p) - (p*k_r*((r*k_r/(K_r + p) - p*r*k_r/(K_r + p)^2)*(s*(1 - r)*k_m/(K_m + s) - p*r*(1 + p)*k_r/(K_r + p)) + p^2*r^2*k_r^2/(K_r + p)^2)/(K_r + p) - p*r*k_r*(r*k_r/(K_r + p) - p*r*k_r/(K_r + p)^2)*(-s*k_m/(K_m + s) - p*(1 + p)*k_r/(K_r + p))/(K_r + p)), -p^3*r*v*k_r^3/(K_r + p)^3]

#------------------- Test Symbolics -------------------------#

using Symbolics

@variables s, p, r, v, k_r, k_m, K_r, K_m

F0_sym = [-(k_m * s / (K_m + s)) * (1-r)*v, (k_m * s / (K_m + s))*(1-r) - (k_r * p / (K_r + p)) *r * (p+1), -(k_r * p / (K_r + p))*r^2, (k_r * p / (K_r + p))*r*v]
F1_sym = [0, 0, (k_r * p / (K_r + p))*r, 0]

F0_prime_sym = Symbolics.jacobian(F0_sym, [s, p, r, v])
F1_prime_sym = Symbolics.jacobian(F1_sym, [s, p, r, v])

F01_sym = F1_prime_sym * F0_sym - F0_prime_sym * F1_sym

F01_prime_sym = Symbolics.jacobian(F01_sym, [s, p, r, v])

F101_sym = F01_prime_sym * F1_sym - F1_prime_sym * F01_sym

#Test avec valeurs :

s = 1
p = 1
r = 1
v = 1
k_r = 1
k_m = 1
K_r = 1
K_m = 1

F0test_sym = [-(k_m * s / (K_m + s)) * (1-r)*v, (k_m * s / (K_m + s))*(1-r) - (k_r * p / (K_r + p)) *r * (p+1), -(k_r * p / (K_r + p))*r^2, (k_r * p / (K_r + p))*r*v]
F1test_sym = [0, 0, (k_r * p / (K_r + p))*r, 0]

sF01_sym = string(F01_sym)
sF01_sym = string(F101_sym)

F01test_sym = [-((k_m*k_r*p*r*s*v) / ((K_m + s)*(K_r + p))), -((k_r*p*r*((-k_m*s) / (K_m + s) + (-k_r*p*(1 + p)) / (K_r + p))) / (K_r + p)), (-k_r*p*(r^2)*((k_r*p) / (K_r + p))) / (K_r + p) + (2(k_r^2)*(p^2)*(r^2)) / ((K_r + p)^2) + ((k_r*r) / (K_r + p) + (-k_r*p*r) / ((K_r + p)^2))*((k_m*(1 - r)*s) / (K_m + s) + (-k_r*p*(1 + p)*r) / (K_r + p)), -(((k_r^2)*(p^2)*r*v) / ((K_r + p)^2))]
F101test_sym = [(-k_m*(k_r^2)*(p^2)*r*s*v) / ((K_m + s)*((K_r + p)^2)), (-(k_r^2)*(p^2)*r*((-k_m*s) / (K_m + s) + (-k_r*p*(1 + p)) / (K_r + p))) / ((K_r + p)^2), (-k_r*p*((-k_r*p*(r^2)*((k_r*p) / (K_r + p))) / (K_r + p) + (2(k_r^2)*(p^2)*(r^2)) / ((K_r + p)^2) + ((k_r*r) / (K_r + p) + (-k_r*p*r) / ((K_r + p)^2))*((k_m*(1 - r)*s) / (K_m + s) + (-k_r*p*(1 + p)*r) / (K_r + p)))) / (K_r + p) + (k_r*p*r*((-2(k_r^2)*(p^2)*r) / ((K_r + p)^2) + (4(k_r^2)*(p^2)*r) / ((K_r + p)^2) + ((-k_m*s) / (K_m + s) + (-k_r*p*(1 + p)) / (K_r + p))*((k_r*r) / (K_r + p) + (-k_r*p*r) / ((K_r + p)^2)) + ((k_m*(1 - r)*s) / (K_m + s) + (-k_r*p*(1 + p)*r) / (K_r + p))*((-k_r*p) / ((K_r + p)^2) + k_r / (K_r + p)))) / (K_r + p) + ((k_r*r) / (K_r + p) + (-k_r*p*r) / ((K_r + p)^2))*((k_r*p*r*((-k_m*s) / (K_m + s) + (-k_r*p*(1 + p)) / (K_r + p))) / (K_r + p)), (-(k_r^3)*(p^3)*r*v) / ((K_r + p)^3)]

#------------------- Comparaison -------------------------#
if F0test == F0test_sym
    println("F0 == F0_sym")
else
    println("F0 != F0_sym")
end

if F1test == F1test_sym
    println("F1 == F1_sym")
else
    println("F1 != F1_sym")
end

if F01test == F01test_sym
    println("F01test == F01test_sym")
else
    println("F01test != F01test_sym")
end

if F101test == F101test_sym
    println("F101test == F101test_sym")
else
    println("F101test != F101test_sym")
end