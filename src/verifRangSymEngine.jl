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

#--------------------------------------------

mat_F1_F01 = [F1 F01]

# Sous matrice de mat_F1_F01
sousDet1_F1_F01 = [mat_F1_F01[1] mat_F1_F01[5] ; mat_F1_F01[2] mat_F1_F01[6]]

sousDet2_F1_F01 = [mat_F1_F01[1] mat_F1_F01[5] ; mat_F1_F01[3] mat_F1_F01[7]]

sousDet3_F1_F01 = [mat_F1_F01[1] mat_F1_F01[5] ; mat_F1_F01[4] mat_F1_F01[8]]

sousDet4_F1_F01 = [mat_F1_F01[2] mat_F1_F01[6] ; mat_F1_F01[3] mat_F1_F01[7]]

sousDet5_F1_F01 = [mat_F1_F01[2] mat_F1_F01[6] ; mat_F1_F01[4] mat_F1_F01[8]]

sousDet6_F1_F01 = [mat_F1_F01[3] mat_F1_F01[7] ; mat_F1_F01[4] mat_F1_F01[8]]

using LinearAlgebra

mineur1_F1_F01 = det(sousDet1_F1_F01)
mineur2_F1_F01 = det(sousDet2_F1_F01)
mineur3_F1_F01 = det(sousDet3_F1_F01)
mineur4_F1_F01 = det(sousDet4_F1_F01)
mineur5_F1_F01 = det(sousDet5_F1_F01)
mineur6_F1_F01 = det(sousDet6_F1_F01)

# Calcul du rang de la matrice F1 F01 et F101
mat_F1_F01_F101 = [F1 F01 F101]

# Sous matrice de mat2
sousDet1_F1_F01_F101 = [mat_F1_F01_F101[1] mat_F1_F01_F101[5] mat_F1_F01_F101[9] ; mat_F1_F01_F101[2] mat_F1_F01_F101[6] mat_F1_F01_F101[10] ; mat_F1_F01_F101[3] mat_F1_F01_F101[7] mat_F1_F01_F101[11]]

sousDet2_F1_F01_F101 = [mat_F1_F01_F101[1] mat_F1_F01_F101[5] mat_F1_F01_F101[9] ; mat_F1_F01_F101[2] mat_F1_F01_F101[6] mat_F1_F01_F101[10] ; mat_F1_F01_F101[4] mat_F1_F01_F101[8] mat_F1_F01_F101[12]]

sousDet3_F1_F01_F101 = [mat_F1_F01_F101[1] mat_F1_F01_F101[5] mat_F1_F01_F101[9] ; mat_F1_F01_F101[3] mat_F1_F01_F101[7] mat_F1_F01_F101[11] ; mat_F1_F01_F101[4] mat_F1_F01_F101[8] mat_F1_F01_F101[12]]

sousDet4_F1_F01_F101 = [mat_F1_F01_F101[2] mat_F1_F01_F101[6] mat_F1_F01_F101[10] ; mat_F1_F01_F101[3] mat_F1_F01_F101[7] mat_F1_F01_F101[11] ; mat_F1_F01_F101[4] mat_F1_F01_F101[8] mat_F1_F01_F101[12]]

using LinearAlgebra

mineur1_F1_F01_F101 = det(sousDet1_F1_F01_F101)
mineur2_F1_F01_F101 = det(sousDet2_F1_F01_F101)
mineur3_F1_F01_F101 = det(sousDet3_F1_F01_F101)
mineur4_F1_F01_F101 = det(sousDet4_F1_F01_F101)