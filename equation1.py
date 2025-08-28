from sympy import symbols, Eq, sin, solve, Sum, pi, Matrix, factorial, sqrt

import time

start_time = time.time()

N = 5
K = 10
b = 0.6
L = 2 * N * b
M = 60.21
EI = 6.38 * 10**6
V = 33.33
P0 = 140000
nsmal = 10

p_bar = 2 * P0 / (M * L)
omega_n = lambda n: (n**2 * pi**2 / L**2) * sqrt(EI / M)
omega_bar = lambda n: (n * pi * V) / L


a_ij = Matrix([[symbols(f'a_{i}_{2 * p - 1}') for p in range(1, K + 1)] for i in range(1, 2 * N)])


i = symbols('i')
print(a_ij)

def compute_C_n_2p_1(n, p, M, L, N, b):
    return Sum((2 / (M * L)) * a_ij[i - 1, p - 1] * sin(n * i * b * pi / L),
                  (i, 1, 2 * N - 1)).doit()

equations = []

for k in range(1, 1 + K):
    for l in range(1, 2 * N):
        result_l = 0
        for n in range(1, nsmal+1):
            for p in range(1, 1 + k):
                term3 = 0

                C_n_2p_1 = compute_C_n_2p_1(n, p, M, L, N, b)

                term1 = ((-1) ** (k + p) * C_n_2p_1 * (omega_n(n)/100) ** (2 * (k - p)) * factorial(2 * p - 1)) / factorial(2 * k + 1)
                term3 += term1

            term3 *= sin(n * pi * l / (2 * N))
            result_l += term3

            term2 = (-1) * (((-1) ** (k + 1)) / factorial(2 * k + 1)) * (
                        (p_bar * omega_bar(n) * ((omega_n(n)/100) ** (2 * k) - omega_bar(n) ** (2 * k))) / (
                        (omega_n(n)/100) ** 2 - omega_bar(n) ** 2))

            term2 *= sin(n * pi * l / (2 * N))
            result_l += term2

        equations.append(Eq(result_l, 0))

print(f"Number of equations: {len(equations)}")
print(f"Number of variables: {a_ij.shape[0] * a_ij.shape[1]}")


with open('../equations.txt', 'w') as file:
    for eq in equations:
        file.write(str(eq) + '\n')
solution = solve(equations, a_ij)

with open('../aij_solutions_with_N_K.txt', 'w') as file:
    output_file = '../aij_solutions_with_N_K.txt'
    with open(output_file, 'w') as file:
        file.write(f"N = {N}\n")
        file.write(f"K = {K}\n")
        file.write(f"b = {b}\n")
        file.write(f"L = {L}\n")
        file.write(f"V = {V}\n")
        file.write("Solutions for a_ij:\n")
        for sol in solution:
            file.write(f"{sol} = {solution[sol]}\n")

end_time = time.time()
execution_time = end_time - start_time

print(execution_time)
