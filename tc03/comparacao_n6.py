import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from scipy.integrate import quad


# ============================================================
# Problema
# ============================================================

Omega = [0.0, 1.0]
n = 6  # número de pontos internos

def ue(x):
    return (x**2) * np.exp(x) + np.cos(3 * x)

def f(x):
    return ((x**2 + 4*x + 2) * np.exp(x) - 9 * np.cos(3 * x))


# ============================================================
# Solução analítica
# ============================================================

x_exact = np.linspace(Omega[0], Omega[1], 400)
u_exact = ue(x_exact)

u0 = ue(Omega[0])
u1 = ue(Omega[1])


# ============================================================
# FDM
# ============================================================

h_fdm = 1 / (n + 1)
x_fdm = np.linspace(Omega[0] + h_fdm, Omega[1] - h_fdm, n)

main_diag = 2 * np.ones(n)
upper_diag = -1 * np.ones(n - 1)
lower_diag = -1 * np.ones(n - 1)

A = diags(
    [lower_diag, main_diag, upper_diag],
    offsets=[-1, 0, 1],
    format="csr"
)

b = -h_fdm**2 * f(x_fdm)
b[0] += u0
b[-1] += u1

u_fdm = spsolve(A, b)

rmse_fdm = np.sqrt(np.mean((u_fdm - ue(x_fdm))**2))


# ============================================================
# FEM
# Mesma discretização do FDM:
# n pontos internos -> n+2 nós -> n+1 elementos
# ============================================================

n_nodes = n + 2
n_elements = n_nodes - 1

x_nodes = np.linspace(Omega[0], Omega[1], n_nodes)

def N1(x, xa, xb):
    return (xb - x) / (xb - xa)

def N2(x, xa, xb):
    return (x - xa) / (xb - xa)

K = np.zeros((n_nodes, n_nodes))
F = np.zeros(n_nodes)

for e in range(n_elements):
    xa = x_nodes[e]
    xb = x_nodes[e + 1]
    he = xb - xa

    Ke = (1.0 / he) * np.array([
        [1.0, -1.0],
        [-1.0, 1.0]
    ])

    Fe = np.zeros(2)
    Fe[0] = -quad(lambda x: N1(x, xa, xb) * f(x), xa, xb)[0]
    Fe[1] = -quad(lambda x: N2(x, xa, xb) * f(x), xa, xb)[0]

    nodes = [e, e + 1]

    for a in range(2):
        A_idx = nodes[a]
        F[A_idx] += Fe[a]

        for b_idx in range(2):
            B_idx = nodes[b_idx]
            K[A_idx, B_idx] += Ke[a, b_idx]

# aplicar Dirichlet
F = F - K[:, 0] * u0
F = F - K[:, -1] * u1

K[0, :] = 0.0
K[:, 0] = 0.0
K[0, 0] = 1.0
F[0] = u0

K[-1, :] = 0.0
K[:, -1] = 0.0
K[-1, -1] = 1.0
F[-1] = u1

u_fem_nodes = np.linalg.solve(K, F)

# interpolação da solução FEM para plot contínuo
x_fem_plot = np.linspace(Omega[0], Omega[1], 400)
u_fem_plot = np.zeros_like(x_fem_plot)

for i, xp in enumerate(x_fem_plot):
    if np.isclose(xp, Omega[1]):
        u_fem_plot[i] = u_fem_nodes[-1]
        continue

    e = np.searchsorted(x_nodes, xp, side="right") - 1
    e = max(0, min(e, n_elements - 1))

    xa = x_nodes[e]
    xb = x_nodes[e + 1]

    u_fem_plot[i] = (
        N1(xp, xa, xb) * u_fem_nodes[e]
        + N2(xp, xa, xb) * u_fem_nodes[e + 1]
    )

x_fem_internal = x_nodes[1:-1]
u_fem_internal = u_fem_nodes[1:-1]

rmse_fem = np.sqrt(np.mean((u_fem_internal - ue(x_fem_internal))**2))


# ============================================================
# Gráfico combinado
# ============================================================

plt.figure(figsize=(10, 6))

# solução analítica
plt.plot(x_exact, u_exact, linewidth=2, label="Solução analítica")

# FEM
plt.plot(x_fem_plot, u_fem_plot, linewidth=2, label="FEM")
plt.scatter(x_nodes, u_fem_nodes, s=45, label="Nós FEM")

# FDM
plt.scatter(x_fdm, u_fdm, s=60, label="FDM")

plt.title("Comparação entre solução analítica, FEM e FDM - n = 6")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.grid(True)
plt.legend()

plt.text(
    0.03, 0.95,
    f"RMSE FDM = {rmse_fdm:.6e}\nRMSE FEM = {rmse_fem:.6e}",
    transform=plt.gca().transAxes,
    verticalalignment="top",
    bbox=dict(boxstyle="round", facecolor="white", edgecolor="black")
)

plt.savefig("comparacao_n6.png", dpi=300, bbox_inches="tight")
plt.close()


# ============================================================
# Saída no terminal
# ============================================================

print("RMSE - FDM:")
print(rmse_fdm)

print("\nRMSE - FEM:")
print(rmse_fem)

print("\nFigura gerada:")
print("- comparacao_n6.png")