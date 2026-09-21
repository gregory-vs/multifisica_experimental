import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# Parâmetros do problema
# ============================================================

L = 1.0
Va = 10.0
Vb = 0.0

# O enunciado pede n ímpar para que exista um nó exatamente
# na interface z = L/2. Aqui são utilizados 9 nós.
n_nodes = 9

# Malha uniforme
z = np.linspace(0.0, L, n_nodes)

# Interface entre os dielétricos
interface = L / 2.0

if not np.any(np.isclose(z, interface)):
    raise ValueError(
        "A malha escolhida não possui um nó em z = L/2. "
        "Escolha um número ímpar de nós."
    )

h = z[1] - z[0]

print("Parâmetros:")
print(f"L = {L}")
print(f"Va = {Va}")
print(f"Vb = {Vb}")
print(f"Número de nós = {n_nodes}")
print(f"h = {h}")
print(f"Interface = {interface}")

print("\nNós da malha:")
print(z)


# ============================================================
# Permissividade relativa
# ============================================================

def epsilon_r(position):
    """
    Permissividade relativa definida pelo enunciado:

        epsilon_r = 3, para z < L/2
        epsilon_r = 1, para z >= L/2
    """
    position = np.asarray(position)

    return np.where(
        position < interface,
        3.0,
        1.0
    )


# ============================================================
# Solução pelo Método dos Elementos Finitos (FEM)
# ============================================================

def solve_fem():
    """
    FEM 1D com elementos lineares.

    Para cada elemento:

        K_e = epsilon_e / h_e * [[1, -1],
                                 [-1, 1]]

    A permissividade de cada elemento é determinada
    pelo ponto médio do elemento.
    """

    K = np.zeros((n_nodes, n_nodes))
    F = np.zeros(n_nodes)

    # Montagem elemento por elemento
    for e in range(n_nodes - 1):

        z1 = z[e]
        z2 = z[e + 1]

        he = z2 - z1
        midpoint = 0.5 * (z1 + z2)

        eps_e = float(epsilon_r(midpoint))

        Ke = (eps_e / he) * np.array([
            [1.0, -1.0],
            [-1.0, 1.0]
        ])

        nodes = [e, e + 1]

        for a in range(2):
            for b in range(2):
                K[nodes[a], nodes[b]] += Ke[a, b]

    # --------------------------------------------------------
    # Condições de contorno de Dirichlet
    # --------------------------------------------------------

    F = F - K[:, 0] * Va
    F = F - K[:, -1] * Vb

    # phi(0) = Va
    K[0, :] = 0.0
    K[:, 0] = 0.0
    K[0, 0] = 1.0
    F[0] = Va

    # phi(L) = Vb
    K[-1, :] = 0.0
    K[:, -1] = 0.0
    K[-1, -1] = 1.0
    F[-1] = Vb

    phi = np.linalg.solve(K, F)

    return phi, K, F


# ============================================================
# Solução pelo Método das Diferenças Finitas (FDM)
# ============================================================

def solve_fdm():
    """
    Discretização conservativa de

        d/dz (epsilon_r dphi/dz) = 0

    usando a permissividade nas faces da malha:

        -eps_(i-1/2) phi_(i-1)
        +(eps_(i-1/2)+eps_(i+1/2)) phi_i
        -eps_(i+1/2) phi_(i+1) = 0

    Como a interface coincide com um nó, cada face pertence
    integralmente a apenas um dos dois materiais.
    """

    A = np.zeros((n_nodes, n_nodes))
    b = np.zeros(n_nodes)

    # Condição de contorno esquerda
    A[0, 0] = 1.0
    b[0] = Va

    # Pontos internos
    for i in range(1, n_nodes - 1):

        z_left_face = 0.5 * (z[i - 1] + z[i])
        z_right_face = 0.5 * (z[i] + z[i + 1])

        eps_left = float(epsilon_r(z_left_face))
        eps_right = float(epsilon_r(z_right_face))

        A[i, i - 1] = -eps_left
        A[i, i] = eps_left + eps_right
        A[i, i + 1] = -eps_right

    # Condição de contorno direita
    A[-1, -1] = 1.0
    b[-1] = Vb

    phi = np.linalg.solve(A, b)

    return phi, A, b


# ============================================================
# Cálculo das soluções
# ============================================================

phi_fem, K_fem, F_fem = solve_fem()
phi_fdm, A_fdm, b_fdm = solve_fdm()


# ============================================================
# Resultados numéricos
# ============================================================

interface_index = np.where(np.isclose(z, interface))[0][0]

print("\nPotencial FEM:")
print(phi_fem)

print("\nPotencial FDM:")
print(phi_fdm)

print("\nPotencial na interface z = L/2:")
print(f"FEM = {phi_fem[interface_index]:.10f}")
print(f"FDM = {phi_fdm[interface_index]:.10f}")

max_difference = np.max(np.abs(phi_fem - phi_fdm))

print("\nDiferença máxima entre FEM e FDM:")
print(f"{max_difference:.10e}")


# ============================================================
# Fluxo elétrico por elemento/intervalo
# ============================================================

flux_fem = np.zeros(n_nodes - 1)
flux_fdm = np.zeros(n_nodes - 1)

for e in range(n_nodes - 1):

    midpoint = 0.5 * (z[e] + z[e + 1])
    eps_e = float(epsilon_r(midpoint))

    flux_fem[e] = eps_e * (
        phi_fem[e + 1] - phi_fem[e]
    ) / h

    flux_fdm[e] = eps_e * (
        phi_fdm[e + 1] - phi_fdm[e]
    ) / h

print("\nFluxo epsilon_r * dphi/dz por elemento (FEM):")
print(flux_fem)

print("\nFluxo epsilon_r * dphi/dz por intervalo (FDM):")
print(flux_fdm)


# ============================================================
# Gráfico comparativo
# ============================================================

plt.figure(figsize=(8, 5))

plt.plot(
    z,
    phi_fem,
    marker="o",
    linewidth=2,
    label="FEM"
)

plt.plot(
    z,
    phi_fdm,
    marker="s",
    linestyle="--",
    linewidth=2,
    label="FDM"
)

plt.axvline(
    interface,
    linestyle=":",
    label="Interface dos dielétricos"
)

plt.xlabel("z")
plt.ylabel(r"$\phi(z)$")
plt.title("Potencial elétrico no capacitor com dois dielétricos")

plt.grid(True)
plt.legend()

plt.savefig(
    "comparacao_capacitor.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# Gráfico da permissividade relativa
# ============================================================

z_plot = np.linspace(0.0, L, 500)

plt.figure(figsize=(8, 4))

plt.step(
    z_plot,
    epsilon_r(z_plot),
    where="post"
)

plt.axvline(
    interface,
    linestyle=":"
)

plt.xlabel("z")
plt.ylabel(r"$\varepsilon_r(z)$")
plt.title("Distribuição da permissividade relativa")

plt.grid(True)

plt.savefig(
    "permissividade_capacitor.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


print("\nImagens geradas:")
print("- comparacao_capacitor.png")
print("- permissividade_capacitor.png")