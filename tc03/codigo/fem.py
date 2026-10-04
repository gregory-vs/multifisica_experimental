import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


# ============================================================
# Parâmetros do problema
# ============================================================

Omega = [0.0, 1.0]

# Número de pontos internos.
# Mesma convenção utilizada no fdm.py.
n = 6

# Número total de nós e de elementos
n_nodes = n + 2
n_elements = n_nodes - 1

# Coordenadas dos nós
x_nodes = np.linspace(Omega[0], Omega[1], n_nodes)

# Espaçamento uniforme
h = (Omega[1] - Omega[0]) / n_elements

print(f"h = {h}")
print(f"Número de pontos internos = {n}")
print(f"Número total de nós = {n_nodes}")
print(f"Número de elementos = {n_elements}")

print("\nNós usados na discretização:")
print(x_nodes)


# ============================================================
# Solução manufaturada e função fonte
# ============================================================

def ue(x):
    """
    Solução analítica manufaturada:
    u_e(x) = x² e^x + cos(3x)
    """
    return (x**2) * np.exp(x) + np.cos(3 * x)


def f(x):
    """
    Função fonte:
    f(x) = u_e''(x)
    """
    return (
        (x**2 + 4*x + 2) * np.exp(x)
        - 9 * np.cos(3 * x)
    )


# ============================================================
# Condições de contorno
# ============================================================

u0 = ue(Omega[0])
u1 = ue(Omega[1])

print("\nCondições de contorno:")
print(f"u(0) = {u0}")
print(f"u(1) = {u1}")


# ============================================================
# Funções de forma locais
# ============================================================

def N1(x, xa, xb):
    """
    Função de forma associada ao nó esquerdo do elemento.
    """
    return (xb - x) / (xb - xa)


def N2(x, xa, xb):
    """
    Função de forma associada ao nó direito do elemento.
    """
    return (x - xa) / (xb - xa)


# ============================================================
# Inicialização da matriz global e vetor global
# ============================================================

K = np.zeros((n_nodes, n_nodes))
F = np.zeros(n_nodes)


# ============================================================
# Montagem elemento por elemento
# ============================================================

for e in range(n_elements):

    xa = x_nodes[e]
    xb = x_nodes[e + 1]

    he = xb - xa

    # --------------------------------------------------------
    # Matriz de rigidez local
    #
    # K_e = integral (dN/dx)^T (dN/dx) dx
    # --------------------------------------------------------

    Ke = (1.0 / he) * np.array([
        [1.0, -1.0],
        [-1.0, 1.0]
    ])

    # --------------------------------------------------------
    # Vetor de forças local
    #
    # Como:
    #
    #     u'' = f
    #
    # a forma fraca é:
    #
    #     integral w' u' dx = - integral w f dx
    #
    # Portanto o vetor possui sinal negativo.
    # --------------------------------------------------------

    Fe = np.zeros(2)

    Fe[0] = -quad(
        lambda x: N1(x, xa, xb) * f(x),
        xa,
        xb
    )[0]

    Fe[1] = -quad(
        lambda x: N2(x, xa, xb) * f(x),
        xa,
        xb
    )[0]

    # Nós globais correspondentes ao elemento
    nodes = [e, e + 1]

    # --------------------------------------------------------
    # Montagem global
    # --------------------------------------------------------

    for a in range(2):

        A = nodes[a]

        F[A] += Fe[a]

        for b in range(2):

            B = nodes[b]

            K[A, B] += Ke[a, b]


print("\nMatriz global K antes das condições de contorno:")
print(K)

print("\nVetor global F antes das condições de contorno:")
print(F)


# ============================================================
# Aplicação das condições de contorno de Dirichlet
# ============================================================

# Antes de alterar a matriz, transferimos as contribuições
# dos valores prescritos para o vetor do lado direito.

F = F - K[:, 0] * u0
F = F - K[:, -1] * u1


# Primeira condição: u(0) = u0
K[0, :] = 0.0
K[:, 0] = 0.0
K[0, 0] = 1.0

F[0] = u0


# Segunda condição: u(1) = u1
K[-1, :] = 0.0
K[:, -1] = 0.0
K[-1, -1] = 1.0

F[-1] = u1


print("\nMatriz global K após as condições de contorno:")
print(K)

print("\nVetor global F após as condições de contorno:")
print(F)


# ============================================================
# Solução do sistema linear
# ============================================================

u_fem_nodes = np.linalg.solve(K, F)

print("\nSolução FEM nos nós:")
print(u_fem_nodes)

print("\nSolução exata nos nós:")
print(ue(x_nodes))


# ============================================================
# Cálculo do RMSE
# ============================================================

# Para manter a comparação coerente com o FDM,
# o RMSE é calculado nos mesmos 6 pontos internos.

u_fem_internal = u_fem_nodes[1:-1]
x_internal = x_nodes[1:-1]

RMSE = np.sqrt(
    np.mean(
        (u_fem_internal - ue(x_internal))**2
    )
)

print("\nErro quadrático médio (RMSE) - FEM:")
print(RMSE)


# ============================================================
# Interpolação da solução FEM
# ============================================================

x_plot = np.linspace(Omega[0], Omega[1], 300)
u_fem_plot = np.zeros_like(x_plot)

for i, xp in enumerate(x_plot):

    # Tratamento do último ponto do domínio
    if np.isclose(xp, Omega[1]):
        u_fem_plot[i] = u_fem_nodes[-1]
        continue

    # Identifica o elemento que contém xp
    e = np.searchsorted(x_nodes, xp, side="right") - 1

    e = max(0, min(e, n_elements - 1))

    xa = x_nodes[e]
    xb = x_nodes[e + 1]

    u_fem_plot[i] = (
        N1(xp, xa, xb) * u_fem_nodes[e]
        +
        N2(xp, xa, xb) * u_fem_nodes[e + 1]
    )


# ============================================================
# Gráfico da solução FEM
# ============================================================

plt.figure(figsize=(8, 4))

plt.plot(
    x_plot,
    ue(x_plot),
    label="Solução analítica"
)

plt.plot(
    x_plot,
    u_fem_plot,
    label="FEM"
)

plt.scatter(
    x_nodes,
    u_fem_nodes,
    label="Nós FEM"
)

plt.title(
    f"Solução pelo Método dos Elementos Finitos - n = {n}"
)

plt.xlabel("x")
plt.ylabel("u(x)")

plt.grid(True)
plt.legend()

plt.text(
    0.05,
    0.90,
    f"RMSE = {RMSE:.6e}",
    transform=plt.gca().transAxes,
    bbox=dict(
        boxstyle="round",
        facecolor="white",
        edgecolor="black"
    )
)

plt.savefig(
    "solucao_fem.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# Visualização da matriz de rigidez
# ============================================================

plt.figure(figsize=(5, 5))

plt.spy(K)

plt.title("Estrutura da matriz do FEM")

plt.savefig(
    "matriz_fem.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


print("\nFiguras geradas:")
print("- solucao_fem.png")
print("- matriz_fem.png")