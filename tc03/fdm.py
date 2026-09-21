import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve


# ============================================================
# Parâmetros do problema
# ============================================================

Omega = [0, 1]

# Número de pontos internos
n = 6

# Espaçamento da malha
h = 1 / (n + 1)

print(f"h = {h}")

# Pontos internos
xj = np.linspace(Omega[0] + h, Omega[1] - h, n)


# ============================================================
# Solução manufaturada e função fonte
# ============================================================

def ue(x):
    """
    Solução analítica manufaturada.
    """
    return (x**2) * np.exp(x) + np.cos(3 * x)


def fx(x):
    """
    Função fonte f(x) = u_e''(x).
    """
    return ((x**2 + 4*x + 2) * np.exp(x)
            - 9 * np.cos(3 * x))


# Malha refinada apenas para plot da solução analítica
xe = np.linspace(Omega[0], Omega[1], 300)


# ============================================================
# Condições de contorno
# ============================================================

u0 = ue(Omega[0])
u1 = ue(Omega[1])

print("\nCondições de contorno:")
print(f"u(0) = {u0}")
print(f"u(1) = {u1}")


# ============================================================
# Gráfico da função fonte
# ============================================================

plt.figure(figsize=(8, 4))

plt.plot(xe, fx(xe), label="f(x)")
plt.scatter(xj, fx(xj), label="Pontos internos")

plt.title("Função fonte")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.grid(True)
plt.legend()

plt.savefig(
    "funcao_fonte_fdm.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# Construção da matriz do método de diferenças finitas
# ============================================================

main_diag = 2 * np.ones(n)
upper_diag = -1 * np.ones(n - 1)
lower_diag = -1 * np.ones(n - 1)

A = diags(
    [lower_diag, main_diag, upper_diag],
    offsets=[-1, 0, 1],
    format="csr"
)

print("\nSparse Matrix A:")
print(A)


# ============================================================
# Vetor do lado direito
# ============================================================

# A matriz foi construída na forma:
#
# 2*u_i - u_(i-1) - u_(i+1) = -h²*f(x_i)
#
# portanto:
b = -h**2 * fx(xj)

# Incorporação das condições de contorno
b[0] += u0
b[-1] += u1

print("\nVector b:")
print(b)


# ============================================================
# Visualização da estrutura da matriz
# ============================================================

plt.figure(figsize=(5, 5))
plt.spy(A)

plt.title("Estrutura da matriz do FDM")

plt.savefig(
    "matriz_fdm.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# Solução do sistema linear
# ============================================================

vj = spsolve(A, b)


# ============================================================
# Cálculo do RMSE
# ============================================================

RMSE = np.sqrt(
    np.mean(
        (vj - ue(xj))**2
    )
)

print("\nErro quadrático médio (RMSE) - FDM:")
print(RMSE)


# ============================================================
# Valores numéricos
# ============================================================

print("\nPontos internos x_j:")
print(xj)

print("\nSolução FDM nos pontos internos:")
print(vj)

print("\nSolução exata nos pontos internos:")
print(ue(xj))


# ============================================================
# Gráfico da solução FDM
# ============================================================

plt.figure(figsize=(8, 4))

# Solução analítica
plt.plot(
    xe,
    ue(xe),
    label="Solução analítica"
)

# Solução FDM
plt.scatter(
    xj,
    vj,
    label="FDM"
)

# Condições de contorno
plt.scatter(
    [Omega[0], Omega[1]],
    [u0, u1],
    label="Condições de contorno"
)

plt.title(
    f"Solução pelo Método das Diferenças Finitas - n = {n}"
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
    "solucao_fdm.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


print("\nFiguras geradas:")
print("- funcao_fonte_fdm.png")
print("- matriz_fdm.png")
print("- solucao_fdm.png")